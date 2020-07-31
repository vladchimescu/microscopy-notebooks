#!/usr/bin/env python3
'''
Script for fitting unsupervised segmentation-free model
on CellPainting data
'''
import matplotlib as mpl
mpl.use('Agg')
import javabridge
import bioformats as bf
import skimage
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pickle
import pandas as pd
import seaborn as sn
from sklearn.feature_selection import VarianceThreshold
from skimage.util import img_as_ubyte
from bioimg import read_image, threshold_img
from bioimg import plot_channels, combine_channels
from bioimg import SegfreeProfiler

drugs = ['vinblastine sulfate',
         'dorsomorphin',
         'amthamine',
         'carboxin',
         'resveratrol',
         'KN-93',
         '5-iodotubercidin',
         'harmol', 'scoulerine',
         'LY-294002',
         'VEGF receptor 2 kinase inhibitor IV',
         'A-23187',
         'capsazepine', 'crustecdysone',
         'butein', 'SB 218078',
         'hinokitiol', 'SRC kinase inhibitor II']

col_params={'colors': ['blue', 
                       'yellow',
                       'green', 
                       'white', 
                       'red'],
            'blend': [2, 0.5, 1, 0.5, 1],
            'gamma': [0.5, 0.8,0.6,0.5,0.7]}

def minmax_scale(a):
    return (a - a.min()) / (a.max() - a.min())

def normalize_channels(img):
    return np.stack([minmax_scale(img[:,:,i]) for i in range(img.shape[-1]) ], axis=-1)

def load_drug_images(path, chemannot, drug, which=0):
    '''Retrun list of CellPainting image file names for the selected drug
    '''
    cmpd_annot = chemannot[chemannot['CPD_NAME'] == drug].reset_index(drop=True)
    plate = cmpd_annot['plateID'][which]
    well = cmpd_annot['well_position'][which]
    
    imgpath = os.path.join(path, str(plate))
    imgfiles = os.listdir(imgpath)
    cmpd_imgs = [f for f in imgfiles if "_" + well + "_" in f]
    
    ctrl_annot = chemannot.loc[np.logical_and(chemannot['broad_sample'] == 'DMSO',
                               chemannot['plateID'] == plate),:].reset_index(drop=True)
    well = ctrl_annot['well_position'][which]
    ctrl_imgs = [f for f in imgfiles if "_" + well + "_" in f]
    
    return imgpath, cmpd_imgs, ctrl_imgs

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)
    
    path = 'data/cytodata/datasets/CDRPBIO-BBBC036-Bray/images/CDRPBIO-BBBC036-Bray'
    chemannot = pd.read_csv("data/Bray-metadata/plate_annot.txt", sep='\t')

    # remove noisy wells
    noisywells = pd.read_csv('data/Bray-noisywells.csv', index_col=0)
    noisywells = noisywells[noisywells['count'] > 3].reset_index(drop=True)

    chemannot = pd.merge(left=chemannot,
                         right=noisywells[['plateID', 'well_position']],
                         how='outer', indicator=True)
    chemannot = (chemannot[chemannot['_merge'] == 'left_only'].
                 drop(columns='_merge').reset_index(drop=True))

    outdir = 'figures/cellpainting/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    imglist = []
    titles = []
    for d in drugs:
        for i in range(3):
            imgpath, cmpd_imgs, ctrl_imgs = load_drug_images(path=path, chemannot=chemannot,
                                                        drug=d,
                                                        which=i)
            fview = 's3'
            imgs = [read_image(fname=os.path.join(imgpath, f), verbose=False) for f in cmpd_imgs if fview in f]
            dmso_imgs = [read_image(fname=os.path.join(imgpath, f), verbose=False) for f in ctrl_imgs if fview in f]

            imglist.append(np.stack(imgs, axis=-1))
            titles.append(d)
            imglist.append(np.stack(dmso_imgs, axis=-1))
            titles.append('DMSO')

    hoechst = [img[:,:,0] for img in imglist]
    # threshold and convert to ubyte
    hoechst = [threshold_img(img, method='otsu') for img in hoechst]
    hoechst = [img_as_ubyte(img) for img in hoechst]

    segf = SegfreeProfiler(tile_size=(10,10),
                       n_block_types=10,
                       n_supblock_types=10)
    nucl_prof = segf.fit_transform(hoechst)
    pickle.dump(segf, open("segf_nuclei.pkl", "wb"))
    # segf = pickle.load(open("segf_nuclei.pkl", "rb"))

    # and segmentation-free profiles for all channels combined
    segf = SegfreeProfiler(tile_size=(10,10),
                       n_block_types=20,
                       n_supblock_types=20)
    imgs_norm = [normalize_channels(img) for img in imglist ]
    cell_prof = segf.fit_transform(imgs_norm)
    pickle.dump(segf, open("segf_cells.pkl", "wb"))
    # segf = pickle.load(open("segf_cells.pkl", "rb"))

    plt.plot(np.cumsum(segf.pca.explained_variance_ratio_), linewidth=3)
    sn.despine()
    plt.axhline(y=1, color='black', linestyle=':')
    plt.xlabel('Number of principal components')
    plt.ylabel('Cumulative explained variance')
    plt.savefig(os.path.join(outdir, 'cumvariance-PCA-segfree.pdf'))

    nchan = 5
    eigentiles = segf.pca.components_.reshape((segf.n_components, *segf.tile_size, nchan))
    plot_channels([combine_channels([eigentiles[j,:,:,i] for i in range(5)], **col_params) for j in range(50)],
              nrow=5, ncol=10, scale_x=2, scale_y=2)
    plt.savefig(os.path.join(outdir, 'eigentiles.pdf'),
                bbox_inches='tight')

    # combine nucl_prof and cell_prof
    nucl_prof.columns = ['-'.join(['nuclei', col]) for col in nucl_prof.columns.values]
    cell_prof.columns = ['-'.join(['cell', col])
                               for col in cell_prof.columns.values]

    segf_prof = pd.concat([nucl_prof, cell_prof], axis=1)
    segf_prof.index = titles
    sel = VarianceThreshold(threshold=1e-4).fit(segf_prof)
    sn.clustermap(segf_prof.loc[:,sel.get_support()])
    plt.savefig(os.path.join(outdir, 'fit-heatmap.pdf'))

    javabridge.kill_vm()
