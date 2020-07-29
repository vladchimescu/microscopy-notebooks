#!/usr/bin/env python3
'''
Script for fitting unsupervised segmentation-free model
on BiTE data
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
import re
import random
import pickle
import pandas as pd
import seaborn as sn
from sklearn.feature_selection import VarianceThreshold
from skimage.util import img_as_ubyte
from bioimg import load_image_series, threshold_img
from bioimg import plot_channels, combine_channels
from bioimg import SegfreeProfiler

col_params = dict(colors=['blue', 'red', 
                          'yellow', 'green'],
                             blend = [2, 0.8, 0.8, 0.8],
                             gamma = [0.6, 0.8, 0.8, 0.8])

plates = ['BiTE-Screen_Tag3_2019-01-15__2019-01-15T11_17_12-Measurement 1',
          'BiTE-Screen_Tag4_2019-01-17__2019-01-17T10_50_57-Measurement 1',
          'BiTE-Screen_Tag8_2019-02-14__2019-02-14T10_23_55-Measurement 1',
          'BiTE-Screen_Tag15_2019-04-30_Plate1__2019-04-30T10_18_30-Measurement 1',
          'BiTE-Screen_Tag15_2019-04-30_Plate2__2019-04-30T11_17_13-Measurement 1']

def minmax_scale(a):
    return (a - a.min()) / (a.max() - a.min())

def normalize_channels(img):
    return np.stack([minmax_scale(img[:,:,i]) for i in range(img.shape[-1]) ], axis=-1)

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)

    # edit this path
    # path = 'data/Tobias/newscreen/'
    path = sys.argv[1]
    annot = pd.read_csv('data/BiTE/bite_annot.csv')

    imglist = []
    titles = []
    random.seed(2707)
    for plate in plates:
        pl = re.search('(.+)__(.+)', plate).group(1)
        pl = pl.replace('Screen_', '')
        if 'Tag15' not in pl:
            pl = re.search('(.+)_(.+)', pl).group(1)
        plate_annot = annot[annot['plate']==pl]
        platedir = os.path.join(path, plate, 'Images')
        fnames = [f for f in os.listdir(platedir) if '.tiff' in f]
        all_wells = list(set([re.search('r[0-9]+c[0-9]+',f).group(0) for f in fnames]))
        all_wells.sort()
        well_ind = random.sample(range(96), 10)
        for w in well_ind:
            well = all_wells[w]
            if plate_annot[plate_annot['well']==well].size:
                titles.append(plate_annot[plate_annot['well']==well]['comb'].values[0])
            else:
                titles.append('NA')
            well_imgs = [f for f in fnames if well in f and 'ch2' not in f]
            fview = 'f0'+str(random.choice(range(4)) + 1) + 'p'
            wellpos = [f for f in well_imgs if fview in f]
            imgseries = load_image_series(path=platedir, imgfiles=wellpos)
            imgseries = np.stack(imgseries)
            imgseries = imgseries.reshape((8, 4, 2160,2160))
            mipseries = np.amax(imgseries, axis=0)
            mipseries = mipseries.swapaxes(0,-1)
            imglist.append(mipseries)
    
    outdir = 'figures/segfree/bite/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    hoechst = [img[:,:,0] for img in imglist]
    # threshold and convert to ubyte
    hoechst = [threshold_img(img, method='otsu') for img in hoechst]
    hoechst = [img_as_ubyte(img) for img in hoechst]

    segf = SegfreeProfiler(tile_size=(10,10),
                       n_block_types=20,
                       n_supblock_types=20)
    nucl_prof = segf.fit_transform(hoechst)
    pickle.dump(segf, open("segf_bite_nuclei.pkl", "wb"))
    # segf = pickle.load(open("segf_bite_nuclei.pkl", "rb"))

    # and segmentation-free profiles for all channels combined
    segf = SegfreeProfiler(tile_size=(10,10),
                       n_block_types=20,
                       n_supblock_types=20)
    imgs_norm = [normalize_channels(img) for img in imglist ]
    cell_prof = segf.fit_transform(imgs_norm)
    pickle.dump(segf, open("segf_bite_cells.pkl", "wb"))
    # segf = pickle.load(open("segf_bite_cells.pkl", "rb"))

    plt.plot(np.cumsum(segf.pca.explained_variance_ratio_), linewidth=3)
    sn.despine()
    plt.axhline(y=1, color='black', linestyle=':')
    plt.xlabel('Number of principal components')
    plt.ylabel('Cumulative explained variance')
    plt.savefig(os.path.join(outdir, 'cumvariance-PCA-segfree.pdf'))

    nchan = 4
    eigentiles = segf.pca.components_.reshape((segf.n_components,
                                               *segf.tile_size, nchan))
    plot_channels([combine_channels([eigentiles[j,:,:,i] for i in range(nchan)], **col_params) for j in range(50)],
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
