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
import re
import pandas as pd
import seaborn as sn
from sklearn.feature_selection import VarianceThreshold
from skimage.util import img_as_ubyte
from bioimg import read_image, threshold_img
from bioimg import plot_channels, combine_channels
from bioimg import SegfreeProfiler

def minmax_scale(a):
    return (a - a.min()) / (a.max() - a.min())

def normalize_channels(img):
    return np.stack([minmax_scale(img[:,:,i]) for i in range(img.shape[-1]) ], axis=-1)

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)
    
    path = 'data/cytodata/datasets/CDRPBIO-BBBC036-Bray/images/CDRPBIO-BBBC036-Bray'
    plate = str(sys.argv[1])
    print("Processing plate: " + plate)

    platedir = os.path.join(path, plate)
    fnames = [f for f in os.listdir(platedir) if '.tif' in f]
    all_wells = list(set([re.search('_[a-z][0-9]+_',f).group(0) for f in fnames]))
    all_wells.sort()
    assert(len(all_wells) == 384)

    well_id = int(sys.argv[2]) - 1
    well = all_wells[well_id]

    well_imgs = [f for f in fnames if well in f]
    assert(len(well_imgs) == 30)

    imglist = []
    for i in range(6):
        fview = '_s'+ str(i+1) + '_'
        imgs = [read_image(fname=os.path.join(platedir, f),
                           verbose=False) for f in well_imgs if fview in f]
        imglist.append(np.stack(imgs, axis=-1))

    outdir = os.path.join('data/segfree/cellpainting/', plate)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    hoechst = [img[:,:,0] for img in imglist]
    # threshold and convert to ubyte
    hoechst = [threshold_img(img, method='otsu') for img in hoechst]
    hoechst = [img_as_ubyte(img) for img in hoechst]

    segf = pickle.load(open("segf_nuclei.pkl", "rb"))
    nucl_prof = segf.transform(hoechst)
    

    # and segmentation-free profiles for all channels
    imgs_norm = [normalize_channels(img) for img in imglist ]
    segf = pickle.load(open("segf_cells.pkl", "rb"))
    cell_prof = segf.transform(imgs_norm)
    

    # combine nucl_prof and cell_prof
    nucl_prof.columns = ['-'.join(['nuclei', col]) for col in nucl_prof.columns.values]
    cell_prof.columns = ['-'.join(['cell', col])
                               for col in cell_prof.columns.values]

    segf_prof = pd.concat([nucl_prof, cell_prof], axis=1)
    df_out = pd.DataFrame(segf_prof.agg('median')).T
    df_out.to_csv(os.path.join(outdir, well.replace('_', '')+'.csv' ), index=False)

    javabridge.kill_vm()
