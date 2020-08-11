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
import re
import random
import pickle
import pandas as pd
import seaborn as sn
from sklearn.feature_selection import VarianceThreshold
from skimage.util import img_as_ubyte
from bioimg import load_image_series
from bioimg import plot_channels, combine_channels
from bioimg import SegfreeProfiler

def get_all_wells(path):
    fnames = [f for f in os.listdir(path) if 'tif' in f]
    fnames.sort()
    all_wells = list(set([re.search('--(W[0-9]+)--(.+)', f).group(1) for f in fnames]))
    return fnames, all_wells

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)
    
    datadir = 'data/organoids-sylwia/CL34/'

    paths = [os.path.join(datadir, 'CL34A_4x_2p1mm_100umsteps_day0_001/data'),
             os.path.join(datadir, 'CL34B_4x_2p1mm_100umsteps_day0_001/data'),
             os.path.join(datadir, 'CL34A_4x_2p1mm_100umsteps_day3_001/data'),
             os.path.join(datadir, 'CL34B_4x_2p1mm_100umsteps_day3_001/data')]
   
    outdir = 'figures/segfree/organoids'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    random.seed(1108)

    imglist = []
    titles = []
    for path in paths:
        fnames, all_wells = get_all_wells(path=path)
        # randomly sample 10 wells
        sel_wells = random.sample(all_wells, k=5)
        for w in sel_wells:
            print("Processing well: %s" % w)
            well_files = [f for f in fnames if w in f]
            imgstack = load_image_series(path=path, imgfiles=[w for w in well_files if 'P00001' in w])
            imgstack = imgstack.swapaxes(0,-1)
            imglist.append(imgstack)
            titles.append(re.search('(.+)(--W[0-9]+)', well_files[0]).group(1))

    # 3D profiling
    segf = SegfreeProfiler(tile_size=(20,20),
                       n_block_types=20,
                       n_supblock_types=50)
    org_prof = segf.fit_transform(imglist)
    pickle.dump(segf, open("segf_organoids.pkl", "wb"))

    plt.plot(np.cumsum(segf.pca.explained_variance_ratio_), linewidth=3)
    sn.despine()
    plt.axhline(y=1, color='black', linestyle=':')
    plt.xlabel('Number of principal components')
    plt.ylabel('Cumulative explained variance')
    plt.savefig(os.path.join(outdir, 'cumvariance-PCA-segfree.pdf'))

    # numpber of optical sections
    nstacks = imglist[0].shape[-1]
    eigentiles = segf.pca.components_.reshape((segf.n_components,
                                               *segf.tile_size, nstacks))
    plot_channels([np.max(eigentiles[i], axis=-1) for i in range(segf.n_components)],
              nrow=5, ncol=10, scale_x=2, scale_y=2)
    plt.savefig(os.path.join(outdir, 'eigentiles.pdf'),
                bbox_inches='tight')

    org_prof.index = titles
    sel = VarianceThreshold(threshold=1e-4).fit(org_prof)
    sn.clustermap(org_prof.loc[:,sel.get_support()])
    plt.savefig(os.path.join(outdir, 'fit-heatmap.pdf'))

    javabridge.kill_vm()
