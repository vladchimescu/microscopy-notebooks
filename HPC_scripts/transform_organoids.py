#!/usr/bin/env python3
'''
Script for fitting unsupervised segmentation-free model
on 3D organoid compound screen images
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
    plate = sys.argv[1]
    well_id = int(sys.argv[2]) - 1

    segf = pickle.load(open("segf_organoids.pkl", "rb"))

    if plate == 'A':
        paths = [os.path.join(datadir, 'CL34A_4x_2p1mm_100umsteps_day0_001/data'),
                 os.path.join(datadir, 'CL34A_4x_2p1mm_100umsteps_day3_001/data')]
    if plate == 'B':
        paths = [os.path.join(datadir, 'CL34B_4x_2p1mm_100umsteps_day0_001/data'),
                 os.path.join(datadir, 'CL34B_4x_2p1mm_100umsteps_day3_001/data')]
   
    outdir = 'data/segfree/organoids'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    imglist = []
    for path in paths:
        fnames, all_wells = get_all_wells(path=path)
        # randomly sample 10 wells
        well = all_wells[well_id]
        well_files = [f for f in fnames if well in f]
        for i in range(4):
            imgstack = load_image_series(path=path, imgfiles=[w for w in well_files if 'P0000' + str(i+1) in w])
            imgstack = imgstack.swapaxes(0,-1)
            imglist.append(imgstack)

    org_prof = segf.transform(imglist)
    cond = re.search('(.+)(--W[0-9]+)', well_files[0]).group(1)
    org_prof.index = [cond + '-day0']*4 + [cond + '-day3']*4

    df_med = org_prof.groupby(org_prof.index).agg('median')
    df_med.columns = ['-'.join(['median', col]) for col in df_med.columns.values]

    df_mean = org_prof.groupby(org_prof.index).agg('mean')
    df_mean.columns = ['-'.join(['mean', col]) for col in df_mean.columns.values]

    df_min = org_prof.groupby(org_prof.index).agg('min')
    df_min.columns = ['-'.join(['min', col]) for col in df_min.columns.values]

    df_max = org_prof.groupby(org_prof.index).agg('max')
    df_max.columns = ['-'.join(['max', col]) for col in df_max.columns.values]

    df_out = pd.concat([df_med, df_mean, df_max, df_min], axis=1)
    df_out.to_csv(os.path.join(outdir, 'plate' + plate + '-' + well+'.csv' ))
    
    javabridge.kill_vm()
