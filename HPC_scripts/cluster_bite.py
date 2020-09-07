#!/usr/bin/env python3
'''
Morphological profiling of BiTE coculture
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
from skimage.morphology import binary_dilation, disk
from skimage.feature import shape_index, blob_log
from skimage.measure import label, regionprops_table
from bioimg.segment.cv_methods import filter_segm
from bioimg import read_bbox
from bioimg import load_image_series, threshold_img, ImgX

if __name__ == "__main__":
    javabridge.start_vm(class_path=bf.JARS)

    #path = "/g/huber/users/vkim/gitlab/microscopy/data/Tobias/newscreen/"
    path = sys.argv[1]
    plate = sys.argv[2]
    print("Processing plate: " + str(plate))
    
    for f in os.listdir(path):
        if re.search(plate + "_", f):
            screen_id = f

    platedir = os.path.join(path, screen_id, 'Images')
    print("Image path: " + str(platedir))

    well_id = int(sys.argv[3])-1

    fnames = [f for f in os.listdir(platedir) if '.tiff' in f]
    all_wells = list(set([re.search('r[0-9]+c[0-9]+',f).group(0) for f in fnames]))
    all_wells.sort()
    
    well = all_wells[well_id]
    well_imgs = [f for f in fnames if well in f and 'ch2' not in f]

    
    disk_size = 5
    bounds = {'area': (1400, np.inf)}
    
    imgdata = []
    for i in range(1,5):
        wellpos = [f for f in well_imgs if 'f0' + str(i)+ 'p' in f]
        imgseries = load_image_series(path=platedir, imgfiles=wellpos)
        imgseries = imgseries.reshape((8, 4, 2160,2160))
        mipseries = np.amax(imgseries, axis=0)
        # get calcein channel
        ca = mipseries[3]
        # threshold the Calcein channel
        ca_th = threshold_img(ca, method='otsu', binary=True)
        ca_dil = binary_dilation(ca_th, disk(disk_size))

        segm = filter_segm(img=ca,
                    labels=label(ca_dil, connectivity=1),
                    bounds=bounds)
        if not np.all(segm==0):
            feats = regionprops_table(segm,
                           intensity_image=ca,
                          properties=['bbox'])
            df = pd.DataFrame(feats)
            df = df.rename(columns={'bbox-0': 'ymin',
                               'bbox-1': 'xmin',
                               'bbox-2': 'ymax',
                               'bbox-3': 'xmax'})
            imax = ca.shape[0] - 1
            bbox = read_bbox(df=df, rmax=imax, cmax=imax, pad=0)

            imgx = ImgX(img=mipseries.swapaxes(0,-1),
                        bbox=bbox,
                        n_chan=['Hoechst', 'PE',
                        'APC', 'Calcein'])
            imgx.compute_props()
            img_df = imgx.get_df().copy()
            imgdata.append(img_df)

    X_df = pd.concat(imgdata).reset_index(drop=True)

    outdir = "data/clust-BiTE-profiles/BiTE-" +  plate
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    X_df.to_csv(os.path.join(outdir, well + '.csv'),
                 index=False)     
    javabridge.kill_vm()
