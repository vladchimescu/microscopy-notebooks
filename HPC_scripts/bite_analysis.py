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
from skimage.feature import shape_index, blob_log
from bioimg import load_image_series, threshold_img, ImgX

def nantonum(img, pad=-1):
    img_r = np.copy(img)
    img_r[np.isnan(img_r)] = pad
    return img_r

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

    gamma = 0.5
    pad = 5
    
    imgdata = []
    for i in range(1,5):
        wellpos = [f for f in well_imgs if 'f0' + str(i)+ 'p' in f]
        imgseries = load_image_series(path=platedir, imgfiles=wellpos)
        imgseries = imgseries.reshape((8, 4, 2160,2160))
        mipseries = np.amax(imgseries, axis=0)
        hoechst = mipseries[0]
        
        img_th = threshold_img(hoechst**gamma, method='otsu', binary=False)
        img_s = shape_index(img_th)
        img_enh = nantonum(img_s, pad=-1)
        # run blob detection on the shape-index enhanced image
        blobs = blob_log(img_enh,
                         min_sigma=10,
                         max_sigma=14,
                         threshold=0.05)
        if len(blobs):
            bbox = np.stack([np.array([bl[1] - bl[2] - pad,
                               bl[1] + bl[2] + pad,
                               bl[0] - bl[2] - pad,
                               bl[0] + bl[2] + pad]) for bl in blobs])
            imax = hoechst.shape[0] - 1
            bbox[bbox < 0] = 0
            bbox[bbox > imax ] = imax
            bbox = bbox.astype(int)
            imgx = ImgX(img=mipseries.swapaxes(0,-1),
                        bbox=bbox,
                        n_chan=['Hoechst', 'PE',
                        'APC', 'Calcein'])
            imgx.compute_props()
            img_df = imgx.get_df().copy()
            imgdata.append(img_df)

    X_df = pd.concat(imgdata).reset_index(drop=True)

    outdir = "data/BiTE-profiles/BiTE-" +  plate
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    X_df.to_csv(os.path.join(outdir, well + '.csv'),
                 index=False)     
    javabridge.kill_vm()
