#!/usr/bin/env python3

import javabridge
import bioformats as bf
import numpy as np
import pandas as pd
import os
import sys
import h5py
from skimage.filters import threshold_otsu
from bioimg import ImgX, load_imgstack, read_bbox

def get_train_instance(path, fname,
                       columns=['ymin', 'xmin', 'ymax', 'xmax'],
                       which=[1,2],
                       pad=0):
    imgstack = load_imgstack(fname=os.path.join(path, fname + ".png"),
                             verbose=False)
    img = np.squeeze(imgstack)
    df = pd.read_csv(os.path.join(path, fname + ".csv"))
    df = df[np.isin(df['class'], which)].reset_index(drop=True)
    rmax, cmax, _ = img.shape
    bbox = read_bbox(df=df, rmax=rmax,
                     cmax=cmax, columns=columns,
                     pad=pad)
    return img, bbox

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)

     # path to the image data
    path = sys.argv[1]
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[2]
    print("Processing plate: " + str(plate))
    platedir = os.path.join(path, plate)
    imgs = [f.replace('.csv', '') for f in os.listdir(platedir) if '.csv' in f]
    # load plate annotation table
    drug_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                          sep='\t')
    drug_df = drug_df.sort_values(['well', 'Culture']).reset_index(drop=True)

    # DMSO wells
    dmso = drug_df[drug_df['Drug']=='DMSO']

    # separately for mono- and coculture
    dmso_mono = dmso[dmso['Culture']=='Mono-culture']['well'].values
    dmso_co = dmso[dmso['Culture']=='Co-culture']['well'].values

    mono_well_imgs = [img for img in imgs for well in dmso_mono if well in img]
    mono_well_imgs.sort()

    co_well_imgs = [img for img in imgs for well in dmso_co if well in img]
    co_well_imgs.sort()

    # find only the lysosomal threshold based on monoculture controls
    ly_thresh = []
    for im in mono_well_imgs:
            img, bbox = get_train_instance(path=platedir,
                                           fname=im)
            img_ly = img[:,:,0]
            ly_thresh += [threshold_otsu(img_ly[x[2]:x[3], x[0]:x[1]]) for x in bbox]
    thresh_mono = np.percentile(ly_thresh, 25)

    # now find the lysosomal threshold based coculture DMSO wells
    # find only the lysosomal threshold
    ly_thresh = []
    for im in co_well_imgs:
            img, bbox = get_train_instance(path=platedir,
                                           fname=im)
            img_ly = img[:,:,0]
            ly_thresh += [threshold_otsu(img_ly[x[2]:x[3], x[0]:x[1]]) for x in bbox]
    thresh_co = np.percentile(ly_thresh, 25)

    thresh_df = pd.DataFrame({'plate': [plate],
              'thresh_mono': [thresh_mono],
              'thresh_co': [thresh_co]})

    outdir = 'data/coculture_metafiles/thresholds'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    thresh_df.to_csv(os.path.join(outdir, plate+'.csv'),
                     index=False)

    javabridge.kill_vm()
