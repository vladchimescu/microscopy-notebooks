#!/usr//bin/env python3
'''
Script for single-cell morphology feature
extraction from viable leukemia cells
'''
import javabridge
import bioformats as bf
import numpy as np
import pandas as pd
import os
import sys
import h5py
from bioimg import ImgX, load_imgstack, read_bbox

def get_train_instance(path, fname,
                       columns=['ymin', 'xmin', 'ymax', 'xmax'],
                       pad=0):
    imgstack = load_imgstack(fname=os.path.join(path, fname + ".png"),
                             verbose=False)
    img = np.squeeze(imgstack)
    df = pd.read_csv(os.path.join(path, fname + ".csv"))
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
    wellnum = int(sys.argv[3]) - 1
    print("Processing plate: " + str(plate))
    platedir = os.path.join(path, plate)
    imgs = [f.replace('.csv', '') for f in os.listdir(platedir) if '.csv' in f]
    # load plate annotation table
    drug_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                          sep='\t')
    drug_df = drug_df.sort_values(['well', 'Culture']).reset_index(drop=True)
    well = drug_df['well'][wellnum]
    well_imgs = [img for img in imgs if well in img]
    well_imgs.sort()

    # load minimum lysosomal intensity (estmated background)
    thresh_df = pd.read_csv('data/coculture_metafiles/thresholds/'+ plate + '.csv')
    bg_thresh = np.min(thresh_df[['thresh_mono', 'thresh_co']].values)

    imgdata = []
    annot = []

    for im in well_imgs:
        df = pd.read_csv(os.path.join(platedir, im + ".csv"))
        labels_df = df['class'].to_frame()
        labels_df['file'] = im

        img, bbox = get_train_instance(path=platedir,
                                       fname=im)
        # initialize 'ImgX' class
        imgx = ImgX(img=img, bbox=bbox,
                    n_chan=['Lysosomal', 'Calcein', 'Hoechst'])
        imgx.params['texture'] = 'both'
        imgx.compute_props()
        img_df = imgx.get_df().copy()

        # set lysosomal area to zero for those
        # bounding boxes that don't have any pixel intensities
        # above the estimated background
        img_ly = img[:,:,0]
        bg_mask = [not np.any(img_ly[x[2]:x[3], x[0]:x[1]] > bg_thresh) for x in bbox]
        if np.any(bg_mask):
            img_df.loc[bg_mask, 'ch-Lysosomal-area'] = 0
        
        if img_df.shape[0] == labels_df.shape[0]:
            annot.append(labels_df)
            imgdata.append(img_df)

    X_df = pd.concat(imgdata).reset_index(drop=True)
    annot_df = pd.concat(annot).reset_index(drop=True)
    X_out = pd.concat([annot_df, X_df], axis=1)

    outdir = os.path.join('imgdata', plate)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    X_out.to_csv(os.path.join(outdir, well + '.csv'),
                 index=False)
    javabridge.kill_vm()
