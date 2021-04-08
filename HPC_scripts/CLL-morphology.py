#!/usr//bin/env python3
'''
Script for single-cell morphology feature
extraction from segmented CLL nuclei
'''
import javabridge
import bioformats as bf
import numpy as np
import pandas as pd
import os
import sys
import re
import h5py
from ast import literal_eval
from bioimg import ImgX, load_image_series, read_bbox

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)

    # path to CLL data
    # path = '/Volumes/gitlab/microscopy/data/Sophie/'
    path = sys.argv[1]
    # plate identifier
    # plate = '180306_Plate1'
    plate = sys.argv[2]
    wellnum = int(sys.argv[3]) - 1
    print("Processing plate: " + str(plate))

    imgdir = os.path.join(path, 'CLL')
    evdir = os.path.join(path, 'Evaluations')
    
    for f in os.listdir(imgdir):
        if re.search(plate + "_", f):
            screen_id = f

    imgdir = os.path.join(imgdir, screen_id, 'Images')
    evdir = os.path.join(evdir, screen_id)

    evname = [s for s in os.listdir(evdir)
              if re.match('Evaluation', s) is not None]
    evdir = os.path.join(evdir, evname[0])

    # load the bounding box file for the whole plate
    eval_df = pd.read_table(os.path.join(evdir,
                                         'Objects_Population - CLL cell nuclei.txt'), skiprows=9)
    # remove objects with area less than 10 um^2
    eval_df = eval_df[eval_df['CLL cell nuclei - CLL Nuclei Area [µm²]'] > 10].reset_index(drop=True)

    eval_df['class'] = np.where(eval_df['CLL cell nuclei - CLL Nuclei Area [µm²]'] < 23.8, 1, 2)
    eval_df['Row'] = eval_df['Row'].astype(str)
    eval_df['Column'] = eval_df['Column'].astype(str)
    eval_df['Field'] = eval_df['Field'].astype(str)

    eval_df = eval_df.assign(well=lambda x: 'r' + x['Row'].str.zfill(2) + 'c' + x['Column'].str.zfill(2))
    eval_df = eval_df.assign(wellpos=lambda x: x['well'] + 'f' + x['Field'].str.zfill(2))
    eval_df = eval_df[['well','wellpos', 'Bounding Box', 'X', 'Y', 'class']]
    eval_df[['xmin','ymin','xmax','ymax']] = pd.DataFrame(eval_df['Bounding Box'].apply(lambda x: list(literal_eval(x))).tolist())

    all_wells = eval_df['well'].unique()
    well = all_wells[wellnum]
    fnames = [f for f in os.listdir(imgdir) if '.tiff' in f]

    # load minimum lysosomal intensity (estmated background)
    thresh_df = pd.read_csv('data/coculture_metafiles/thresholds/'+ plate + '.csv')
    bg_thresh = np.min(thresh_df[['thresh_mono', 'thresh_co']].values)

    # apply gamma correction
    gamma = 0.3
   
    imgdata = []
    for fview in range(1,4):
       wellpos =  well + 'f' + str(fview).zfill(2)
       wfiles = [f for f in fnames if wellpos in f and '(2)' not in f]
       imgstack = load_image_series(path=imgdir,
                                    imgfiles=[w for w in wfiles if 'ch1' in w])
       hoechst = np.max(imgstack, axis=0)
       hoechst = hoechst**gamma
       
       imgstack = load_image_series(path=imgdir,
                                    imgfiles=[w for w in wfiles if 'ch2' in w])
       ly = np.max(imgstack, axis=0)
       ly = ly**gamma

       well_df = eval_df[eval_df['wellpos']==wellpos]
       if well_df.shape[0]:
           rmax, cmax = hoechst.shape
           bbox = read_bbox(df=well_df, rmax=rmax,
                            cmax=cmax,
                            columns=['ymin','xmin','ymax','xmax'],
                            pad=5)

           imgx = ImgX(img=np.stack([hoechst, ly], axis=-1), 
                bbox=bbox,
                n_chan=['Hoechst', 'Lysosomal'])
           imgx.params['texture'] = 'both'
           imgx.compute_props()
           img_df = imgx.get_df().copy()
           img_df['class'] = well_df['class'].values
           bg_mask = [not np.any(ly[x[2]:x[3], x[0]:x[1]] > bg_thresh) for x in bbox]
           if np.any(bg_mask):
               img_df.loc[bg_mask, 'ch-Lysosomal-area'] = 0
           imgdata.append(img_df)


    X_df = pd.concat(imgdata).reset_index(drop=True)
    
    outdir = os.path.join('CLLdata', plate)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    X_df.to_csv(os.path.join(outdir, well + '.csv'),
                 index=False)
    javabridge.kill_vm()
