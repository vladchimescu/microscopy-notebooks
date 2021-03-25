#!/usr/bin/env python3
import javabridge
import bioformats as bf
import numpy as np
import pandas as pd
import os
import sys
import re
import h5py
from skimage.filters import threshold_otsu
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
    print("Processing plate: " + str(plate))

    # load the screen map
    screenmap = pd.read_csv('data/coculture_metafiles/screen_map.txt',
                            sep='\t')
    # load plate annotation table
    lt = screenmap[screenmap['plate']==plate]['lt'].values[0]
    layout_file = 'CLL_' + lt + '_layout.txt'
    annot_df = pd.read_csv(os.path.join('data/coculture_metafiles', layout_file),
                           sep='\t')

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

    # subset only to viable CLL nuclei
    eval_df = eval_df[eval_df['class']==2].reset_index(drop=True)

    all_wells = eval_df['well'].unique()
    fnames = [f for f in os.listdir(imgdir) if '.tiff' in f]

    # DMSO wells
    dmso = annot_df[annot_df['Drug']=='DMSO']

    dmso_mono = dmso[dmso['Culture']=='Mono-culture']['well'].values
    dmso_co = dmso[dmso['Culture']=='Co-culture']['well'].values

    gamma = 0.3
    # find only the lysosomal threshold based on monoculture controls
    ly_thresh = []
    for well in dmso_mono[:3]:
        for fview in range(1,4):
            wellpos =  well + 'f' + str(fview).zfill(2)
            wfiles = [f for f in fnames if wellpos in f and '(2)' not in f]
            imgstack = load_image_series(path=imgdir,
                                        imgfiles=[w for w in wfiles if 'ch2' in w])
            img_ly = np.max(imgstack, axis=0)
            img_ly = img_ly**gamma
            well_df = eval_df[eval_df['wellpos']==wellpos]
            if well_df.shape[0]:
                rmax, cmax = img_ly.shape
                bbox = read_bbox(df=well_df, rmax=rmax,
                                cmax=cmax,
                                columns=['ymin','xmin','ymax','xmax'],
                                pad=5)
                ly_thresh += [threshold_otsu(img_ly[x[2]:x[3], x[0]:x[1]]) for x in bbox]
    thresh_mono = np.percentile(ly_thresh, 25)

    # find only the lysosomal threshold based on monoculture controls
    ly_thresh = []
    for well in dmso_co[:3]:
        for fview in range(1,4):
            wellpos =  well + 'f' + str(fview).zfill(2)
            wfiles = [f for f in fnames if wellpos in f and '(2)' not in f]
            imgstack = load_image_series(path=imgdir,
                                        imgfiles=[w for w in wfiles if 'ch2' in w])
            img_ly = np.max(imgstack, axis=0)
            img_ly = img_ly**gamma
            well_df = eval_df[eval_df['wellpos']==wellpos]
            if well_df.shape[0]:
                rmax, cmax = img_ly.shape
                bbox = read_bbox(df=well_df, rmax=rmax,
                                cmax=cmax,
                                columns=['ymin','xmin','ymax','xmax'],
                                pad=5)
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
