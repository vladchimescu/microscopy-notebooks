#!/usr/bin/env python3
'''
Script for aggregating single-cell profiles
'''
import numpy as np
import pandas as pd
import os
import sys
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from bioimg.singlecell import preprocess_data
from bioimg.singlecell import select_features, scale_data
from bioimg.singlecell import aggregate_profiles

def load_CLL_cells(platedir, wells, annot, which=[1,2]):
    imgdf = []
    for w in wells:
        if os.path.isfile(os.path.join(platedir, w+'.csv')):
            df = pd.read_csv(os.path.join(platedir, w+'.csv'))
            df['well'] = w
            imgdf.append(df[np.isin(df['class'], which)])
    imgdf = pd.concat(imgdf).reset_index(drop=True)
    labels = imgdf[['class', 'well']]
    imgdf = imgdf.drop(['class', 'well'], axis=1)
    labels['class'] = labels['class'].apply(lambda x: 'Viable' if x == 2 else 'Apoptotic')    
    labels = pd.merge(labels, annot, on='well')
    return imgdf, labels

def normalize_by_control(imgdf, img_annot):
    ctrl_df = imgdf[img_annot['Drug']=='DMSO']
    # center and scale by control wells
    scaler = StandardScaler().fit(ctrl_df)
    imgdf_scaled = scale_data(imgdf, scaler=scaler)
    return imgdf_scaled


if __name__ == '__main__':
    # path to the image data
    path = 'CLLdata/'
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[1]
    print("Processing plate: " + str(plate))

    screenmap = pd.read_csv('data/coculture_metafiles/screen_map.txt',
                            sep='\t')
    # load plate annotation table
    lt = screenmap[screenmap['plate']==plate]['lt'].values[0]
    layout_file = 'CLL_' + lt + '_layout.txt'
    annot_df = pd.read_csv(os.path.join('data/coculture_metafiles', layout_file),
                           sep='\t')
    platedir = os.path.join(path, plate)

    # load all wells in which cells were detected
    all_wells = [f.replace('.csv', '') for f
                 in os.listdir(platedir) if '.csv' in f]    
    all_wells = np.array(all_wells)

    imgdf, img_annot = load_CLL_cells(platedir=platedir, 
                                        wells=all_wells,
                                       annot=annot_df,
                                      which=2)
    sel = VarianceThreshold(threshold=1e-12).fit(imgdf)
    imgdf = preprocess_data(df=imgdf, sel=sel, glog=True)

    mono_df = imgdf[img_annot['Culture']=='Mono-culture']
    mono_annot =  img_annot[img_annot['Culture']=='Mono-culture']
    co_df = imgdf[img_annot['Culture']=='Co-culture']
    co_annot = img_annot[img_annot['Culture']=='Co-culture']

    # normalize mono- and cocultures separately
    mono_scaled = normalize_by_control(imgdf=mono_df, img_annot=mono_annot)
    co_scaled = normalize_by_control(imgdf=co_df, img_annot=co_annot)
    
    imgdf_scaled = pd.concat([mono_scaled, co_scaled]).reset_index(drop=True)
    img_annot = pd.concat([mono_annot, co_annot]).reset_index(drop=True)

    ## aggregated profiles
    img_prof = aggregate_profiles(imgdf_scaled, img_annot)

    outdir = 'data/CLL_profiles/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    img_prof.to_csv(os.path.join(outdir, plate+'.csv'), index=False)
