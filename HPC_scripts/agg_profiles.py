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

def load_viable_cells(platedir, wells, annot):
    imgdf = []
    for w in wells:
        df = pd.read_csv(os.path.join(platedir, w+'.csv'))
        imgdf.append(df[df['class'] == 2])
    imgdf = pd.concat(imgdf).reset_index(drop=True)
    labels = imgdf[['class', 'file']]
    imgdf = imgdf.drop(['class', 'file'], axis=1)
    labels['well'] = labels['file'].replace(regex=r'f[0-9].+', value='')
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
    path = 'imgdata/'
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[1]
    print("Processing plate: " + str(plate))

    patannot = pd.read_csv('data/coculture_metafiles/patannot.txt',
                      sep='\t')
    platedir = os.path.join(path, plate)

    # load plate annotation table
    if patannot[patannot['plate']== plate]['Diagnosis'].values[0] == 'HCL':
        annot_df = pd.read_csv('data/coculture_metafiles/HCL_plate_layout.txt',
                       sep='\t')
    else:
        annot_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                           sep='\t')
    

    # load all wells in which cells were detected
    all_wells = [f.replace('.csv', '') for f
                 in os.listdir(platedir) if '.csv' in f]    
    all_wells = np.array(all_wells)

    imgdf, img_annot = load_viable_cells(platedir=platedir, 
                                        wells=all_wells,
                                       annot=annot_df)
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

    outdir = 'data/coculture_profiles/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    img_prof.to_csv(os.path.join(outdir, plate+'.csv'), index=False)
