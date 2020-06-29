#!/usr//bin/env python3
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


if __name__ == '__main__':
    # path to the image data
    path = 'imgdata/'
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[1]
    print("Processing plate: " + str(plate))
    # load plate annotation table
    annot_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                           sep='\t')
    platedir = os.path.join(path, plate)

    # load all wells
    all_wells = annot_df['well'].values

    imgdf, img_annot = load_viable_cells(platedir=platedir, 
                                        wells=all_wells,
                                       annot=annot_df)
    ctrl_df = imgdf[img_annot['Drug']=='DMSO']
    sel = VarianceThreshold(threshold=1e-12).fit(ctrl_df)
    ctrl_df = select_features(df=ctrl_df, sel=sel)

    imgdf = preprocess_data(df=imgdf, sel=sel, glog=True)
    # center and scale by control wells
    scaler = StandardScaler().fit(ctrl_df)
    imgdf_scaled = scale_data(imgdf, scaler=scaler)

    ## aggregated profiles
    img_prof = aggregate_profiles(imgdf_scaled, img_annot)

    outdir = 'data/coculture_profiles/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    img_prof.to_csv(os.path.join(outdir, plate+'.csv'), index=False)
