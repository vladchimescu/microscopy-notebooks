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


if __name__ == '__main__':
    # path to the image data
    path = 'CLLdata/'
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[1]
    print("Processing plate: " + str(plate))
    # load plate annotation table
    annot_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                           sep='\t')
    platedir = os.path.join(path, plate)

    # load all wells in which cells were detected
    all_wells = [f.replace('.csv', '') for f
                 in os.listdir(platedir) if '.csv' in f]    
    all_wells = np.array(all_wells)

    imgdf, img_annot = load_CLL_cells(platedir=platedir, 
                                        wells=all_wells,
                                       annot=annot_df)

    viab_df = img_annot.groupby(['class', 'well', 'Drug', 'conc', 'Culture'],as_index=False).agg('size').reset_index(name='n')
    # compute raw viability

    outdir = 'data/CLL_viab/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    viab_df.to_csv(os.path.join(outdir, plate+'.csv'), index=False)
