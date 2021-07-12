#!/usr/bin/env python3
'''
Script for aggregating single-cell profiles of BiTE data
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

def load_bite(platedir, wells, annot, thresh=600):
    imgdf = []
    for w in wells:
        if os.path.isfile(os.path.join(platedir, w+'.csv')):
            df = pd.read_csv(os.path.join(platedir, w+'.csv'))
            df['well'] = w
            imgdf.append(df)
    imgdf = pd.concat(imgdf)
    labels = imgdf[['well']]
    imgdf = imgdf.drop(['well'], axis=1) 
    labels = pd.merge(left=labels, right=annot,
                      on='well', how='left')

    # cutoff for viability based on Calcein pixel area
    viab = (imgdf['ch-Calcein-area'] > thresh).values

    imgdf = imgdf[viab].reset_index(drop=True)
    labels = labels[viab].reset_index(drop=True)
    
    return imgdf, labels

# manually selected features
sel = ['ch-Calcein-area', 
       'ch-Calcein-mean_intensity',
       'ch-Calcein-solidity',
       'ch-Calcein-extent',
       'ch-Calcein-filled_area',
       'ch-PE-mean_intensity',
       'ch-PE-area',
       'ch-PE-solidity',
       'ch-PE-extent',
       'ch-PE-filled_area',
       'ch-APC-area',
       'ch-APC-mean_intensity',
       'ch-APC-solidity',
       'ch-APC-extent',
       'ch-APC-filled_area',
       'ch-Hoechst-area',
       'ch-Hoechst-mean_intensity',
       'ch-Hoechst-solidity',
       'ch-Hoechst-extent',
       'ch-Hoechst-filled_area']

if __name__ == '__main__':
    # path to the image data
    path = 'data/clust-BiTE-profiles/'
    # plate identifier (e.g. '180528_Plate3')
    plate = 'BiTE-' + sys.argv[1]
    print("Processing plate: " + str(plate))
    platedir = os.path.join(path, plate)

    # load plate annotation
    annot_df = pd.read_csv('data/BiTE/bite_annot.csv')
    annot_df = annot_df[annot_df['plate']==plate]
    annot_df = annot_df.reset_index(drop=True)

    # load all wells in which cells were detected
    all_wells = [f.replace('.csv', '') for f
                 in os.listdir(platedir) if '.csv' in f]    
    all_wells = np.array(all_wells)

    imgdf, img_annot = load_bite(platedir=platedir, 
                                        wells=all_wells,
                                       annot=annot_df)
    imgdf = preprocess_data(df=imgdf, sel=sel, glog=False)

    ## aggregated profiles
    img_prof = aggregate_profiles(imgdf, img_annot)

    # add cluster count in each well (unnormalized count)
    clust_count = img_annot.groupby('well').size().reset_index()
    clust_count.rename(columns={0 : 'count'}, inplace=True)
    img_prof = pd.merge(img_prof, clust_count, on='well')
    

    outdir = 'data/clustbite_mean_profiles_unnormalized/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    img_prof.to_csv(os.path.join(outdir, plate+'.csv'), index=False)
