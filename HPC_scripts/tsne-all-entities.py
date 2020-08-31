#!/usr//bin/env python3
'''
Script for plotting T-SNE and UMAP based
on image features
'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
import json
import seaborn as sn
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from bioimg.singlecell import preprocess_data
from sklearn.preprocessing import StandardScaler
from bioimg.singlecell import scale_data, check_data
from bioimg.singlecell import plot_dimred, facet_dimred, facet_density

feat_subset = ['ch-Hoechst-mean_intensity',
               'ch-Lysosomal-mean_intensity',
               'ch-Hoechst-eccentricity',
               'ch-Lysosomal-InfoMeas1-d7-0',
               'ch-Hoechst-solidity',
               'ch-Hoechst-InfoMeas1-d7-3',
               'ch-Lysosomal-Contrast-d7-3',
               'ch-Lysosomal-area',
               'ch-Lysosomal-extent']

def load_nonCLL_cells(platedir, wells, annot, which=[1,2]):
    imgdf = []
    for w in wells:
        if os.path.exists(os.path.join(platedir, w+'.csv')):
            df = pd.read_csv(os.path.join(platedir, w+'.csv'))
            imgdf.append(df[np.isin(df['class'], which)])
    imgdf = pd.concat(imgdf).reset_index(drop=True)
    labels = imgdf[['class', 'file']]
    imgdf = imgdf.drop(['class', 'file'], axis=1)
    labels['well'] = labels['file'].replace(regex=r'f[0-9].+', value='') 
    labels = pd.merge(labels, annot, on='well')
    labels['class'] = labels['class'].apply(
        lambda x: 'Viable' if x == 2 else 'Apoptotic')
    return imgdf, labels

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
    # path to the image data (different for CLL and non-CLL entities)
    path = sys.argv[1]
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[2]
    print("Processing plate: " + str(plate))

    if 'CLL' in path:
        load_cells = load_CLL_cells
    else:
        load_cells = load_nonCLL_cells
    
    # load plate annotation table
    annot_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                           sep='\t')
    platedir = os.path.join(path, plate)
    dmso = annot_df[annot_df.Drug == 'DMSO'].reset_index(drop=True)
    dmso_wells = dmso['well'].unique()

    with open('data/featselect.json') as file:
        featdict_noncll = json.load(file)

    with open('data/CLL-featselect.json') as file:
        featdict_cll = json.load(file)

    # union of selected features in CLL and non-CLL entities
    sel_feat = list(set(featdict_noncll['residcor'] + featdict_cll['residcor']))
    # remove Calcein features
    sel = [f for f in sel_feat if 'Calcein' not in f]
    sel.sort()

    '''T-SNE of control wells for both viable and apoptotic cells
       ----------------------------------------------------------
    '''
    imgdf, labels = load_cells(platedir=platedir, wells=dmso_wells,
                               annot=annot_df)
    imgdf = preprocess_data(df=imgdf, sel=sel, glog=True)
    scaler = StandardScaler().fit(X=imgdf)
    imgdf_scaled = scale_data(imgdf, scaler=scaler)
    pcs = PCA(n_components=20).fit_transform(imgdf_scaled)
    X_tsne = TSNE(n_components=2, random_state=21, perplexity=50).fit_transform(pcs)
    X_df = pd.concat(
        [pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2']), labels], axis=1)

    outdir = 'figures/coculture/tsnes/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plot_dimred(X_df, 
            hue='Culture',
            style='class',
            title='DMSO control wells',
            style_order=['Viable', 'Apoptotic'])
    plt.savefig(os.path.join(outdir, plate + '-DMSO' + '.pdf'),
                bbox_inches='tight')

    '''T-SNE of viable cells in control wells
       --------------------------------------
    '''
    # which=2 selects only 'Viable' cancer cells
    ctrl_df, ctrl_annot = load_cells(platedir=platedir,
                                     wells=dmso_wells,
                                     annot=annot_df, which=2)
    ctrl_df = preprocess_data(df=ctrl_df, sel=sel, glog=True)

    scaler = StandardScaler().fit(X=ctrl_df)
    ctrl_scaled = scale_data(ctrl_df, scaler=scaler)

    pcs = PCA(n_components=15).fit_transform(ctrl_scaled)
    X_tsne = TSNE(n_components=2, random_state=34, perplexity=30).fit_transform(pcs)
    X_ctrl = pd.concat([pd.DataFrame(X_tsne,
                                     columns=['tsne1', 'tsne2']),
                    ctrl_annot], axis=1)

    plot_dimred(X_ctrl, hue='Culture',
           title='Viable cells in control wells')
    plt.savefig(os.path.join(outdir, plate + '-DMSO-viable' + '.pdf'),
                bbox_inches='tight')

    '''T-SNE of viable cells (DMSO) colored by morphological features
       --------------------------------------------------------------
    '''
    X_ctrl = pd.concat([X_ctrl, ctrl_scaled.loc[:,feat_subset]], axis=1)
    facet_dimred(X_ctrl, feat_subset=feat_subset,
            nrows=3, ncols=3)
    plt.savefig(os.path.join(outdir, plate + '-DMSO-viab-features' + '.pdf'),
                bbox_inches='tight')

    '''Density plots of features in DMSO wells
       ---------------------------------------
    '''
    colsub = ['Culture'] + feat_subset
    X_long = pd.melt(X_ctrl[colsub], 
                     id_vars=['Culture'], 
                     value_vars=feat_subset,
                     var_name='feature', 
                     value_name='val')
    facet_density(X_long, feat_column='feature',
             ncols=3, hue='Culture',
             sharey=False, size=(10,8),
             aspect=1.4)
    plt.savefig(os.path.join(outdir, plate + '-DMSO-featdist' + '.pdf'),
                bbox_inches='tight')
