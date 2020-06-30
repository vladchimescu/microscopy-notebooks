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
from bioimg.singlecell import plot_dimred, facet_dimred, facet_density, facet_boxplot

feat_subset = ['ch-Hoechst-mean_intensity',
               'ch-Lysosomal-mean_intensity',
               'ch-Calcein-convex_area',
               'ch-Hoechst-SumAverage-d3-0',
               'ch-Calcein-eccentricity',
               'ch-Lysosomal-InfoMeas1-d7-0',
               'ch-Hoechst-weighted_moments-0-1',
               'ch-Hoechst-InfoMeas1-d7-3',
               'ch-Lysosomal-Contrast-d7-3']

drug_sel = ['Tofacitinib', 'Midostaurin',
            'Ganetespib', 'Lenalidomide',
            'Pyridone 6', 'UMI-77',
            'Bafilomycin A1',
            'Quizartinib', 'Hydroxychloroquine',
            'Fludarabine', 'Vorinostat',
            'Thioguanine', 'Nutlin 3a',
            'Palbociclib', 'Carfilzomib',
            'JQ1', 'Cytarabine',
            'BAY61-3606', 'Venetoclax',
            'Ixazomib']

drugs_boxplot = ['Tofacitinib', 'Midostaurin',
                   'Ganetespib', 'Lenalidomide',
                   'Pyridone 6', 'UMI-77',
                   'Bafilomycin A1', 
                   'Fludarabine', 'Vorinostat',
                   'Palbociclib', 'Carfilzomib',
                   'JQ1', 'BAY61-3606', 'Venetoclax',
                   'Ixazomib']

feat_boxplot = ['ch-Hoechst-mean_intensity', 'ch-Lysosomal-mean_intensity',
            'ch-Calcein-convex_area', 'ch-Calcein-eccentricity',
            'ch-Lysosomal-Contrast-d7-3', 'ch-Hoechst-InfoMeas1-d7-3']

def load_cells(platedir, wells, annot, which=[1,2]):
    imgdf = []
    for w in wells:
        df = pd.read_csv(os.path.join(platedir, w+'.csv'))
        imgdf.append(df[np.isin(df['class'], which)])
    imgdf = pd.concat(imgdf).reset_index(drop=True)
    labels = imgdf[['class', 'file']]
    imgdf = imgdf.drop(['class', 'file'], axis=1)
    labels['well'] = labels['file'].replace(regex=r'f[0-9].+', value='') 
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
    dmso = annot_df[annot_df.Drug == 'DMSO'].reset_index(drop=True)
    dmso_wells = dmso['well'].unique()

    with open('data/featselect.json') as file:
        featdict = json.load(file)

    '''T-SNE of control wells for both viable and apoptotic cells
       ----------------------------------------------------------
    '''
    imgdf, labels = load_cells(platedir=platedir, wells=dmso_wells,
                               annot=annot_df)
    labels['class'] = labels['class'].apply(
        lambda x: 'Viable' if x == 2 else 'Apoptotic')
    #labels = pd.merge(labels, dmso, on='well')
    sel = featdict['repcor']
    imgdf = preprocess_data(df=imgdf, sel=sel, glog=True)
    scaler = StandardScaler().fit(X=imgdf)
    imgdf_scaled = scale_data(imgdf, scaler=scaler)
    pcs = PCA(n_components=20).fit_transform(imgdf_scaled)
    X_tsne = TSNE(n_components=2, random_state=21, perplexity=50).fit_transform(pcs)
    X_df = pd.concat(
        [pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2']), labels], axis=1)
    #X_df = pd.concat([X_df, Xfeat], axis=1)

    outdir = 'figures/dimreduction'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plot_dimred(X_df, 
            hue='Culture',
            style='class',
            title='DMSO control wells',
            style_order=['Viable', 'Apoptotic'])
    # plt.title('DMSO control wells')
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

    # subset to features selected by residual correlation method
    X_subset = ctrl_scaled[featdict['residcor']]
    pcs = PCA(n_components=15).fit_transform(X_subset)
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

    drugs = annot_df[np.isin(annot_df.Drug, drug_sel)].reset_index(drop=True)
    drug_wells = drugs['well'].unique()

    # which=2 selects only 'Viable' cancer cells
    drug_df, drug_annot = load_cells(platedir=platedir,
                                     wells=drug_wells,
                                     annot=annot_df, which=2)

    drug_df = preprocess_data(df=drug_df, sel=sel, glog=True)
    # scale by control wells
    drug_scaled = scale_data(drug_df, scaler=scaler)

    '''T-SNE of viable cells in drug wells
       -----------------------------------
    '''
    pcs = PCA(n_components=15).fit_transform(drug_scaled[featdict['residcor']])
    X_tsne = TSNE(n_components=2, random_state=21, perplexity=50).fit_transform(pcs)

    X_df = pd.concat([pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2']),
                      drug_annot], axis=1)
    X_df = pd.concat([X_df, drug_scaled.loc[:,feat_subset]], axis=1)

    plot_dimred(X_df, hue='Culture',
           title='Viable cells in drug-treated wells')
    plt.savefig(os.path.join(outdir,
                             plate + '-drugs-viab-by-culture' + '.pdf'),
                bbox_inches='tight')

    '''T-SNE of viable cells (drugs) colored by features
       -------------------------------------------------
    '''
    facet_dimred(X_df, feat_subset=feat_subset,
            nrows=3, ncols=3)
    plt.savefig(os.path.join(outdir,
                             plate + '-drugs-viab-features' + '.pdf'),
                bbox_inches='tight')

    '''Density plots of features in drug wells
       ---------------------------------------
    '''
    colsub = ['Culture'] + feat_subset
    X_long = pd.melt(X_df[colsub], 
                     id_vars=['Culture'], 
                     value_vars=feat_subset,
                     var_name='feature', 
                     value_name='val')
    facet_density(X_long, feat_column='feature',
             ncols=3, hue='Culture', sharey=False,
             aspect=1.4)
    plt.savefig(os.path.join(outdir, plate + '-drugs-featdist' + '.pdf'),
                bbox_inches='tight')

    '''Boxplots stratified by drug and culture
       ---------------------------------------
    '''
    colsub = ['Drug', 'Culture'] + feat_subset
    X_long = pd.melt(X_df[colsub], 
                     id_vars=['Drug', 'Culture'], 
                     value_vars=feat_subset,
                     var_name='feature', 
                     value_name='val')
    
    X_long = X_long.loc[np.isin(X_long['Drug'], drugs_boxplot),:]
    X_long = X_long.loc[np.isin(X_long['feature'], feat_boxplot),:]
    facet_boxplot(X_long, x='Drug',
              y='val', feat_column='feature',
              ncols=3,
              nrows=2, hue='Culture')
    plt.savefig(os.path.join(outdir, plate + '-drugs-boxplot' + '.pdf'),
                bbox_inches='tight')


    '''T-SNE colored by drug
       ---------------------
    '''
    drug_chunks = [drug_sel[i:i + 5] for i in range(0, len(drug_sel), 5)]
    fig, ax = plt.subplots(ncols=2, nrows=2,
                           figsize=(14, 12))
    sn.set(font_scale=1.2)
    sn.set_style('white')
    sn.despine()
    for r in range(2):
        for c in range(2):
            sn.scatterplot(x='tsne1', y='tsne2',
                           data=X_df[np.isin(
                               X_df['Drug'], drug_chunks[r*2+c])],
                           hue='Drug',
                           s=40, alpha=0.8, ax=ax[r, c])
            ax[r, c].legend(loc='lower right', bbox_to_anchor=(1.4, 0.7))
            ax[r, c].set_xlabel('TSNE 1')
            ax[r, c].set_ylabel('TSNE 2')
    fig.subplots_adjust(wspace=0.5)
    fig.savefig(os.path.join(outdir,
                             plate + '-drugs-viab-by-drug' + '.pdf'),
                bbox_inches='tight')
