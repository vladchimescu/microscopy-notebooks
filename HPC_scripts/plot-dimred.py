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
import seaborn as sn
from sklearn.manifold import TSNE
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def scale_columns(df):
    return (df-df.mean())/df.std()


feat_subset = ['ch-Calcein-area',
               'ch-Calcein-eccentricity',
               'ch-Calcein-mean_intensity',
               'ch-Hoechst-area',
               'ch-Hoechst-eccentricity',
               'ch-Hoechst-mean_intensity',
               'ch-Lysosomal-area',
               'ch-Lysosomal-eccentricity',
               'ch-Lysosomal-mean_intensity']

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

    imgdf = []
    for w in dmso_wells:
        imgdf.append(pd.read_csv(os.path.join(platedir, w+'.csv')))

    imgdf = pd.concat(imgdf).reset_index(drop=True)
    labels = imgdf[['class', 'file']]
    imgdf = imgdf.drop(['class', 'file'], axis=1)
    labels['well'] = labels['file'].replace(regex=r'f[0-9].+', value='')
    labels['class'] = labels['class'].apply(
        lambda x: 'Viable' if x == 2 else 'Apoptotic')
    labels = pd.merge(labels, dmso, on='well')
    imgdf_scaled = scale_columns(imgdf)
    Xfeat = imgdf_scaled.loc[:, feat_subset]
    X_tsne = TSNE(n_components=2, random_state=21,
                  perplexity=50).fit_transform(imgdf_scaled)
    X_df = pd.concat(
        [pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2']), labels], axis=1)
    X_df = pd.concat([X_df, Xfeat], axis=1)

    outdir = 'figures/dimreduction'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fig, ax = plt.subplots(figsize=(9, 8))
    sn.set(font_scale=1.4)
    sn.set_style('white')
    sn.despine()
    sn.scatterplot(x='tsne1', y='tsne2', data=X_df,
                   hue='Culture', style='class',
                   s=40, alpha=0.7)
    plt.legend(loc='lower right', bbox_to_anchor=(1.2, 0.05))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('TSNE 1')
    plt.ylabel('TSNE 2')
    # plt.title('DMSO control wells')
    fig.savefig(os.path.join(outdir, plate + '-DMSO' + '.pdf'),
                bbox_inches='tight')

    img_viab = imgdf.iloc[np.where(labels['class'] == 'Viable')[0], :]
    img_viab = scale_columns(img_viab).reset_index(drop=True)
    labels_viab = (labels[labels['class'] == 'Viable'].
                   reset_index(drop=True))

    X_tsne = TSNE(n_components=2, random_state=21,
                  perplexity=30).fit_transform(img_viab)
    X_viab = pd.concat(
        [pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2']), labels_viab],
        axis=1)
    X_viab = pd.concat([X_viab, img_viab.loc[:, feat_subset]], axis=1)

    fig, ax = plt.subplots(figsize=(7, 7))
    sn.set(font_scale=1.3)
    sn.set_style('white')
    sn.despine()
    sn.scatterplot(x='tsne1', y='tsne2', data=X_viab,
                   hue='Culture',
                   s=60, alpha=0.5)
    plt.legend(loc='lower right', bbox_to_anchor=(1.3, 0.8))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('TSNE 1')
    plt.ylabel('TSNE 2')
    # plt.title('Viable cells')
    fig.savefig(os.path.join(outdir, plate + '-DMSO-viable' + '.pdf'),
                bbox_inches='tight')

    cmap = sn.diverging_palette(240, 15, as_cmap=True)
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(16, 14))
    sn.despine()
    for r in range(3):
        for c in range(3):
            f = feat_subset[3*r+c]
            ax[r, c].set_title(f.replace('ch-', ''))
            maxval = np.max([np.abs(np.percentile(X_viab[f].values, 0.1)),
                             np.abs(np.percentile(X_viab[f].values, 0.9))])
            sc = ax[r, c].scatter(x=X_viab['tsne1'].values,
                                  y=X_viab['tsne2'].values,
                                  c=X_viab[f].values,
                                  cmap=cmap, alpha=0.5,
                                  vmin=-maxval, vmax=maxval)
            cbaxes = inset_axes(ax[r, c], width="3%", height="45%", loc=1)
            plt.colorbar(sc, ax=ax[r, c], cax=cbaxes)
    fig.text(0.5, 0.04, 'TSNE 1', ha='center')
    fig.text(0.04, 0.5, 'TSNE 2', va='center', rotation='vertical')
    fig.subplots_adjust(wspace=0.4)
    fig.savefig(os.path.join(outdir, plate + '-DMSO-viab-features' + '.pdf'),
                bbox_inches='tight')

    drugs = annot_df[np.isin(annot_df.Drug, drug_sel)].reset_index(drop=True)
    drug_wells = drugs['well'].unique()
    imgdf = []
    for w in drug_wells:
        df = pd.read_csv(os.path.join(platedir, w+'.csv'))
        imgdf.append(df[df['class'] == 2])
    imgdf = pd.concat(imgdf).reset_index(drop=True)
    labels = imgdf[['class', 'file']]
    imgdf = imgdf.drop(['class', 'file'], axis=1)
    labels['well'] = labels['file'].replace(regex=r'f[0-9].+', value='')
    labels['class'] = labels['class'].apply(
        lambda x: 'Viable' if x == 2 else 'Apoptotic')
    labels = pd.merge(labels, drugs, on='well')
    X_drug = scale_columns(imgdf)
    Xfeat = X_drug.loc[:, feat_subset]
    X_tsne = TSNE(n_components=2, random_state=21,
                  perplexity=50).fit_transform(X_drug)
    X_df = pd.concat(
        [pd.DataFrame(X_tsne, columns=['tsne1', 'tsne2']), labels], axis=1)
    X_df = pd.concat([X_df, Xfeat], axis=1)

    fig, ax = plt.subplots(figsize=(9, 8))
    sn.set(font_scale=1.4)
    sn.set_style('white')
    sn.despine()
    sn.scatterplot(x='tsne1', y='tsne2', data=X_df,
                   hue='Culture',
                   s=40, alpha=0.8)
    plt.legend(loc='lower right', bbox_to_anchor=(1.2, 0.05))
    plt.xlabel('TSNE 1')
    plt.ylabel('TSNE 2')
    # plt.title('Drug-treated wells')
    fig.savefig(os.path.join(outdir,
                             plate + '-drugs-viab-by-culture' + '.pdf'),
                bbox_inches='tight')

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

    cmap = sn.diverging_palette(240, 15, as_cmap=True)
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(16, 14))
    sn.despine()
    for r in range(3):
        for c in range(3):
            f = feat_subset[3*r+c]
            ax[r, c].set_title(f.replace('ch-', ''))
            maxval = np.max([np.abs(np.percentile(X_df[f].values, 0.1)),
                             np.abs(np.percentile(X_df[f].values, 0.9))])
            sc = ax[r, c].scatter(x=X_df['tsne1'].values,
                                  y=X_df['tsne2'].values,
                                  c=X_df[f].values,
                                  cmap=cmap, alpha=0.5,
                                  vmin=-maxval, vmax=maxval)
            cbaxes = inset_axes(ax[r, c], width="3%", height="45%",
                                loc='upper right')
            plt.colorbar(sc, ax=ax[r, c], cax=cbaxes)

    fig.subplots_adjust(wspace=0.4)
    fig.text(0.5, 0.04, 'TSNE 1', ha='center')
    fig.text(0.04, 0.5, 'TSNE 2', va='center', rotation='vertical')
    fig.savefig(os.path.join(outdir,
                             plate + '-drugs-viab-features' + '.pdf'),
                bbox_inches='tight')
