#!/usr/bin/env python3
"""
Script for generating initial segmentation maps
in leukemia entities in HPC environment
"""
import javabridge
import bioformats as bf
import skimage
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import pandas as pd
import os
import sys
from skimage.morphology import disk, binary_erosion
from skimage.morphology import white_tophat
from skimage.morphology import remove_small_objects
from skimage.measure import label
from skimage.measure import regionprops
sys.path.append('..')
from base.utils import load_imgstack
from transform.process import threshold_img
from segment.cv_methods import filter_segm
from base.future_versions import regionprops_table


def segment_connected_components(img):
    '''Initial segmentation for leukemia coculture
       -------------------------------------------
       Use simple connected component labelling to
       segment the nuclei of leukemia and stroma cells.
       The initial segmentation map preserves only regions
       within the specified area and size range

       Parameters
       ----------
       img : array
           Intensity image

       Returns
       -------
       segm_out : array
           Image array with labelled regions
    '''
    img_th = binary_erosion(threshold_img(img, method='otsu',
                                          binary=True), disk(5))
    segm = label(img_th, connectivity=1)
    # change the bound values if necessary
    bounds = {'area': (500, 6000), 'perimeter': (100, 1000)}
    segm1 = filter_segm(img=img, labels=segm, bounds=bounds)
    # change 'bounds' values if necessary
    big = filter_segm(img=img, labels=segm, bounds={'area': (6000, np.inf)}) +\
        filter_segm(img=img, labels=segm, bounds={
                    'perimeter': (1000, np.inf)})
    if len(big[big != 0]) == 0:
        return segm1
    big_obj = img*np.isin(segm, np.unique(big[big != 0]))
    img_tophat = white_tophat(big_obj, disk(25))
    # bright spots from the large regions that were filtered out
    # in the previous step
    segm2 = remove_small_objects(threshold_img(img_tophat,
                                               method='yen', binary=True),
                                 min_size=500)
    # only non-background pixels
    segm1 = (segm1 != 0)
    segm_out = label(np.logical_or(segm1, segm2))
    return segm_out


def get_feature_table(segm, img):
    '''Return feature table for labelled regions
       -----------------------------------------
       Computes region properties for the labelled regions
       and returns all features except 'image', 'filled_image',
       'convex_image', 'extent' and 'coords'.

       Parameters
       ----------
       segm : array
           Labelled image array
       img : array
           Intensity image

       Returns
       -------
       feat_df : DataFrame
           Table with region properties with 'image' and 'filled_image'
           excluded.
    '''
    feats_out = regionprops(label_image=segm, intensity_image=img)
    keys = [k for k in feats_out[0]]
    exclude = ['convex_image', 'coords', 'extent',
               'filled_image', 'image']
    selected_keys = list(set(keys) - set(exclude))
    # sort by key lexicographically
    selected_keys.sort()

    feat_dict = regionprops_table(segm,
                                  intensity_image=img,
                                  properties=selected_keys)
    feat_df = pd.DataFrame(feat_dict)
    return feat_df


if __name__ == "__main__":
    javabridge.start_vm(class_path=bf.JARS)

    # path to the image data
    path = sys.argv[1]
    # plate identifier (e.g. '180528_Plate3')
    plate = sys.argv[2]
    print("Processing plate: " + str(plate))

    # image name
    fname = sys.argv[3] + '.tiff'

    imgname = os.path.join(path, plate, fname)
    imgstack = load_imgstack(fname=imgname)

    # remove a 'dummy' z-axis
    img = np.squeeze(imgstack)
    # for gamma correction
    gamma = 0.3
    # nuclei
    hoechst = img[:, :, 0]**gamma

    # obtain initial segmentation map for
    # leukemia and stroma nuclei
    segm_out = segment_connected_components(img=hoechst)

    # get the region peroperties
    feat_df = get_feature_table(segm=segm_out, img=hoechst)
    fout = os.path.join(path, plate, fname.replace('.tiff', '.csv'))
    feat_df.to_csv(fout, index=False)

    javabridge.kill_vm()
