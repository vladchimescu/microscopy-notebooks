#!/usr//bin/env python3
'''
Script for generating MIP images for non-CLL diseases
'''
# %%
import javabridge
import bioformats as bf
import numpy as np
import pandas as pd
import os
import sys
import re
import h5py
from ast import literal_eval
from bioimg import ImgX, load_image_series, write_image

# %%

if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)

    # path to CLL data
    # path = '/Volumes/gitlab/microscopy/data/Sophie/'
    path = sys.argv[1]
    # plate identifier
    # plate = '180306_Plate1'
    plate = sys.argv[2]
    wellnum = int(sys.argv[3]) - 1
    print("Processing plate: " + str(plate))

    #imgdir = os.path.join(path, 'CLL')
    imgdir = path
    
    for f in os.listdir(imgdir):
        if re.search(plate + "_", f):
            screen_id = f

    imgdir = os.path.join(imgdir, screen_id, 'Images')

    # load plate annotation table
    annot_df = pd.read_csv('data/AML_trainset/drugannot.txt',
                           sep='\t')
    annot_df = annot_df.sort_values('well').reset_index(drop=True)

    all_wells = annot_df['well'].unique()
    well = all_wells[wellnum]
    fnames = [f for f in os.listdir(imgdir) if '.tiff' in f]

    outdir = os.path.join('data/MIP_images', plate)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for fview in range(1,4):
       wellpos =  well + 'f' + str(fview).zfill(2)
       wfiles = [f for f in fnames if wellpos in f and '(2)' not in f]
       imgstack = load_image_series(path=imgdir,
                                    imgfiles=[w for w in wfiles if 'ch1' in w])
       hoechst_mip = np.max(imgstack, axis=0)
       write_image(img=hoechst_mip, path=os.path.join(outdir, wellpos + '-ch1.tiff'))

       imgstack = load_image_series(path=imgdir,
                                    imgfiles=[w for w in wfiles if 'ch2' in w])
       ly_mip = np.max(imgstack, axis=0)
       write_image(img=ly_mip, path=os.path.join(outdir, wellpos + '-ch2.tiff'))

       imgstack = load_image_series(path=imgdir,
                                    imgfiles=[w for w in wfiles if 'ch3' in w])
       ca_mip = np.max(imgstack, axis=0)
       write_image(img=ca_mip, path=os.path.join(outdir, wellpos + '-ch3.tiff'))

    javabridge.kill_vm()
