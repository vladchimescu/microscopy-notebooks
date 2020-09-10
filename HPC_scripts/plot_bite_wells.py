#!/usr/bin/env/python3
'''Script for outputting BiTE well images
'''

import matplotlib.pyplot as plt
import numpy as np
import javabridge
import bioformats as bf
import pandas as pd
import seaborn as sn
import random
import sys
import os
import re
import csv
from bioimg import load_image_series
from bioimg import combine_channels, plot_channels

def load_bite_image(well, fview=1):
    well_imgs = [f for f in fnames if well in f and 'ch2' not in f]
    wellpos = [f for f in well_imgs if 'f0'+ str(fview)+ 'p' in f]
    imgseries = load_image_series(path=platedir, imgfiles=wellpos)
    imgseries = imgseries.reshape((8, 4, 2160,2160))
    mipseries = np.amax(imgseries, axis=0)
    return mipseries

col_params = dict(colors=['red', 
                          'yellow',
                          'green'],
                             blend = [ 1, 1.2, 1.5],
                             gamma = [0.4, 0.3, 0.6])


if __name__ == '__main__':
    javabridge.start_vm(class_path=bf.JARS)
    
    path = sys.argv[1]
    well_num = int(sys.argv[2]) - 1

    annot_df = pd.read_csv('data/BiTE/bite_annot.csv')
    samples = (annot_df[['plate', 'PatientID']].
               drop_duplicates().
               reset_index(drop=True))
    topwells = pd.read_csv('data/bite-topwells.csv')
    s = topwells['PatientID'].values[well_num]
    drug = topwells['comb'].values[well_num]

    plate = samples[samples['PatientID']==s]['plate'].values[0]
    plate = plate.replace('BiTE-', '')
    print("Processing plate: " + str(plate))
    
    for f in os.listdir(path):
        if re.search(plate + "_", f):
            screen_id = f

    platedir = os.path.join(path, screen_id, 'Images')
    print("Image path: " + str(platedir))

    fnames = [f for f in os.listdir(platedir) if '.tiff' in f]
    annot_samp = annot_df[annot_df['PatientID']==s]

    ctrl_well = annot_samp[annot_samp['comb']=='DMSO']['well'].values[1]
    ctrl_imglist = []
    for i in range(1,5):
        mipseries = load_bite_image(well=ctrl_well, fview=i)
        ctrl_image = combine_channels([mipseries[i] for i in range(1,4)],
                             **col_params)
        ctrl_imglist.append(ctrl_image)

    comb_well = annot_samp[annot_samp['comb']==drug]['well'].values[0]
    bite_imglist = []
    for i in range(1,5):
        mipseries = load_bite_image(well=comb_well, fview=i)
        bite_image = combine_channels([mipseries[i] for i in range(1,4)],
                                 **col_params)
        bite_imglist.append(bite_image)
        
    outdir = "figures/BiTE-topwells"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    imglist = ctrl_imglist + bite_imglist
    titles = ['Control']*4 + [drug]*4
    plot_channels(imglist,
                  nrow=2, ncol=4,
                  titles=titles)
    plt.savefig(os.path.join(outdir, s + '-' + drug + '.pdf'),
                bbox_inches='tight')

    javabridge.kill_vm()
