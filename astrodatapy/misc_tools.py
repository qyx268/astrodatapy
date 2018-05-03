#!/usr/bin/env python               
# -*- coding: utf-8 -*-             
import numpy as np

def filename_to_redshift(fname):
    interger = fname[fname.rfind('z')+1:fname.rfind('pt')]    
    decimal  = fname[fname.rfind('pt')+2:fname.rfind('.dat')] 
    return float(interger) + float(decimal) / 10**len(decimal)

def available_redshifts(feature):
    import os
    import fnmatch
    fileDir = os.path.dirname(os.path.abspath(__file__))
    all_fnames = []
    for root, dirnames, filenames in os.walk(fileDir+'/data/%s'%feature):
        for filename in fnmatch.filter(filenames, '*'):
            if 'dat' in filename:
                all_fnames.append(filename)

    #all_fnames = glob.glob('data/%s/**/*'%feature)
    if len(all_fnames)>0:
        #all_fnames = [f[f.rfind('z'):] for f in all_fnames]
        return np.array([filename_to_redshift(fname) for fname in set(all_fnames)])
    else:
        return np.array([])
