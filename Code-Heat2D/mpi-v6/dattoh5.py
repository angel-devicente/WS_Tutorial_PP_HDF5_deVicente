#!/usr/bin/env python3

# HDF5 uses C storage convention (as Python), while Fortran storage convention
# is different, so we have to transpose the arrays in this script and treat ny
# as the first dimension as opposed to Fortran, where nx is the first dimension.
#
# See https://support.hdfgroup.org/HDF5/doc1.6/UG/12_Dataspaces.html (C versus
# Fortran Dataspaces)

import h5py
import numpy as np

def create_h5(datfile,h5file):
    fh5 = h5py.File(h5file, 'w')
    
    with open(datfile) as f:
        l = f.readline().split()
        nx = int(l[1])
        ny = int(l[2])

        fh5.attrs['nx'] = np.int32([nx])
        fh5.attrs['ny'] = np.int32([ny])
        
        temp = np.zeros([ny,nx])
        cnt = 0
        while l:
            l = f.readline().split()
            if len(l) == 0:
                break
            else:
                for j in range(ny):
                    temp[j,cnt] = float(l[j])
                cnt += 1
    
        fh5['T'] = temp
    
    fh5.close()

    
    
