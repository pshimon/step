#! python
########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
import numpy as np
ARRAY_SHAPE_LENGTH=8
Chr_LBL	=0
Flt_LBL	=16
Dbl_LBL=32
Cflt_LBL=48
Cdbl_LBL=64
Int_LBL	=80

def read_header(filename):
    fid = open(filename, 'rb')
    h=np.fromfile(fid,np.int32,ARRAY_SHAPE_LENGTH)
    r=h[0]&15
    print('rank=',r)
    t=h[0]-r
    shp=h[r:0:-1]
    print('shape',shp)
    if t==Chr_LBL:
        print('char')
    elif t==Flt_LBL:
        print('float')
    elif t==Dbl_LBL:
        print('double')
    elif t==Cflt_LBL:
        print('complex float')
    elif t==Cdbl_LBL:
        print('complex double')
    elif t==Int_LBL:
        print('integer')
    else :
        print('unknown type')
    return h

def read_databuf(filename):
    fid = open(filename, 'rb')
    h=np.fromfile(fid,np.int32,ARRAY_SHAPE_LENGTH)
    r=h[0]&15
    t=h[0]-r
    shp=h[r:0:-1]
    n=1
    for nn in shp:
        n=n*nn

    if t==Chr_LBL:
        a = np.fromfile(fid,np.int8,n)
    elif t==Flt_LBL:
        a = np.fromfile(fid,np.float32,n)
    elif t==Dbl_LBL:
        a = np.fromfile(fid,np.float64,n)
    elif t==Cflt_LBL:
        a = np.fromfile(fid,np.complex64,n)
    elif t==Cdbl_LBL:
        a = np.fromfile(fid,np.complex128,n)
    elif t==Int_LBL:
        a=np.fromfile(fid,np.int32,n)
    else :
        print('unknown type')
        return h

    return np.reshape(a,shp)


