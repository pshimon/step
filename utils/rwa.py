#! python
########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
import numpy as np
def write_arr2(a,filename):
    (n2,n1)=a.shape
    nn=np.zeros(2,int)
    nn[0]=n1
    nn[1]=n2
    fid = open(filename, 'wb')
    nn.tofile(fid)
    a.tofile(fid)
    fid.close()


def read_arr2_flt(filename):
    fid = open(filename, 'rb')
    nn=np.fromfile(fid,int,2)
    n1=nn[0]
    n2=nn[1]
    a = np.fromfile(fid,np.float32,n1*n2)
    return np.reshape(a, (n2, n1))



def read_arr2_cmpl(filename):
    fid = open(filename, 'rb')
    nn=np.fromfile(fid,int,2)
    n1=nn[0]
    n2=nn[1]
    a = np.fromfile(fid,np.complex64,n1*n2)
    return np.reshape(a, (n2, n1))


 
