#! python
########################################################################
# Shimon Panfil: Industrial Physics and Simulations                   ##
# http://industrialphys.com                                           ##
# THE SOFTWARE IS PROVIDED "AS IS", USE IT AT YOUR OWN RISK           ##
########################################################################
import numpy as np
import matplotlib.pyplot as plt

def read_arr2_flt(filename):
    fid = open(filename, 'rb')
    nn=np.fromfile(fid,np.int32,2)
    n1=nn[0]
    n2=nn[1]
    a = np.fromfile(fid,np.float32,n1*n2)
    return np.reshape(a, (n2, n1))

def read_arr2_dbl(filename):
    fid = open(filename, 'rb')
    nn=np.fromfile(fid,np.int32,2)
    n1=nn[0]
    n2=nn[1]
    a = np.fromfile(fid,np.float64,n1*n2)
    return np.reshape(a, (n2, n1))



lp=read_arr2_dbl('lp.dbl')
rk=read_arr2_dbl('rk.dbl')
t=lp[:,2]
erk=0.5*rk[:,1]**2-np.cos(rk[:,0])
elp=0.5*lp[:,1]**2-np.cos(lp[:,0])
plt.plot(t,elp,t,erk)
plt.show()


