from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss
import numpy as np

import os
import argparse

#routines
def select_subsample(insize, outsize):

    index = numpy.zeros(insize, dtype=bool)
    valid = numpy.random.choice(range(insize), replace=False, size=outsize)
    index[valid] = True
    return index

def A_16(data, p=0):
    if p == 0:
        w1 = w2 = 1.
    else:
        w1 = eboss.bias_weight(data['Z'], eboss.fidcosmo, ell=0)
        w2 = eboss.fnl_weight(data['Z'], p=p)
        
    return (data['NZ'] * data['FKPWeight']**2 * data['Weight']**2 * w1 * w2).sum().compute()

def A_14(data, p=0):
    if p == 0:
        w1 = w2 = 1.
    else:
        w1 = eboss.bias_weight(data['Z'], eboss.fidcosmo, ell=0)
        w2 = eboss.fnl_weight(data['Z'], p=p)
        
    return (data['NZ'] * data['FKPWeight']**2 * data['Weight'] * w1 * w2).sum().compute()

def print_As(data):
    print("p=0.0: A_16",  A_16(data, p=0))
    print("p=0.0: A_14",  A_14(data, p=0))
    
    print("p=1.0: A_16",  A_16(data, p=1.0))
    print("p=1.0: A_14",  A_14(data, p=1.0))
    
    print("p=1.6: A_16",  A_16(data, p=1.6))
    print("p=1.6: A_14",  A_14(data, p=1.6))
    
    print("p=3.0: A_16",  A_16(data, p=3.0))
    print("p=3.0: A_14",  A_14(data, p=3.0))

#NGC
#data
data_N = eboss.read_data('N', 'dr16')

eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("Data NGC")

print_As(data_N)

#randoms
data_N = eboss.read_randoms('N', 'dr16')

eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

#print("All randoms")
#
#print_As(data_N)
#
#ns = 2200000
#print(ns, "randoms: seed + i")
#
#for i in [1,2,3,4,5]:
#
#    numpy.random.seed(42*(data_N.comm.rank+1000+(i)))
#    valid = select_subsample(data_N.size, int(ns//data_N.comm.size))
#    data = data_N[valid]
#    
#    print("Sample", i)
#    print_As(data)
#
#
#ns = 2200000
#print(ns, "randoms: seed + 2 * i")
#
#for i in [1,2,3,4,5]:
#
#    numpy.random.seed(42*(data_N.comm.rank+1000+(i*2)))
#    valid = select_subsample(data_N.size, int(ns//data_N.comm.size))
#    data = data_N[valid]
#    
#    print("Sample", i)
#    print_As(data)
    
#SGC
#data
data_S = eboss.read_data('S', 'dr16')

eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("Data SGC")

print_As(data_S)

#randoms
#data_S = eboss.read_randoms('S', 'dr16')
#
#eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)
#
#print("All randoms")
#
#print_As(data_S)
#
#ns = 1440000
#print(ns, "randoms: seed + i")
#
#for i in [1,2,3,4,5]:
#
#    numpy.random.seed(42*(data_S.comm.rank+1000+(i)))
#    valid = select_subsample(data_S.size, int(ns//data_S.comm.size))
#    data = data_S[valid]
#    
#    print("Sample", i)
#    print_As(data)
#
#
#ns = 1440000
#print(ns, "randoms: seed + 2 * i")
#
#for i in [1,2,3,4,5]:
#
#    numpy.random.seed(42*(data_S.comm.rank+1000+(i*2)))
#    valid = select_subsample(data_S.size, int(ns//data_S.comm.size))
#    data = data_S[valid]
#    
#    print("Sample", i)
#    print_As(data)