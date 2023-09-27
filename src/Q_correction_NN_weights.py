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

def A_16(data):
    return (data['NZ'] * data['FKPWeight']**2 * data['Weight']**2).sum().compute()

def A_14(data):
    return (data['NZ'] * data['FKPWeight']**2 * data['Weight']).sum().compute()

def Q_correction(data):
    return ((data['Weight'].sum())**2 - (data['Weight']**2).sum()).compute()

def alpha(data, total):
    d = data['Weight'].sum().compute()
    return total/d

def print_results(sub_random, tot_data, tot_random):
    print("none", Q_correction(sub_random))

    print("A_16",  A_16(sub_random))
    print("A_14",  A_14(sub_random))
    print("alpha", alpha(sub_random, tot_data))
    print("alpha randoms", alpha(sub_random, tot_random))
    
#main
#NGC
#data
data_N = eboss.read_data('N', 'dr16', NN_weights=True)

data_N = eboss.trim_redshift_range(data_N, zmin=0.8, zmax=2.2)

eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("Data NGC")

print("none", Q_correction(data_N))

print("A_16",  A_16(data_N))
print("A_14",  A_14(data_N))
print(eboss.compute_effective_quantities(data_N, eboss.fidcosmo, P0=3e4))

tot_data = data_N['Weight'].sum().compute()

#randoms
data_N = eboss.read_randoms('N', 'dr16', NN_weights=True)

data_N = eboss.trim_redshift_range(data_N, zmin=0.8, zmax=2.2)

eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("All randoms")

#print("none", Q_correction(data_N))
#
#print("A_16",  A_16(data_N))
#print("A_14",  A_14(data_N))
#print(eboss.compute_effective_quantities(data_N, eboss.fidcosmo, P0=3e4))

tot = data_N['Weight'].sum().compute()

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
#    print_results(data, tot_data, tot)


ns = 2200000
print(ns, "randoms: seed + 2 * i")

for i in [1,2,3,4,5]:

    numpy.random.seed(42*(data_N.comm.rank+1000+(i*2)))
    valid = select_subsample(data_N.size, int(ns//data_N.comm.size))
    data = data_N[valid]
    
    print("Sample", i)
    print_results(data, tot_data, tot)
    
#SGC
#data
data_S = eboss.read_data('S', 'dr16', NN_weights=True)

data_S = eboss.trim_redshift_range(data_S, zmin=0.8, zmax=2.2)

eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("Data SGC")

print("none", Q_correction(data_S))

print("A_16",  A_16(data_S))
print("A_14",  A_14(data_S))
print(eboss.compute_effective_quantities(data_S, eboss.fidcosmo, P0=3e4))

tot_data = data_S['Weight'].sum().compute()

#randoms
data_S = eboss.read_randoms('S', 'dr16', NN_weights=True)

data_S = eboss.trim_redshift_range(data_S, zmin=0.8, zmax=2.2)

eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

#print("All randoms")
#
#print("none", Q_correction(data_S))
#
#print("A_16",  A_16(data_S))
#print("A_14",  A_14(data_S))
#print(eboss.compute_effective_quantities(data_S, eboss.fidcosmo, P0=3e4))

tot = data_S['Weight'].sum().compute()

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
#    print_results(data, tot_data, tot)
#

ns = 1440000
print(ns, "randoms: seed + 2 * i")

for i in [1,2,3,4,5]:

    numpy.random.seed(42*(data_S.comm.rank+1000+(i*2)))
    valid = select_subsample(data_S.size, int(ns//data_S.comm.size))
    data = data_S[valid]
    
    print("Sample", i)
    print_results(data, tot_data, tot)