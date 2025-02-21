from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss
import numpy as np

import os
import argparse

#routines
def I22(data):    
    return (data['NZ'] * data['FKPWeight']**2 * data['Weight']**2).sum().compute()

def I33(data):
    return (data['NZ']**2 * data['FKPWeight']**3 * data['Weight']**3).sum().compute()

def print_Is(data, alpha=1):
    print("I22",  I22(data) * alpha)
    print("I33",  I33(data) * alpha)
    

#NGC
#data
data_N = eboss.read_data('N', 'dr16')
eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

tot_data = data_N['Weight'].sum().compute()

print("Data NGC")

print_Is(data_N)

#randoms
data_N = eboss.read_randoms('N', 'dr16', NN_weights=False)
eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("Randoms NGC")
tot_random = data_N['Weight'].sum().compute()

print_Is(data_N, tot_data/tot_random)
    
#SGC
#data
data_S = eboss.read_data('S', 'dr16')
eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

tot_data = data_S['Weight'].sum().compute()

print("Data SGC")

print_Is(data_S)

#randoms
data_S = eboss.read_randoms('S', 'dr16', NN_weights=False)
eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

print("Randoms SGC")
tot_random = data_S['Weight'].sum().compute()

print_Is(data_S, tot_data/tot_random)
