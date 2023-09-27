import numpy as np
import eboss_qso.measurements as eboss
from eboss_qso import EBOSS_FITS
import os

import argparse

k_file = os.path.join(EBOSS_FITS, 'input', 'models', 'k_eff_old.dat')

k = np.loadtxt(k_file)

def main(ns):
    
    if ns.p == 0:
        p = None
    else:
        p = ns.p
    
    data = eboss.load_data_spectra('dr16', ns.sample, p=p, ell=0, focal_weights=False, NN_weights=ns.NN_weights)
    pk = data.poles['power_0'].real
    lenk = len(pk)
    k_eff = k[:lenk]
    
    print(pk[0])
    
    data_f = 'input_NN_weights' if ns.NN_weights else 'input'
    
    folder_data = os.path.join(EBOSS_FITS, data_f, 'data', 'spectra')
    sample = 'NGC' if ns.sample == 'N' else 'SGC'
    folder_sample = os.path.join(folder_data, sample)
    folder_p = os.path.join(folder_sample, f'{ns.p:2.1f}')
    
    pk_file = os.path.join(folder_p, 'Pk_data.dat')
    header = f"Pk for p={ns.p}, {sample} \n k(h/Mpc) Pk"
    
    if ns.norm:
        A = data.attrs['randoms.norm']
        #window_f = 'subsample_mean_NN_weights' if ns.NN_weights else 'subsample_mean'
        window_f = 'subsample_mean'
        folder_window = os.path.join(EBOSS_FITS, 'input', 'data', 'window', 'dr16', window_f)
        window = os.path.join(folder_window, f'poles_dr16-QSO-{ns.sample}-p_{ns.p:2.1f}.dat')
        Q = np.loadtxt(window)
        Q00 = Q[0,1]
        pk *= A / Q00
        header += f"*A/Q0(0) {A:5.5f} / {Q00:5.5f}"
    
    np.savetxt(pk_file, np.stack([k_eff,pk]).transpose(), header=header)
    
    print(pk_file, 'ready!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write in a .dat file the P(k) and k_eff")
    
    h = 'the sample, either North or South'
    parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)
    
    h = 'the value of p to use'
    parser.add_argument('--p', type = float, help = h, choices=[0., 1., 1.6, 3.0], required = True)
    
    h = 'flag renormalize the Pk to Q0(0)'
    parser.add_argument('--norm', action = 'store_true')
    
    h = 'whether to use the catalogs with NN systot weights'
    parser.add_argument('--NN-weights', action='store_true', help=h)
    
    args = parser.parse_args()
    
    main(args)