import numpy as np
import eboss_qso.measurements as eboss
from eboss_qso import EBOSS_FITS
import os

import argparse

def main(ns):
    
    pk_mocks = np.zeros(2000)
    pk_mocks = pk_mocks.reshape(2, 1000)
    lenpk1 = 1000
    
    mock_box = np.arange(1, ns.N + 1, 1)
    
    if ns.norm:
        folder_window = os.path.join(EBOSS_FITS, 'input', 'data', 'window', 'dr16', 'subsample_mean')
        window = os.path.join(folder_window, f'poles_dr16-QSO-{ns.sample}-p_{ns.p:2.1f}.dat')
        Q = np.loadtxt(window)
        Q00 = Q[0,1]
    else:
        Q00 = 1.
    
    if ns.p == 0:
        p = None
    else:
        p = ns.p
    
    for box in mock_box:
        mock = eboss.load_ezmock_spectra(ns.sample, p=p, box=box, subtract_shot_noise=True, average=False, ell=0)
        
        pk = mock.poles['power_0'].real
        lenpk2 = len(pk)
        lenpk = min(lenpk2, lenpk1)
        lenpk1 = lenpk
        #print(lenpk)
        
        if ns.norm:
            A = mock.attrs['randoms.norm']
        else:
            A = 1.
        pk *= A / Q00
        
        pk_mocks = np.vstack((pk_mocks[:,:lenpk], pk[:lenpk]))
        
    print(pk_mocks.shape)
    
    pk_mocks = pk_mocks[2:,:]
    
    print(pk_mocks.shape)
    
    cov = np.cov(pk_mocks.transpose())
    
    print(cov.shape)
    
    folder_data = os.path.join(EBOSS_FITS, 'input', 'data', 'spectra')
    sample = 'NGC' if ns.sample == 'N' else 'SGC'
    folder_sample = os.path.join(folder_data, sample)
    folder_p = os.path.join(folder_sample, f'{ns.p:2.1f}')
    
    mock_file = os.path.join(folder_p, 'Pk_mocks.dat')
    header = f"Pk of mocks for p={ns.p}, {sample}"
    if ns.norm:
        header += " Pk*A/Q0(0)"
    np.savetxt(mock_file, pk_mocks.transpose(), header=header)
    
    cov_file = os.path.join(folder_p, 'covariance.dat')
    header = f"Pk covarinace for p={ns.p}, {sample}"
    if ns.norm:
        header += " Pk*A/Q0(0)"
    np.savetxt(cov_file, cov, header=header)
    
    print(cov_file, 'ready!')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write in a .dat file the P(k) and k_eff")
    
    h = 'the sample, either North or South'
    parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)
    
    h = 'the value of p to use'
    parser.add_argument('--p', type = float, help = h, choices=[0., 1., 1.6, 3.0], required = True)
    
    h = 'the total number of mocks'
    parser.add_argument('--N', type = int, help = h, default = 1000)
    
    h = 'flag renormalize the Pk to Q0(0)'
    parser.add_argument('--norm', action = 'store_true')
    
    args = parser.parse_args()
    
    main(args)

