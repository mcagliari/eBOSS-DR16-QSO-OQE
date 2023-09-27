import numpy as np
import eboss_qso.measurements as eboss
from eboss_qso import EBOSS_FITS
import os

import argparse

def main(ns):
    
    #container for Pks computed with mock+random mock
    pk_mocks = np.zeros(2000)
    pk_mocks = pk_mocks.reshape(2, 1000)
    #container for Pks computed with mock+random data
    pk_mocks_ric = np.zeros(2000)
    pk_mocks_ric = pk_mocks_ric.reshape(2, 1000)
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
        mock = eboss.load_ezmock_spectra(ns.sample, p=p, box=box, subtract_shot_noise=True, average=False, ell=0, eztype='complete')
        mock_ric = eboss.load_ezmock_spectra(ns.sample, p=p, box=box, subtract_shot_noise=True, average=False, ell=0, eztype='shuffle')
        
        pk = mock.poles['power_0'].real
        pk_ric = mock_ric.poles['power_0'].real
        lenpk2 = len(pk)
        lenpk3 = len(pk_ric)
        lenpk = min(lenpk3, lenpk2, lenpk1)
        lenpk1 = lenpk
        #print(lenpk)
        
        if ns.norm:
            A = mock.attrs['randoms.norm']
            A_ric = mock_ric.attrs['randoms.norm']
        else:
            A = 1.
            A_ric = 1.
        pk *= A / Q00
        pk_ric *= A_ric / Q00
        
        pk_mocks = np.vstack((pk_mocks[:,:lenpk], pk[:lenpk]))
        pk_mocks_ric = np.vstack((pk_mocks_ric[:,:lenpk], pk_ric[:lenpk]))
        
    print(pk_mocks.shape, pk_mocks_ric.shape)
    
    pk_mocks = pk_mocks[2:,:]
    pk_mocks_ric = pk_mocks_ric[2:,:]
    
    print(pk_mocks.shape, pk_mocks_ric.shape)
    
    mean_mock = np.mean(pk_mocks, axis=0)
    mean_mock_ric = np.mean(pk_mocks_ric, axis=0)
    
    Wkric = (mean_mock - mean_mock_ric) / mean_mock
    
    print(mean_mock.shape, mean_mock_ric.shape)
    
    folder_ric = os.path.join(EBOSS_FITS, 'input', 'data', 'window', 'Wkric')
    sample = 'NGC' if ns.sample == 'N' else 'SGC'
    folder_sample = os.path.join(folder_ric, sample)
    folder_p = os.path.join(folder_sample, f'{ns.p:2.1f}')
    
    ric_file = os.path.join(folder_p, 'Wkric.dat')
    mm_file = os.path.join(folder_p, 'mean_mock.dat')
    mmr_file = os.path.join(folder_p, 'mean_mock_ric.dat')
    
    header = f"W_ric(k) for p={ns.p}, {sample}"
    if ns.norm:
        header += " Pk*A/Q0(0)"
    np.savetxt(ric_file, Wkric, header=header)
    
    header = f"mean mock P(k) for p={ns.p}, {sample}"
    if ns.norm:
        header += " Pk*A/Q0(0)"
    np.savetxt(mm_file, mean_mock, header=header)
    
    header = f"mean mock RIC P(k) for p={ns.p}, {sample}"
    if ns.norm:
        header += " Pk*A/Q0(0)"
    np.savetxt(mmr_file, mean_mock_ric, header=header)
    
    print(ric_file, 'ready!')
    

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

