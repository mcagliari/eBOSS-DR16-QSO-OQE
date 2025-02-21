from nbodykit.lab import *
from nbodykit import setup_logging

import numpy as np

import os
import argparse

def get_Q00(sample):
    if sample == 'N':
        return 1.620437571133163335e+00
    elif sample == 'S':
        return 9.087155237345634351e-01

def get_alpha16(w_data, w_random):
    """
    Compute alpha_16 = (sum_data w_c_data) / (sum_random w_c_random).
    
    Parameters
    ----------
    w_data: array-like
        data completeness weights
    w_random: array-like
        random completeness weights
    """
    data = w_data.sum().compute()
    random = w_random.sum().compute()
    
    return data/random    
    
def get_beta(Q00, I22):
    """
    Compute beta = Q_0(0) / I22.
    """
    return Q00/I22
    
def get_Iab(catalog, a, b):
    """
    Compute Iab = sum_i nz_i^(a-1) * wc_i^b * wFKP_i^b
    
    Parameters
    ----------
    catalog : nbodykit FITSCatalog
        catalog for which we compute Iab
    """
    if a == 1:
        return (catalog['FKPWeight']**b * catalog['Weight']**b).sum().compute()
    
    return (catalog['NZ']**(a-1) * catalog['FKPWeight']**b * catalog['Weight']**b).sum().compute()
    
    
def compute_coefficients(sample, eboss_code=False, folder=None, mock=False, mock_type=None, box=None, folder_r=None, NN_weights=False):
    """
    Read the data and random catalog and compute [alpha_16, beta, I22_random, I33_random, I23_random, I13_data, I13_random, I12_data, I12_random] for a given pair of data and random catalogs
    
    Parameters
    ----------
    sample : str, 'N' or 'S'
        the sample to load
    eboss_code : bool, default=False
        if uses the loading routines from eBOSS-QSO-PNG (True) or from read_data_mocks.py (False) 
    folder : str, default=None
        folder containing the data catalog, not required if eboss_code=True
    mock : bool, default=False
        if the catalogs are the real data (False) or an EZmock (True)
    mock_type : str, 'realistic', 'complete', 'shuffle', default=None
        the EZmock type. Required if mock=True
    box : int, default=None
        EZmock box number. Required if mock=True
    folder_r : str
        folder containing the random catalog if different from the folder containig the data catalog
    """
    
    assert sample in ['N', 'S']
    if mock:
        assert mock_type is not None, "Please specify a mock type"
        assert mock_type in ['realistic', 'complete', 'shuffle']
        assert box is not None, "A mock box number is required"
    
    #load catalogs
    if eboss_code:
        import eboss_qso.measurements as eboss
        
        if not mock:
            #load data
            data = eboss.read_data(sample, 'dr16', NN_weights=NN_weights)
            eboss.finalize_data(data, eboss.fidcosmo, 'dr16', P0_FKP=3e4)
            
            #load randoms
            randoms = eboss.read_randoms(sample, 'dr16', NN_weights=NN_weights)
            eboss.finalize_data(randoms, eboss.fidcosmo, 'dr16', P0_FKP=3e4)
        else:
            if mock_type == 'realistic':
                # load data
                data = eboss.read_ezmock_data(box, sample)
                eboss.finalize_ezmock(data, eboss.ezmock_cosmo, P0_FKP=3e4)
                
                # load randoms
                randoms = eboss.read_new_ezmock_randoms(box, sample)
                eboss.finalize_ezmock(randoms, eboss.ezmock_cosmo, P0_FKP=3e4)

            else:
                shuffle = True if mock_type=='shuffle' else False
                # load data
                data = eboss.read_ezmock_complete_data(box, sample)
                eboss.finalize_ezmock_complete(data, eboss.ezmock_cosmo, P0_FKP=3e4)
                
                # load randoms
                randoms = eboss.read_new_ezmock_complete_randoms(box, sample, shuffle=shuffle)
                eboss.finalize_ezmock_complete(randoms, eboss.ezmock_cosmo, P0_FKP=3e4)            
    else:
        import read_data_mocks as rdm
        folder_r = folder if folder_r is None else folder_r
        
        if not mock:
            #load data
            data = rdm.read_data(sample, folder)
            rdm.finalize_data(data, P0_FKP=3e4)
            
            #load randoms
            randoms = rdm.read_randoms(sample, folder_r)
            rdm.finalize_data(randoms, P0_FKP=3e4)
        else:
            if mock_type == 'realistic':
                # load data
                data = rdm.read_ezmock_data(box, sample, folder)
                rdm.finalize_ezmock(data, P0_FKP=3e4)
                
                # load randoms
                randoms = rdm.read_ezmock_randoms(box, sample, folder_r)
                rdm.finalize_ezmock(randoms, P0_FKP=3e4)

            else:
                shuffle = True if mock_type=='shuffle' else False
                # load data
                data = rdm.read_ezmock_complete_data(box, sample, folder)
                rdm.finalize_ezmock_complete(data, P0_FKP=3e4)
                
                # load randoms
                randoms = rdm.read_ezmock_complete_randoms(box, sample, folder_r, shuffle=shuffle)
                rdm.finalize_ezmock_complete(randoms, P0_FKP=3e4)
    
    # Computing the coefficients
    alpha_16 = get_alpha16(data['Weight'], randoms['Weight'])
    
    Q00 = get_Q00(sample)
    I22_random = get_Iab(randoms, 2, 2)
    beta = get_beta(Q00, alpha_16 * I22_random)
    
    I33_random = get_Iab(randoms, 3, 3)
    
    I23_random = get_Iab(randoms, 2, 3)
    
    I13_data = get_Iab(data, 1, 3)
    I13_random = get_Iab(randoms, 1, 3)
    
    I12_data = get_Iab(data, 1, 2)
    I12_random = get_Iab(randoms, 1, 2)
    
    return [alpha_16, beta, I22_random, I33_random, I23_random, I13_data, I13_random, I12_data, I12_random]

def main(ns):
    
    alphi = []
    if not ns.mocks:
        name = 'data'
        folder = None if ns.folder is None else str(ns.folder)
        folder_r = None if ns.folder_r is None else str(ns.folder_r)
            
        alphi.append(compute_coefficients(ns.sample, 
                                          eboss_code=ns.eboss_code, 
                                          folder=folder, 
                                          folder_r=folder_r,
                                          NN_weights=ns.NN_weights))
    else:
        first_m = 1 if ns.m_start is None else ns.m_start
        last_m = 1001 if ns.m_end is None else ns.m_end + 1
        
        assert first_m >= 1, "--m-start must be >= 1"
        assert last_m <= 1001, "--m-end must be <= 1000"
        assert first_m < last_m, "--m-start must be < --m-end"
        
        name = 'EZmocks_' + str(ns.mock_type)
        for i in range(first_m, last_m):
            folder = None if ns.folder is None else str(ns.folder)
            folder_r = None if ns.folder_r is None else str(ns.folder_r)
            
            alphi.append(compute_coefficients(ns.sample, 
                                              eboss_code=ns.eboss_code, 
                                              folder=folder, 
                                              mock=True, 
                                              mock_type=str(ns.mock_type), 
                                              box=i, 
                                              folder_r=folder_r))
            
    alphi = np.array(alphi)
    file_name = f"Iabs_{ns.sample}GC_" + name 
    file_name += "_NN_weights.dat" if ns.NN_weights else ".dat"
    path = os.path.join(ns.output, file_name)
    headers = f'''Coeffiencts for {ns.sample}GC {name} NN weights: {ns.NN_weights} \n NOTE. Iab_r are not yet multiplied by alpha_16 \n alpha_16  beta  I22_r  I33_r  I23_r  I13_d  I13_r  I12_d  I12_r'''
    
    np.savetxt(path, alphi, header=headers)
    print(path, "ready!")
    
    
    
if __name__ == '__main__':
    desc = '''Compute and write on file the coefficients for the Rustico's spectra and bispectra post-processing'''
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str, choices=['N', 'S'], help=h, required=True)
    
    h = 'folder to save the coefficient file'
    group.add_argument('--output', help=h, required=True)

    
    h = 'flag for using mock catalogs'
    parser.add_argument('--mocks', action='store_true', help=h)
    
    h = 'first mock box to analyse, must be >= 1'
    parser.add_argument('--m-start', type=int, default=None, help=h)
    
    
    h = 'last mock box to analyse, must be <= 1000'
    parser.add_argument('--m-end', type=int, default=None, help=h)
    
    h = 'mock catalog type, realistic, complete, or shuffle'
    parser.add_argument('--mock-type', default=None, help=h)
    
    h = 'flag for using the reading routines in eBOSS-QSO-PNG'
    parser.add_argument('--eboss-code', action='store_true', help=h)
    
    h = 'folder containing the data catalogs'
    parser.add_argument('--folder', default=None, help=h)
    
    h = 'folder containing the random catalogs'
    parser.add_argument('--folder-r', default=None, help=h)
    
    h = 'use NN catalogflag'
    parser.add_argument('--NN-weights', action='store_true', help=h)
    
    # and go!
    main(parser.parse_args())
    
