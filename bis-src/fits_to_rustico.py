import numpy as np
from astropy.io import fits
import argparse
import os

#from eboss_qso.measurements import data_dir

def to_rustico(data, P0=3e4):
    data['WEIGHT_FKP'] = 1. / (1 + data['NZ']*P0)
    for_rustico = np.vstack([data['RA'], 
                             data['DEC'],
                             data['Z'],
                             data['WEIGHt_FKP'],
                             data['WEIGHT_CP'],
                             data['WEIGHT_NOZ'],
                             data['WEIGHT_SYSTOT'],
                             data['NZ']]).T
    
    return for_rustico

def main(ns):
    
    filename = f'eBOSS_QSO_clustering_{ns.type}-{ns.sample}GC-vDR16.fits'
    outname = f'eBOSS_QSO_clustering_{ns.type}-{ns.sample}GC-vDR16.txt'
    if ns.NN_weights:
        ctype = 'dat' if ns.type == 'data' else 'ran'
        filename = f'eBOSS_QSO_{ns.sample}GC_v7_2.dat.fits.gz'
        outname = f'eBOSS_QSO_{ns.sample}GC_v7_2.dat.txt'
    version = 'dr16'
    path_base = "path/to/catalogs"
    pathin = os.path.join(path_base, 'data', version, filename)
    pathout = os.path.join(path_base, 'data', outname)
    
    hdul = fits.open(pathin)
    data = hdul[1].data
    
    rustico = to_rustico(data, ns.P0)
    
    np.savetxt(pathout, rustico)
    
    print(ns)
    print("Origninal catalogue:", pathin)
    print(pathout, "ready!")
    
if __name__ == '__main__':
    
    desc = 'Produces a .txt version of eBOSS DR16 catalogs for Rustico inut'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str,
                       choices=['N', 'S'], help=h, required=True)

    h = 'catalog type'
    group.add_argument('--type', type=str, 
                       choices=['data', 'random'], help=h, required=True)
    
    h = 'the P0 FKP version to use'
    parser.add_argument('--P0', type=float, default=3e4, help=h)
    
    h = 'whether to use the catalogs with NN systot weights'
    parser.add_argument('--NN-weights', action='store_true', help=h)

    # and go!
    main(parser.parse_args())
    
    
