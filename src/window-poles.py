import numpy
import os
import hashlib
import argparse
import eboss_qso.measurements as eboss

from eboss_qso.fits.window import compute_window
from eboss_qso.measurements import data_dir
from eboss_qso.measurements.utils import make_hash, get_hashkeys


def find_window_measurement(version, sample, zmin, zmax, p, ell, subsample, sub_number=None):
    """
    Try to find and return a matching window function file.

    Parameters
    ----------
    version : str
        the data version
    sample : 'N','S'
        the sample
    zmin : float
        the minimum redshift
    zmax : float
        the maximum redshift
    p : None, 1.0, 1.6
        FKP-only is p=None
    ell : 0, 2
        the multipole to load
    subsample : int
        number of object in the subsample
    sub_number : int
        iteration number of the paircount file
    """
    from glob import glob

    # the directory holding any window results
    home_dir = os.environ['EBOSS_DIR']
    dirname = os.path.join(data_dir, 'measurements', 'window', version)
    
    if sub_number is not None:
        dirname = os.path.join(dirname, 'subsample_iteration')

    filename = f"RR_eboss_{version}-QSO-{sample}-*.json"
    pattern = os.path.join(dirname, filename)

    # search all file matches
    for f in glob(pattern):
        hashinput = get_hashkeys(f, 'SurveyDataPairCount')

        # compare zmin and zmax
        if sub_number is None:
            x = [zmin, zmax, subsample]
            y = [hashinput[k] for k in ['zmin', 'zmax', 'subsample']]
        else:
            x = [zmin, zmax, subsample, sub_number]
            y = [hashinput[k] for k in ['zmin', 'zmax', 'subsample', 'subsample_number']]
        if numpy.allclose(x, y):
            if hashinput.get('p', -1) == p:
                if p is None:
                    return f
                elif hashinput.get('ell') == ell:
                    return f

    msg = f"no window file match found for pattern '{pattern}'; "
    msg += f" ell={ell}, p={p}"
    raise ValueError(msg)

def write_window(zmin, zmax, sample, version, p, ell, subsample, sub_number=None, NN_weights=False):
    """
    Write the necessary window file.
    """
    
    # the window file
    if sub_number is None:
        meta = {'zmin': zmin, 'zmax': zmax, 'p': p, 'ell': ell, 'subsample': subsample}
    else:
        meta = {'zmin': zmin, 'zmax': zmax, 'p': p, 'ell': ell, 'subsample': subsample, 'subsample_number': sub_number}
    if NN_weights:
        meta['NN weights'] = True
    hashstr = make_hash(meta)
    filename = f"poles_{version}-QSO-{sample}_{hashstr}.dat"
    output_dir = os.path.join(data_dir, 'fits/input/data/window/', version)
    output = os.path.join(output_dir, filename)
    
   
    # make the window
   
    # the name of the (unformatted) window file
    if p == -1.: p = None
    
    window_file = find_window_measurement(version, sample,
                                          zmin, zmax, p, ell, subsample, sub_number)
    print('using window file %s...' % window_file)
    ells = [0, 2, 4, 6, 8]
    compute_window(window_file, ells, output, meta, smin=1e-2,
                   smax=1e4, quiet=False)
    
def write_named_window(zmin, zmax, sample, version, p, ell, subsample, name):
    """
    Write the necessary window file.
    """
    
    # the window file
    meta = {'zmin': zmin, 'zmax': zmax, 'p': p, 'ell': ell, 'subsample': subsample, 'name': name}
    filename = f"poles_{version}-QSO-{sample}_{name}.dat"
    output_dir = os.path.join(data_dir, 'fits/input/data/window/', version)
    output = os.path.join(output_dir, filename)
    
   
    # make the window
   
    # the name of the (unformatted) window file
    
    window_file = f"RR_eboss_{version}-QSO-{sample}-{name}.json"
    window_dir = os.path.join(data_dir, 'measurements/window/', version)
    window_file = os.path.join(window_dir, window_file)
    print('using window file %s...' % window_file)
    ells = [0, 2, 4, 6, 8]
    compute_window(window_file, ells, output, meta, smin=1e-2,
                   smax=1e4, quiet=False)
    
def main(ns):
    
    if ns.name is None:
        if ns.subnumber is None:
            write_window(ns.zmin, ns.zmax, ns.sample, ns.version, ns.p, 0, ns.subsample, NN_weights=ns.NN_weights)
        else:
            for i in range(1, ns.subnumber+1):
                write_window(ns.zmin, ns.zmax, ns.sample, ns.version, ns.p, 0, ns.subsample, i, ns.NN_weights)
    else:
        write_named_window(ns.zmin, ns.zmax, ns.sample, ns.version, ns.p, 0, ns.subsample, ns.name)
    
if __name__ == '__main__':
        
    parser = argparse.ArgumentParser(description="Compute the poles of a given DD window function")
    
    h = 'the sample, either North or South'
    parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)
    
    h = 'the version to load'
    parser.add_argument('--version', type = str, choices = eboss.DATA_VERSIONS, help = h, default = 'dr16')
    
    h = 'minimum redshift'
    parser.add_argument('--zmin', type = float, help = h, default = 0.8)
    
    h = 'maximum redshift'
    parser.add_argument('--zmax', type = float, help = h, default = 2.2)
    
    h = 'the value of p to use'
    parser.add_argument('--p', type = float, help = h, choices=[-1., 0., 1., 1.6, 3.0], required = True)
    
    h = 'subsample over which the window function was computed'
    parser.add_argument('--subsample', type = float, help = h, required = True)
    
    h = 'iteration number of the DD paircount file'
    parser.add_argument('--subnumber', type = int, help = h, default = None)
    
    h = 'identifying name on the DD paircount file if not hash named'
    parser.add_argument('--name', type = str, help = h, default = None)
    
    h = 'whether to use the paircounts from the catalogs with NN systot weights'
    parser.add_argument('--NN-weights', action='store_true', help=h)
    
    args = parser.parse_args()
    
    main(args)