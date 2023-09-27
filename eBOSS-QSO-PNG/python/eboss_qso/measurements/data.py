
from . import data_dir, DATA_VERSIONS
import os
from nbodykit.lab import FITSCatalog

def read_data(sample, version, focal_weights=False, NN_weights=False):
    """
    Read a eBOSS QSO data file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    focal_weights : bool, optional
        whether we are using focal plane redshift error corrections
    """
    assert version in DATA_VERSIONS
    
    if version == 'dr14':
        dr = 'DR14'
    elif version == 'dr16':
        dr = 'DR16'

    # get the file path
    filename = f'eBOSS_QSO_clustering_data-{sample}GC-v{dr}.fits'
    if NN_weights:
        filename = f'eBOSS_QSO_{sample}GC_v7_2.dat.fits.gz'
        version = 'dr16'
    path = os.path.join(data_dir, 'data', version, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    if focal_weights and version == 'dr14': usecols += ['WEIGHT_FOCAL']
    return FITSCatalog(path)[usecols]


def read_randoms(sample, version, NN_weights=True):
    """
    Read a eBOSS QSO randoms file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    """
    assert version in DATA_VERSIONS
    
    if version == 'dr14':
        dr = 'DR14'
    elif version == 'dr16':
        dr = 'DR16'

    # get the file path
    filename = f'eBOSS_QSO_clustering_random-{sample}GC-v{dr}.fits'
    if NN_weights:
        filename = f'eBOSS_QSO_{sample}GC_v7_2.ran.fits.gz'
        version = 'dr16'
    path = os.path.join(data_dir, 'data', version, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    if version == 'dr16': usecols += ['WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]

def finalize_data(s, cosmo, version, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', 'Weight', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    cosmo : Cosmology
        the cosmology parameters
    version : str
        the string specifying which version is in use
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """
    assert version in DATA_VERSIONS
    from nbodykit.transform import SkyToCartesian

    # add the Position column
    s['Position'] = SkyToCartesian(s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)

    # add systematic weights
    if version == 'dr16':
        s['Weight'] = s['WEIGHT_SYSTOT'] * s['WEIGHT_CP'] * s['WEIGHT_NOZ']
    if 'WEIGHT_NOZ' in s and 'WEIGHT_FOCAL' not in s and version == 'dr14':
        s['Weight'] = s['WEIGHT_SYSTOT'] * (s['WEIGHT_NOZ'] + s['WEIGHT_CP'] - 1.)
    elif 'WEIGHT_FOCAL' in s and version == 'dr14':
        s['Weight'] = s['WEIGHT_SYSTOT'] * s['WEIGHT_CP'] * s['WEIGHT_FOCAL']

    # FKP WEIGHT
    if P0_FKP is not None:
        s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
