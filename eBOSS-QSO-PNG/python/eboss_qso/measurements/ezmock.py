
from . import data_dir, MOCK_VERSIONS, MOCK_SUBVERSIONS
import os
from nbodykit.lab import FITSCatalog


def read_ezmock_data(box, sample):
    """
    Read an eBOSS QSO EZ mock data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    """

    # get the file path
    filename = f'EZmock_realistic_eBOSS_QSO_{sample}GC_v7_{box:04d}.dat.fits'
    path = os.path.join(data_dir, 'ezmock', filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]


def read_new_ezmock_randoms(box, sample):
    """
    Read an eBOSS QSO EZ mock data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    """

    # get the file path
    filename = f'EZmock_realistic_eBOSS_QSO_{sample}GC_v7_{box:04d}.ran.fits'
    path = os.path.join(data_dir, 'ezmock', filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]

def read_ezmock_complete_data(box, sample):
    """
    Read an eBOSS QSO EZ mock complete data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    """

    # get the file path
    filename = f'EZmock_complete_eBOSS_QSO_{sample}GC_{box:04d}.dat.fits'
    path = os.path.join(data_dir, 'ezmock', filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    return FITSCatalog(path)[usecols]


def read_new_ezmock_complete_randoms(box, sample, shuffle=False):
    """
    Read an eBOSS QSO EZ mock complete data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    shuffle : bool
        if the shuffled randoms are readed
    """

    # get the file path
    if shuffle:
        filename = f'EZmock_complete_eBOSS_QSO_{sample}GC_{box:04d}.ran.shuf.fits'
    else:
        filename = f'EZmock_complete_eBOSS_QSO_{sample}GC.ran.fits'
        
    path = os.path.join(data_dir, 'ezmock', filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    return FITSCatalog(path)[usecols]


def read_ezmock_randoms(sample, version):
    """
    Read a eBOSS QSO EZ mock randoms file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    version : str
        the string specifying which version to load
    """
    assert version in MOCK_VERSIONS

    sample = 'ngc' if sample == 'N' else 'sgc'

    # get the file path
    filename = f'randx100_QSO_{version}_veto_{sample}.dat'
    path = os.path.join(data_dir, 'mocks', 'ezmock', 'randoms', filename)

    # load the source
    names = ['RA', 'DEC', 'Z', 'WEIGHT_FKP', 'COMP', 'NZ', 'VETO']
    s = CSVCatalog(path, names=names)

    # inverse completeness
    s['INV_COMP'] = 1. / s['COMP']

    # up-weight expected angular completeness
    s['NZ'] *= s['INV_COMP']

    return s


def finalize_ezmock(s, cosmo, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', 'Weight', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    cosmo : Cosmology
        the cosmology parameters
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """
    from nbodykit.transform import SkyToCartesian

    # add the Position column
    s['Position'] = SkyToCartesian(
        s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)

    # add systematic weights
    s['Weight'] = s['WEIGHT_SYSTOT'] * s['WEIGHT_CP'] * s['WEIGHT_NOZ']
    
    # FKP WEIGHT
    if P0_FKP is not None:
        s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
        
def finalize_ezmock_complete(s, cosmo, P0_FKP=3e4):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    cosmo : Cosmology
        the cosmology parameters
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """
    from nbodykit.transform import SkyToCartesian

    # add the Position column
    s['Position'] = SkyToCartesian(
        s['RA'], s['DEC'], s['Z'], cosmo, degrees=True)
    
    s['Weight'] = 1.
    
    # FKP WEIGHT
    s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
        