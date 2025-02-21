import os
from nbodykit.lab import FITSCatalog


##### Data ######


def read_data(sample, folder):
    """
    Read a eBOSS QSO data file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    folder : str
        the folder that contains the data file
    """
    

    # get the file path
    filename = f'eBOSS_QSO_clustering_data-{sample}GC-vDR16.fits'
    path = os.path.join(folder, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]

def read_randoms(sample, folder):
    """
    Read a eBOSS QSO randoms file.

    Parameters
    ----------
    sample : 'N' or 'S'
        the sample to load
    folder : str
        the folder that contains the random file
    """

    # get the file path
    filename = f'eBOSS_QSO_clustering_random-{sample}GC-vDR16.fits'
    path = os.path.join(folder, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]

def finalize_data(s, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', 'Weight', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """
    
    # add systematic weights
    s['Weight'] = s['WEIGHT_SYSTOT'] * s['WEIGHT_CP'] * s['WEIGHT_NOZ']
    

    # FKP WEIGHT
    if P0_FKP is not None:
        s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)

        
###### EZmocks ######
        
        
def read_ezmock_data(box, sample, folder):
    """
    Read an eBOSS QSO realistic EZmock data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    folder : str
        the folder that contains the data file
    """

    # get the file path
    filename = f'EZmock_realistic_eBOSS_QSO_{sample}GC_v7_{box:04d}.dat.fits'
    path = os.path.join(folder, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]

def read_ezmock_randoms(box, sample, folder):
    """
    Read an eBOSS QSO realistic EZmock data file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    folder : str
        the folder that contains the random file
    """

    # get the file path
    filename = f'EZmock_realistic_eBOSS_QSO_{sample}GC_v7_{box:04d}.ran.fits'
    path = os.path.join(folder, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT']
    return FITSCatalog(path)[usecols]

def read_ezmock_complete_data(box, sample, folder):
    """
    Read an eBOSS QSO complete EZmock data file.
    NOTE: this data file is shared between the complete and shuffled EZmock datasets.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    folder : str
        the folder that contains the data file
    """

    # get the file path
    filename = f'EZmock_complete_eBOSS_QSO_{sample}GC_{box:04d}.dat.fits'
    path = os.path.join(folder, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    return FITSCatalog(path)[usecols]


def read_ezmock_complete_randoms(box, sample, folder, shuffle=False):
    """
    Read an eBOSS QSO complete or shuffled EZmock random file.

    Parameters
    ----------
    box : int
        the box number to load
    sample : 'N' or 'S'
        the sample to load
    folder : str
        the folder that contains the random file
    shuffle : bool
        if the shuffled randoms are readed
    """

    # get the file path
    if shuffle:
        filename = f'EZmock_complete_eBOSS_QSO_{sample}GC_{box:04d}.ran.shuf.fits'
    else:
        filename = f'EZmock_complete_eBOSS_QSO_{sample}GC.ran.fits'
        
    path = os.path.join(folder, filename)

    # load the source
    usecols = ['RA', 'DEC', 'Z', 'NZ']
    return FITSCatalog(path)[usecols]

def finalize_ezmock(s, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', 'Weight', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """

    # add systematic weights
    s['Weight'] = s['WEIGHT_SYSTOT'] * s['WEIGHT_CP'] * s['WEIGHT_NOZ']
    
    # FKP WEIGHT
    if P0_FKP is not None:
        s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
        
def finalize_ezmock_complete(s, P0_FKP=None):
    """
    Finalize the creation of a CatalogSource from a data file by
    adding 'Position', and 'FKPWeight'.

    Parameters
    ----------
    s : CatalogSource
        the catalog source object
    P0_FKP : float, optional
        the P0 value to use to for FKPWeights
    """
    from nbodykit.transform import SkyToCartesian
    
    # set systematic weights to 1
    s['Weight'] = 1.
    
    # FKP WEIGHT
    s['FKPWeight'] = 1. / (1 + s['NZ']*P0_FKP)
        
