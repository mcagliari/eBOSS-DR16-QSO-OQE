import os
from nbodykit.cosmology import Cosmology, Planck15
from astropy.cosmology import Planck18_arXiv_v2
import socket

# data and result directories
data_dir = os.environ['EBOSS_DIR']
results_dir = os.path.join(data_dir, 'measurements')

# allowed data versions
DATA_VERSIONS = ['dr14', 'dr16']

# mock subversions
MOCK_VERSIONS = ['v1.8e']
MOCK_SUBVERSIONS = ['reg', 'no', 'fph']

# the fiducial DR12 cosmology
# _h = 0.676
# _Ob0 = 0.022/_h**2
# _Ocdm0 = 0.31 - _Ob0
# fidcosmo = Cosmology(h=_h, Omega_b=_Ob0, Omega_cdm=_Ocdm0, m_ncdm=None, T_cmb=2.7255)

#fidcosmo = Planck15
#DR16 fiducial cosmology Planck18
_kws = {'ln10^{10}A_s':3.043, 'n_s':0.9652, 'k_pivot':0.05, 'tau_reio':0.0540}
fidcosmo = Cosmology.from_astropy(Planck18_arXiv_v2, **_kws)

_h = 0.6777
_Ob0 = 0.048206
_Ocdm0 = 0.307115 - _Ob0
ezmock_cosmo = Cosmology(h=_h, Omega_b=_Ob0, Omega_cdm=_Ocdm0,
                         m_ncdm=None, n_s=0.9611, T_cmb=2.7255).match(sigma8=0.8225)

_h = 0.676
_Ob0 = 0.022/_h**2
_Ocdm0 = 0.31 - _Ob0
qpm_cosmo = Cosmology(h=_h, Omega_b=_Ob0, Omega_cdm=_Ocdm0,
                      m_ncdm=None, n_s=0.97, T_cmb=2.7255).match(sigma8=0.8)

# data functions
from .data import read_data, read_randoms, finalize_data
from .results import save_data_spectra, load_data_spectra

# ezmock functions
from .ezmock import read_ezmock_data, read_ezmock_randoms, read_new_ezmock_randoms, finalize_ezmock
from .ezmock import read_ezmock_complete_data, read_new_ezmock_complete_randoms, finalize_ezmock_complete
from .results import save_ezmock_spectra, load_ezmock_spectra

# qpm function
from .qpm import read_qpm_data, read_qpm_randoms, finalize_qpm
from .results import save_qpm_spectra, load_qpm_spectra

# results
from .results import save_RR_paircount, save_RR_poles

# utilities and weights
from .utils import trim_redshift_range, redshift_range_type, get_hashkeys, \
    compute_effective_redshift, compute_effective_nbar, \
    find_window_measurement, compute_effective_quantities


from .weights import fnl_weight, bias_weight, bias_model
