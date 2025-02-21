from nbodykit.lab import *
import eboss_qso.measurements as eboss
from eboss_qso.measurements.weights import fnl_weight, bias_weight
import numpy as np

def compute_effective_redshift_npower(r, cosmo, npow=1, p=None, P0=3e4, ell=0):
    """
    Compute effective redshift for n^(npow-1)
    """
    assert 'Z' in r.columns
    assert 'NZ' in r.columns

    # weights
    w_fkp = 1. / (1. + r['NZ']*P0)
    w_comp = r['Weight']
    if p is not None:
        w1 = fnl_weight(r['Z'], p=p)
        w2 = bias_weight(r['Z'], cosmo, ell=ell, eva=True)
    else:
        w1 = w2 = 1.0

    # effective redshift
    A = (r['Z'] * r['NZ']**(npow-1) * w_fkp**npow * w_comp**npow).sum().compute()
    norm = (r['NZ']**(npow-1) * w_fkp**npow * w_comp**npow).sum().compute()
    z_eff = A/norm

    return z_eff

def bias_model(z):
    alpha = 0.278
    beta = 2.393
    gamma = 6.565

    b = alpha * ((1 + z)**2 - gamma) + beta
    return b

def main():
    file_output = eboss.results_dir+'/effective_redshifts/effective_redshift_npower.txt'
    P0= 3e4

    with open(file_output, 'w') as f:
        print(f"P0 = {P0:.2e}", file=f)
        
        data_N = eboss.read_data('N', 'dr16')
        data_S = eboss.read_data('S', 'dr16')
        
        data_N = eboss.trim_redshift_range(data_N, zmin=0.8, zmax=2.2)
        data_S = eboss.trim_redshift_range(data_S, zmin=0.8, zmax=2.2)

        eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=P0)
        eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=P0)

        print("DATA: NGC \n", file=f)
        zeff = compute_effective_redshift_npower(data_N, eboss.fidcosmo, npow=2, P0=P0)
        print(f"zeff propto n^2: {zeff}, bias: {bias_model(zeff)}", file=f)
        zeff = compute_effective_redshift_npower(data_N, eboss.fidcosmo, npow=3, P0=P0)
        print(f"zeff propto n^3: {zeff}, bias: {bias_model(zeff)} \n", file=f)

        print("DATA: SGC \n", file=f)
        zeff = compute_effective_redshift_npower(data_S, eboss.fidcosmo, npow=2, P0=P0)
        print(f"zeff propto n^2: {zeff}, bias: {bias_model(zeff)}", file=f)
        zeff = compute_effective_redshift_npower(data_S, eboss.fidcosmo, npow=3, P0=P0)
        print(f"zeff propto n^3: {zeff}, bias: {bias_model(zeff)} \n", file=f)

if __name__ == '__main__':
    main()
