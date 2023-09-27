
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss
from astropy.io import fits
import astropy.utils.data

import os
import argparse

setup_logging()


def main(ns):

    

    for box_num in range(ns.start, ns.stop, ns.step):

        urlpath = f'https://data.sdss.org/sas/dr17/eboss/lss/EZmocks/v1_0_0/realistic/eBOSS_QSO/'
        dataname = f'dat/EZmock_realistic_eBOSS_QSO_{ns.sample}GC_v7_{box_num:04d}.dat.fits.gz'
        randname = f'ran/EZmock_realistic_eBOSS_QSO_{ns.sample}GC_v7_{box_num:04d}.ran.fits.gz'
        
        data_fits = fits.open(urlpath+dataname, cache=False)
        data_path = os.path.join(eboss.data_dir, 'ezmock', dataname[4:-3])
        data_fits.writeto(data_path)
        del data_fits
        
        rand_fits = fits.open(urlpath+randname, cache=False)
        rand_path = os.path.join(eboss.data_dir, 'ezmock', randname[4:-3])
        rand_fits.writeto(rand_path)
        del rand_fits
        
        zmin, zmax = ns.zbins
            
        # load the data
        data = eboss.read_ezmock_data(box_num, ns.sample)
        data = eboss.trim_redshift_range(data, zmin=zmin, zmax=zmax)
        
        eboss.finalize_ezmock(data, eboss.ezmock_cosmo, P0_FKP=ns.P0_FKP)
        
        # load the randoms
        randoms = eboss.read_new_ezmock_randoms(box_num, ns.sample)
        randoms = eboss.trim_redshift_range(randoms, zmin=zmin, zmax=zmax)
        
        eboss.finalize_ezmock(randoms, eboss.ezmock_cosmo, P0_FKP=ns.P0_FKP)

        # combine data and randoms into the FKP source
        fkp = FKPCatalog(data=data, randoms=randoms, BoxPad=0.1)

        # mesh kwargs
        mesh_kwargs = {'Nmesh': 512, 'interlaced': True,
                           'window': 'tsc', 'dtype': 'f8'}

        # compute unweighted results
        if ns.p == 0.:
            unweighted_mesh = fkp.to_mesh(
                    nbar='NZ', fkp_weight='FKPWeight', comp_weight='Weight', **mesh_kwargs)

            # run
            result = ConvolvedFFTPower(first=unweighted_mesh, poles=[
                                           0, 2], dk=0.005, kmin=0.)

            # save
            meta = {'p': None, 'zmin': 0.8,
                        'zmax': 2.2, 'P0_FKP': ns.P0_FKP}
            eboss.save_ezmock_spectra(result, box_num, ns.sample, **meta)

        else:

            # the fnl weight for the second field
            fkp['data/FnlWeight'] = data['FKPWeight'] * \
                eboss.fnl_weight(data['Z'], p=ns.p)
            fkp['randoms/FnlWeight'] = randoms['FKPWeight'] * \
                eboss.fnl_weight(randoms['Z'], p=ns.p)

            # mesh for f_nl weight
            mesh2 = fkp.to_mesh(
                nbar='NZ', fkp_weight='FnlWeight', comp_weight='Weight', **mesh_kwargs)

            # weight ell=0,2, separately
            for ell in [0,]:

                # the bias weight for the first field
                fkp['data/BiasWeight'] = data['FKPWeight'] * \
                    eboss.bias_weight(
                        data['Z'], eboss.ezmock_cosmo, ell=ell)
                fkp['randoms/BiasWeight'] = randoms['FKPWeight'] * \
                    eboss.bias_weight(
                        randoms['Z'], eboss.ezmock_cosmo, ell=ell)

                # convert to mesh
                mesh1 = fkp.to_mesh(
                    nbar='NZ', fkp_weight='BiasWeight', comp_weight='Weight', **mesh_kwargs)

                # run
                result = ConvolvedFFTPower(first=mesh1,
                                           second=mesh2,
                                           poles=[ell],
                                           dk=0.005,
                                           kmin=0.)
                # save
                meta = {'p': ns.p, 'zmin': 0.8,
                        'zmax': 2.2, 'P0_FKP': ns.P0_FKP}
                eboss.save_ezmock_spectra(result, box_num, ns.sample, **meta)
        
        os.system('rm ' + data_path)
        os.system('rm ' + rand_path)


if __name__ == '__main__':

    desc = 'compute the redshift weighted power spectra of the eBOSS QSO EZ mocks'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')


    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str,
                       choices=['N', 'S'], help=h, required=True)
    
    h = 'the value of p to use'
    group.add_argument('--p', type=float, help=h,
                       choices=[0., 1., 1.6, 3.0], required=True)

    h = 'the P0 FKP version to use'
    parser.add_argument('--P0_FKP', type=float, default=3e4, help=h)
    
    h = 'the redshift bins to compute in the form of ZMIN,ZMAX'
    group.add_argument(
        '--zbins', type=eboss.redshift_range_type, nargs='+', help=h)

    h = 'the start box number'
    parser.add_argument('--start', default=1, type=int, help=h)

    h = 'the stop box number'
    parser.add_argument('--stop', default=1000, type=int, help=h)

    h = 'the step box number'
    parser.add_argument('--step', default=1, type=int, help=h)

    #astropy conf
    conf = astropy.utils.data.Conf()
    
    with conf.set_temp('remote_timeout', 90):
      # and go!
      main(parser.parse_args())
