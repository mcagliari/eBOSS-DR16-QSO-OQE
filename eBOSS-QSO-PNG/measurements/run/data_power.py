
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss

import os
import argparse
import numpy as np

setup_logging()


def main(ns):

    # load the data first
    data = eboss.read_data(ns.sample, ns.version,
                           focal_weights=ns.focal_weights, NN_weights=ns.NN_weights)

    # load the randoms
    randoms = eboss.read_randoms(ns.sample, ns.version, NN_weights=ns.NN_weights)

    print(ns.zbins)
    # compute for every redshift bin
    for (zmin, zmax) in ns.zbins:

        # trim redshift range
        d = eboss.trim_redshift_range(data, zmin=zmin, zmax=zmax)
        r = eboss.trim_redshift_range(randoms, zmin=zmin, zmax=zmax)
        
        # finalize columns
        eboss.finalize_data(d, eboss.fidcosmo, ns.version, P0_FKP=ns.P0_FKP)
        eboss.finalize_data(r, eboss.fidcosmo, ns.version, P0_FKP=ns.P0_FKP)

        # compute effective values
        z_eff = eboss.compute_effective_redshift(r)
        nbar_eff = eboss.compute_effective_nbar(r)

        # combine data and randoms into the FKP source
        fkp = FKPCatalog(data=d, randoms=r, BoxPad=0.1)

        # the mesh kwargs to use
        mesh_kwargs = {'Nmesh': 512, 'interlaced': True,
                       'window': 'tsc', 'dtype': 'f8'}

        # compute unweighted results
        if ns.p == 0.:
            unweighted_mesh = fkp.to_mesh(
                nbar='NZ', fkp_weight='FKPWeight', comp_weight='Weight', **mesh_kwargs)

            # run
            result = ConvolvedFFTPower(first=unweighted_mesh, poles=[
                                       0, 2], dk=0.005, kmin=0.)

            # add effective redshift and nbar from randoms
            result.attrs['z_eff'] = z_eff
            result.attrs['nbar_eff'] = nbar_eff

            meta = {'p': None, 'zmin': zmin, 'zmax': zmax, 'P0_FKP': ns.P0_FKP, 'NN weights': ns.NN_weights}
          
            eboss.save_data_spectra(
                result, ns.sample, ns.version, ns.focal_weights, ns.NN_weights, **meta)
        else:

            # the fnl weight for the second field
            fkp['data/FnlWeight'] = d['FKPWeight'] * \
                eboss.fnl_weight(d['Z'], p=ns.p)
            fkp['randoms/FnlWeight'] = r['FKPWeight'] * \
                eboss.fnl_weight(r['Z'], p=ns.p)

            mesh2 = fkp.to_mesh(
                nbar='NZ', fkp_weight='FnlWeight', comp_weight='Weight', **mesh_kwargs)

            # use separate bias weights for ell=0,2
            for ell in [0,]: #we are interested in the monopole

                # the bias weight for the first field
                fkp['data/BiasWeight'] = d['FKPWeight'] * \
                    eboss.bias_weight(d['Z'], eboss.fidcosmo, ell=ell)
                fkp['randoms/BiasWeight'] = r['FKPWeight'] * \
                    eboss.bias_weight(r['Z'], eboss.fidcosmo, ell=ell)

                # convert to mesh
                mesh1 = fkp.to_mesh(
                    nbar='NZ', fkp_weight='BiasWeight', comp_weight='Weight', **mesh_kwargs)

                # run
                result = ConvolvedFFTPower(first=mesh1,
                                           second=mesh2,
                                           poles=[ell],
                                           dk=0.005,
                                           kmin=0.)

                # add effective redshift and nbar from randoms
                result.attrs['z_eff'] = z_eff
                result.attrs['nbar_eff'] = nbar_eff

                # and save
                meta = {'p': ns.p, 'zmin': zmin,
                        'zmax': zmax, 'P0_FKP': ns.P0_FKP}
                if ns.NN_weights:
                    meta['NN weights'] = True
                eboss.save_data_spectra(
                    result, ns.sample, ns.version, ns.focal_weights, ns.NN_weights, **meta)


if __name__ == '__main__':

    desc = 'compute the redshift weighted power spectra of the eBOSS QSO data sample'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str,
                       choices=['N', 'S'], help=h, required=True)

    h = 'the version to load'
    group.add_argument('--version', type=str,
                       choices=eboss.DATA_VERSIONS, help=h, required=True)

    h = 'the redshift bins to compute in the form of ZMIN,ZMAX'
    group.add_argument(
        '--zbins', type=eboss.redshift_range_type, nargs='+', help=h)

    h = 'the value of p to use'
    group.add_argument('--p', type=float, help=h,
                       choices=[0., 1., 1.6, 3.0], required=True)

    h = 'the P0 FKP version to use'
    parser.add_argument('--P0_FKP', type=float, default=3e4, help=h)

    h = 'whether to use focal plane weights'
    parser.add_argument('--focal-weights', action='store_true', help=h)
    
    h = 'whether to use the catalogs with NN systot weights'
    parser.add_argument('--NN-weights', action='store_true', help=h)

    # and go!
    main(parser.parse_args())
