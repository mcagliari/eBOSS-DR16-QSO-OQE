
from nbodykit.lab import *
from nbodykit import setup_logging
import eboss_qso.measurements as eboss
import os
import argparse

setup_logging()


def select_subsample(insize, outsize):

    index = numpy.zeros(insize, dtype=bool)
    valid = numpy.random.choice(range(insize), replace=False, size=outsize)
    index[valid] = True
    return index


def main(ns):

    # load the randoms
    randoms = eboss.read_randoms(ns.sample, ns.version, NN_weights=ns.NN_weights)

    # set the seed
    numpy.random.seed(42*(randoms.comm.rank+1000+(ns.nsub*2))) #original + ns.nsub

    # select a subsample?
    if ns.subsample is not None and randoms.csize > ns.subsample:
        valid = select_subsample(randoms.size, int(
            ns.subsample//randoms.comm.size))
        randoms = randoms[valid]

    # trim redshift range
    r = eboss.trim_redshift_range(randoms, zmin=ns.zmin, zmax=ns.zmax)

    # and finalize
    eboss.finalize_data(r, eboss.fidcosmo, ns.version, P0_FKP=ns.P0_FKP)

    # do the pair count
    redges = numpy.logspace(0,3.75,500) #original numpy.logspace(0,4,500)
    if ns.p == 0:
        r['Weight'] *= r['FKPWeight']
        result = SurveyDataPairCount('2d', r, redges, eboss.fidcosmo,
                                     Nmu=100, ra='RA', dec='DEC', redshift='Z',
                                     weight='Weight')
        meta = {'p': None, 'P0_FKP': ns.P0_FKP, 'z_weighted': False}
    else:
        second = r.copy()
        
        r['Weight'] *= r['FKPWeight'] * \
            eboss.bias_weight(r['Z'], eboss.ezmock_cosmo, ell=ns.ell)

        
        second['Weight'] *= second['FKPWeight'] * \
            eboss.fnl_weight(second['Z'], p=ns.p)
        result = SurveyDataPairCount('2d', r, redges, eboss.fidcosmo,
                                     second=second, Nmu=100, ra='RA', dec='DEC',
                                     redshift='Z', weight='Weight')
        meta = {'p': ns.p, 'P0_FKP': ns.P0_FKP,
                'z_weighted': True, 'ell': ns.ell}

    # save the csize
    result.attrs['N'] = r.csize

    # and save!
    kws = {'subsample': ns.subsample, 'zmin': ns.zmin,
           'zmax': ns.zmax, 'redges_str': 'logspace(0,3.75,500)',
           'ell': ns.ell, 'subsample_number': ns.nsub, 'NN weights': ns.NN_weights}
    
    kws.update(meta)
    eboss.save_RR_paircount(result, ns.sample, ns.version, ns.NN_weights, **kws)


if __name__ == '__main__':

    desc = 'compute the pair count of an eBOSS QSO randoms file'
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str,
                       choices=['N', 'S'], help=h, required=True)

    h = 'the version to load'
    group.add_argument('--version', type=str,
                       choices=eboss.DATA_VERSIONS, help=h, required=True)

    h = 'the minimum redshift to include'
    group.add_argument('--zmin', type=float, help=h, required=True)

    h = 'the maximum redshift to include'
    group.add_argument('--zmax', type=float, help=h, required=True)

    h = 'the desired collective size to subsample to'
    parser.add_argument('--subsample', type=float, help=h)

    h = 'the P0 FKP version to use'
    parser.add_argument('--P0_FKP', type=float, default=3e4, help=h)

    h = 'the value of p to use'
    group.add_argument('--p', type=float, help=h,
                       choices=[0., 1., 1.6, 3.0], required=True)

    h = 'whether to use ell=0 or ell=2 weights'
    group.add_argument('--ell', type=int, help=h,
                       choices=[0, 2], required=True)
    
    h = 'subsample number'
    parser.add_argument('--nsub', type=int, help=h, default=1)
    
    h = 'whether to use the catalog with NN systot weights'
    parser.add_argument('--NN-weights', action='store_true', help=h)

    # and go!
    main(parser.parse_args())
