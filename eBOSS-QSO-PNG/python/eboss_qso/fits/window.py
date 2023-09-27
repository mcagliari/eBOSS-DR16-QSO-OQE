from __future__ import print_function

import argparse
import os
import numpy as np
from scipy.optimize import curve_fit
import functools

from scipy.special import legendre
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.signal import savgol_filter

_epsilon = np.finfo(float).eps

def transition(r, r_transition, b):
    """
    Tanh function used for smoothly transitioning
    between two functions
    """
    return 0.5 + 0.5*np.tanh((r-r_transition)/b)


def to_poles(filename, ells, sstart=0, sstop=3.75, sbin=500):
    """
    Open the file holding the DD(s,mu) results and convert to multipoles

    Parameters
    ----------
    filename : str
        the name of the `npz` file holding the DD(s,mu) array results
    ells : list of int
        the list of multipole numbers to compute
    sstart : float, default=0
    sstop : float, default=3.75
    sbin : int, default=500
        parameters used for the logarithic binning in the paircount (np.logspace(sstart, sstop, sbin))

    Returns
    -------
    W : array_like
        the window multipoles; first column is ``s`` and
        remaining columns are the requested multipoles
    """
    from nbodykit.algorithms import SurveyDataPairCount
    r = SurveyDataPairCount.load(filename)
    d = r.pairs

    mu = np.repeat(d.coords['mu'][None,:], d.data.shape[0], axis=0)
    weights = np.array([(2*ell+1)*legendre(ell)(mu) for ell in ells])
    N_1d = d['wnpairs'].sum(axis=-1)
    poles = (weights*d['wnpairs']).sum(axis=-1).T
    
    toret = np.zeros((len(N_1d), 1+len(ells)))
    toret[:,0] = (d['wnpairs']*d['r']).sum(axis=-1) / d['wnpairs'].sum(axis=-1)
    
    s = np.logspace(sstart, sstop, sbin)
    width = np.diff(s) / toret[:,0]
    
    toret[:,1:] = poles[:,:]/ toret[:,0][:,None]**3 / width[:,np.newaxis]
    return toret

def get_monopole_model(s, W0, smax=200., smin=5):
    """
    Fit an exponential function to the small-scale monopole,
    returning a callable function

    The model fit is

        W0(s) = A0 exp(-(r/alpha)^beta)
    """
    def model_function(r, A0, alpha, beta):
        return A0*np.exp(-(r/alpha)**beta)

    idx = (s <= smax)&(s >= smin)
    p0 = [W0[:10].mean(), 100., 1.]
    coeffs, _ = curve_fit(model_function, s[idx], W0[idx], p0, maxfev=2000)

    kws = dict(zip(['A0', 'alpha', 'beta'], coeffs))
    return functools.partial(model_function, **kws)

def get_lgt0_model(s, W, polyorder=3, smax=50, smin=20):
    """
    Fit a polynomial to the small-scale multipoles for ell > 0
    """
    # polynomial fit
    idx = (s < smax)&(s > smin)
    Wpoly = np.poly1d(np.polyfit(s[idx], W[idx], polyorder))

    # transition
    def model(s, s0=5., b=3):
        switch = transition(s, s0, b)
        return ((1-switch)*0. + switch*Wpoly(s))*(1-np.exp(-(s/1.0)**2))
    return model

def blend_models(W0, models, smin, smax, s_switch):
    """
    Extend the window s range
    """
    s0 = W0[:,0]
    if smin is not None:
        smin = min(smin, s0.min())
    else:
        smin = s0.min()
    if smax is not None:
        smax = max(smax, s0.max())
    else:
        smax = s0.max()

    s = np.logspace(np.log10(smin), np.log10(smax), 1024)
    W = np.zeros((len(s), W0.shape[1]))
    W[:,0] = s[:]

    for i in range(1, W0.shape[1]):

        spl = spline(W0[:,0], savgol_filter(W0[:,i], 31, 5))
        ell = 2*(i-1)

        switch = transition(s, s_switch[i-1], 10.)
        model_lo = models[ell](s) / models[ell](s_switch[i-1]) * spl(s_switch[i-1])
        model_hi = spl(s)
        model_hi[s <= s0.min()] = spl(s0.min())
        model_hi[s >= s0.max()] = spl(s0.max())
        W[:,i] = (1-switch)*model_lo + switch*model_hi

        # the small-scale models
        idx = s < s_switch[i-1]
        #ell = 2*(i-1)
        #W[idx,i] = models[ell](W[:,0][idx]) / models[ell](s_switch[i-1]) * spl(s_switch[i-1])


    return W

def compute_window(filename, ells, output, meta_info, smin=1e-2, smax=1e5, sstart=0, sstop=3.75, sbin=500, quiet=False):
    """
    Compute the smoothed window function correlation function multipoles
    from ``DD(s,mu)`` and save to a plaintext file

    Parameters
    ----------
    filename : str
        the input file holding the ``DD(s,mu)`` pair counts
    ells : list of int
        list of multipoles integers to compute
    output : str
        the name of the file to save the output to
    smin : float, default
        extend output to this minimum separation
    smax : float, default
        extend output to this maximum separation
    sstart : float, default=0
    sstop : float, default=3.75
    sbin : int, default=500
        parameters used for the logarithic binning in the paircount (np.logspace(sstart, sstop, sbin))
    """
    # read in the data and convert to poles
    W = to_poles(filename, ells, sstart=sstart, sstop=sstop, sbin=sbin)

    # remove NaNs
    idx = ~np.isnan(W[:,0])
    W = W[idx]

    # smax for fitting of models
    mono_smax = 100.
    lgt0_smax = 150

    # get the models
    models = {}
    for i, ell in enumerate(ells):
        if ell == 0:
            models[ell] = get_monopole_model(W[:,0], W[:,ell//2+1], smax=mono_smax, smin=10)
        else:
            models[ell] = get_lgt0_model(W[:,0], W[:,ell//2+1], smax=lgt0_smax)

    # blend the small-scale models
    s_switch = [mono_smax] + [lgt0_smax]*(len(ells)-1)
    W = blend_models(W, models, smin, smax, s_switch)

    # normalize
    #W0_norm = W[:,1][0]
    #W[:,1:] /= W0_norm

    # set small values to zero
    for i in range(1, W.shape[1]):
        idx = abs(W[:,i]) <= _epsilon
        W[idx,i] = 0.

    # save the window function
    if not quiet:
        print(f"saving {output}...")
    np.savetxt(output, W, header=str(meta_info))

def __main__():

    desc = 'compute the window multipoles from DD(s,mu)'
    parser = argparse.ArgumentParser(description=desc)

    h = 'the file to read'
    parser.add_argument('filename', type=str, help=h)

    h = 'the multipole numbers'
    parser.add_argument('--ells', required=True, type=int, nargs="+", help=h)

    h = 'the output file name'
    parser.add_argument('--output', '-o', required=True, type=str, help=h)

    h = 'extend to this minimum separation'
    parser.add_argument('--smin', type=float, default=1e-2, help=h)

    h = 'extend to this maximum separation'
    parser.add_argument('--smax', type=float, default=1e4, help=h)

    ns = parser.parse_args()
    compute_window(ns.filename, ns.ells, ns.output, smin=ns.smin, smax=ns.smax)
