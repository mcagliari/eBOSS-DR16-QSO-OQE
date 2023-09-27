def fkp_weight(nz, P0):
    """
    Return the FKP weight.

    .. math::

        w_fkp = (1 + nz * P0)^(-1)
    """
    return 1 / (1 + nz*P0)


def bias_weight(z, cosmo, ell=0):
    """
    Compute the bias redshift weight.

    .. math::

        W = D(z) [ b(z) + f(z)/3 ] for monopole
        W = 2 D(z) f(z) / 3  for quadrupole
    """
    assert ell in [0, 2]

    b = bias_model(z)
    D = z.map_blocks(cosmo.scale_independent_growth_factor)
    f = z.map_blocks(cosmo.scale_independent_growth_rate)

    if ell == 0:
        return D*(b+f/3.)
    elif ell == 2:
        return 2./3.*D*f


def fnl_weight(z, p=1.6):
    """
    Compute the fnl redshift weight.

    .. math::

        W = (b(z) - p),

    with p = 1.0 or 1.6.
    """
    b = bias_model(z)
    weight = b - p
    #if p == 3.0:
    #    weight *= -1
    return weight


def bias_model(z):
    """
    Return the bias fit as a function of redshift from Laurent et al. 2017
    (1705.04718).
    """
    alpha = 0.278
    beta = 2.393
    return alpha * ((1+z)**2 - 6.565) + beta
