import numpy as np
from classy import Class
from scipy import interpolate
from eboss_qso import EBOSS_FITS
import os

import argparse

k_file = os.path.join(EBOSS_FITS, 'input', 'models', 'table_p.dat')

k = np.loadtxt(k_file)

class_params = {
            "output": "dTk mPk", 
            #"extra metric transfer functions": 'yes', 
            "n_s": 0.9652, 
            "gauge": "synchronous", 
            "N_ur": 2.0328, 
            "h": 0.6766, 
            "ln10^{10}A_s": 3.043, 
            "k_pivot": 0.05, 
            "tau_reio": 0.054, 
            "T_cmb": 2.7255, 
            "Omega_cdm": 0.26069, 
            "Omega_b": 0.04897, 
            "Omega_k": 0.0, 
            "N_ncdm": 1, 
            "m_ncdm": 0.06, 
            "P_k_max_h/Mpc": 10.0, 
            "z_max_pk": 100.0
        }


c = Class()
c.set(class_params)
c.compute()

def pk(cosmo, ks, z_pk):
    h = cosmo.h()
    pk = np.array([np.array([cosmo.pk(p, z_pk) for p in k * h]) for k in ks]) * h**3.
    return pk

def tk_raw(cosmo, z_pk, save=False):
    tk = cosmo.get_transfer(z_pk)
    if save:
        tkarray = np.zeros(np.array(tk['d_tot']).shape)
        for n in tk.keys():
            print(n, tkarray.shape)
            if tkarray.shape == np.array(tk['d_tot']).shape:
                tkarray=np.stack((tkarray,np.array(tk[n]))).transpose()
            else:
                tkarray=np.concatenate((tkarray,np.array(tk[n])[:,np.newaxis]), axis=1)
        return tkarray[:,1:], tk.keys()
    return tk

def tk(cosmo, ks, z_pk):
    tk = tk_raw(cosmo=cosmo, z_pk=z_pk, save=False)
    Omega_cdm = cosmo.Omega0_cdm()
    Omega_b = cosmo.Omega_b()
    tkm = Omega_cdm / (Omega_cdm + Omega_b) * tk['d_cdm'] + Omega_b / (Omega_cdm + Omega_b) * tk['d_b']
    Tk = -tkm/tk['k (h/Mpc)']**2
    Tk /= Tk[0]
    
    f = interpolate.interp1d(tk['k (h/Mpc)'], Tk)
    interpTk = f(ks)
    
    return interpTk

def Dz_n_fz(z_pk, Df_file=EBOSS_FITS + "/input/models/Dz_fz_Planck18.txt"):
    z, D, f = np.loadtxt(Df_file, unpack=True)
    
    Doz = interpolate.interp1d(z, D)
    foz = interpolate.interp1d(z, f)
    
    interpD = Doz(z_pk)
    interpf = foz(z_pk)
    return interpD, interpf

def alphatilde(cosmo, ks, Tk, Dz):
    H0 = 100. #the term h**2 in in Tk and Dz
    Om0 = cosmo.Omega0_m()
    deltac = 1.686
    c_light  = 299792.458 #speed of light [km/s]
    alphat = 3. * Om0 * H0**2 * deltac / c_light**2 / ks**2 / Tk / Dz
    return alphat

def main(ns):
    folder_models = os.path.join(EBOSS_FITS, 'input', 'models')
    sample = 'NGC' if ns.sample == 'N' else 'SGC'
    folder_sample = os.path.join(folder_models, sample)
    folder_p = os.path.join(folder_sample, f'{ns.p:2.1f}')
    
    #P(k)
    pok = pk(c, k, ns.zeff)
    pk_file = os.path.join(folder_p, 'Pk_model.dat')
    np.savetxt(pk_file, np.hstack([k[:,0][:,np.newaxis] * 10,pok]), header=f"P(k,p) at z_eff={ns.zeff} (p={ns.p}, {sample}) \n k(h/Mpc) Pk(p)")
    
    #T(k)
    if ns.raw:
        tokraw, keys = tk_raw(c, ns.zeff, save=True)
        tk_raw_file = os.path.join(folder_p, 'Tk_model_raw.txt')
        np.savetxt(tk_raw_file, tokraw, header=f"Tk at z_eff={ns.zeff} (p={ns.p}, {sample}) raw output of classy \n {keys}")
        
    tok = tk(c, k, ns.zeff)
    tk_file = os.path.join(folder_p, 'Tk_model.dat')
    np.savetxt(tk_file, np.hstack([k[:,0][:,np.newaxis] * 10,tok]), header=f"T(k,p) at z_eff={ns.zeff} (p={ns.p}, {sample}) \n k(h/Mpc) Tk(p)")
    
    #Dz fz
    doz, foz = Dz_n_fz(ns.zeff, ns.Dffile)
    dfz_file = os.path.join(folder_p, 'Dz_fz_model.dat')
    np.savetxt(dfz_file, np.array([ns.zeff, doz, foz]), header=f"z_eff D f")
    
    #alphak
    alpha = alphatilde(c, k, tok, doz)
    alpha_file = os.path.join(folder_p, 'alphatildek_model.dat')
    np.savetxt(alpha_file, np.hstack([k[:,0][:,np.newaxis] * 10, alpha]), header=f"alphatilde(k,p) at z_eff={ns.zeff} (p={ns.p}, {sample}) \n k(h/Mpc) alphatildek(p)")
    
    print(f'Model in {folder_p} is ready!')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compute the model P(k), T(k), D, f and alphatilde at a given z_eff")
    
    h = 'effective redshift'
    parser.add_argument('--zeff', type = float, help = h, required = True)
    
    h = 'the sample, either North or South'
    parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)
    
    h = 'the value of p to use'
    parser.add_argument('--p', type = float, help = h, choices=[0., 1., 1.6, 3.0], required = True)
    
    h = 'path to D(z) and f(z) file'
    Df_file = os.path.join(EBOSS_FITS, 'input', 'models', 'Dz_fz_Planck18.txt')
    parser.add_argument('--Dffile', type = str, help = h, default = Df_file)
    
    h = 'flag for raw Tk file'
    parser.add_argument('--raw', action = 'store_true')
    
    args = parser.parse_args()
    
    main(args)