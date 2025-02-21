import numpy as np
import eboss_qso.measurements as eboss
from eboss_qso import EBOSS_DIR
import os

import argparse

def main(ns):
    
    in_base ='input' + ns.rebin
    
    c_file = os.path.join(EBOSS_DIR, "bis-fits", in_base, "data/bispectra", f"{ns.sample}GC", f"Bispectrum_{ns.sample}GC_EZmock_complete.dat")
    Bc = np.loadtxt(c_file)
    
    s_file = os.path.join(EBOSS_DIR, "bis-fits", in_base, "data/bispectra", f"{ns.sample}GC", f"Bispectrum_{ns.sample}GC_EZmock_shuffled.dat")
    Bs = np.loadtxt(s_file)
    
    #check ks
    assert np.all(Bc[:,0] == Bs[:,0]), "k1 grid off!"
    assert np.all(Bc[:,1] == Bs[:,1]), "k2 grid off!"
    assert np.all(Bc[:,2] == Bs[:,2]), "k3 grid off!"
    
    #Compute W_ric
    i = 4 if ns.noise_sub else 3
    W_ric = (Bc[:,i] - Bs[:,i]) / Bs[:,i]
    
    Wk_ric = np.stack([Bc[:,0],
                       Bc[:,1],
                       Bc[:,2],
                       W_ric], axis=1)
    
    SN_name = "_SNsub" if ns.noise_sub else ""
    
    output_file = os.path.join(EBOSS_DIR, "bis-fits", in_base, "data/window/Wkric/", f"{ns.sample}GC", "Wkric" + SN_name + ".dat")
    header = f"Bispectrum RIC correction for sample {ns.sample}GC and SN subtraction {ns.noise_sub} \n k1_eff k2_eff k3_eff Wk_ric"
    np.savetxt(output_file, Wk_ric, header=header)
    
    print(f"{output_file} saved!")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write in a .dat file the bispectrum radial IC correction")
    
    h = 'the sample, either North or South'
    parser.add_argument('--sample', type = str, choices = ['N', 'S'], help = h, required = True)
    
    h = 'flag to use SN subtracted B'
    parser.add_argument('--noise-sub', action = 'store_true', help=h)
    
    h = 'rebin folder'
    parser.add_argument('--rebin', type=str, default="", help=h)
    
    args = parser.parse_args()
    
    main(args)