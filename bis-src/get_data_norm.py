import numpy as np
import argparse
import os

def main(ns):
    ### CHECK NAMES
    if ns.type == 'pk':
        lines = 18
        folder = 'Pk'
        file_start = 'Power_Spectrum_'
    elif ns.type == 'bk':
        lines = 15
        folder = 'Bk'
        file_start = 'Bispectrum_'
        
    folder_measurements = os.path.join('../measurements/bispectra/',  'data', folder)
    
    file_name = file_start + f'eBOSS-QSO-{ns.sample}GC_dk0p005_kmin0p0025_kmax0p08_PCS_step2_Ngrid512.txt'
    
    if ns.NN_weights:
        file_name = file_start + f'eBOSS-QSO-NNweights-{ns.sample}GC_dk0p005_kmin0p0025_kmax0p08_PCS_step2_Ngrid512.txt'
        
    with open(os.path.join(folder_measurements, file_name)) as input_file:
        head = [next(input_file) for _ in range(lines)]
        
    norm = np.array([float(''.join((ch if ch in '0123456789.+-e' else ' ') for ch in head[lines-1]))])
    
    output = os.path.join('../measurements/bispectra/data/', f'{file_start}normalization_eBOSS_QSO_{ns.sample}GC.txt')
    if ns.NN_weights:
        output = os.path.join('../measurements/bispectra/data/', f'{file_start}normalization_eBOSS_QSO_{ns.sample}GC-NN_weights.txt')
    np.savetxt(output, norm, header=f'''Rustico's {ns.type} normalizations for {ns.sample}GC''', fmt='%1.6e')
    
    print(output, 'ready!')
    
    
if __name__=="__main__":
    desc = '''Produces a .txt file containing the Rustico's normalizations'''
    parser = argparse.ArgumentParser(description=desc)

    # required arguments
    group = parser.add_argument_group('required arguments')

    h = 'the sample, either North or South'
    group.add_argument('--sample', type=str,
                       choices=['N', 'S'], help=h, required=True)

    h = 'output type'
    group.add_argument('--type', type=str, 
                       choices=['pk', 'bk'], help=h, required=True)
    
    h = 'use NN weight catalog flag'
    parser.add_argument("--NN-weights", action='store_true', help=h)
    
    # and go!
    main(parser.parse_args())
    