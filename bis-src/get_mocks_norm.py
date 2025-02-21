import numpy as np
import argparse
import os

def main(ns):
    
    if ns.type == 'pk':
        lines = 18
        folder = 'Pk'
        file_start = 'Power_Spectrum_EZmock_'
    elif ns.type == 'bk':
        lines = 15
        folder = 'Bk'
        file_start = 'Bispectrum_EZmock_'
        
    folder_measurements = os.path.join('../measurements/bispectra/mocks/', ns.mock_type, f'{ns.sample}GC', folder)
    
    norms = []
    
    for i in range(1,1001):
        name = f'{ns.mock_type}_eBOSS_QSO_{ns.sample}GC_v7_{i:d}.txt' if ns.mock_type=='realistic' else f'{ns.mock_type}_eBOSS_QSO_{ns.sample}GC_{i:d}.txt'
        file_name = file_start + name
        
        with open(os.path.join(folder_measurements, file_name)) as input_file:
            head = [next(input_file) for _ in range(lines)]
            
        norm = float(''.join((ch if ch in '0123456789.+-e' else ' ') for ch in head[lines-1]))
        norms.append(norm)
        
    norms = np.array(norms)
    
    output = os.path.join('../measurements/bispectra/mocks/', f'{file_start}normalization_{ns.mock_type}_eBOSS_QSO_{ns.sample}GC_v7.txt')
    np.savetxt(output, norms, header=f'''Rustico's {ns.type} normalizations for {ns.mock_type} mocks from 1 to 1000 of {ns.sample}GC''', fmt='%1.6e')
    
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
    
    h = 'mock type'
    parser.add_argument('--mock-type', type=str, 
                       choices=['realistic', 'complete', 'shuffled'], default='realistic', help=h)
    

    # and go!
    main(parser.parse_args())
    
