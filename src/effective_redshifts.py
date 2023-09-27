from nbodykit.lab import *
import eboss_qso.measurements as eboss

def select_subsample(insize, outsize):

    index = numpy.zeros(insize, dtype=bool)
    valid = numpy.random.choice(range(insize), replace=False, size=outsize)
    index[valid] = True
    return index

data_N = eboss.read_data('N', 'dr16')
data_S = eboss.read_data('S', 'dr16')

eboss.finalize_data(data_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)
eboss.finalize_data(data_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

with open(eboss.results_dir+'/effective_redshifts_dr16.txt', 'w') as f:
    print('DATA\n', file=f)
    print('NGC\n', file=f)
    print(f'''z_min = {min(data_N['Z'].compute())}, z_max = {max(data_N['Z'].compute())}\n''', file=f)
    print(f'''(z_eff, nbar_eff, norm)\n''', file=f)
    print(f'no z cut NGC data FKP: {eboss.compute_effective_quantities(data_N, eboss.fidcosmo, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=1.0: {eboss.compute_effective_quantities(data_N, eboss.fidcosmo, p=1.0, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=1.6: {eboss.compute_effective_quantities(data_N, eboss.fidcosmo, p=1.6, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=3.0: {eboss.compute_effective_quantities(data_N, eboss.fidcosmo, p=3.0, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=0.0: {eboss.compute_effective_quantities(data_N, eboss.fidcosmo, p=0.0, P0=3e4)}\n', file=f)
    
   
    
    print('SGC\n', file=f)
    print(f'''z_min = {min(data_S['Z'].compute())}, z_max = {max(data_S['Z'].compute())}\n''', file=f)
    print(f'''(z_eff, nbar_eff, norm)\n''', file=f)
    print(f'no z cut SGC data FKP: {eboss.compute_effective_quantities(data_S, eboss.fidcosmo, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=1.0: {eboss.compute_effective_quantities(data_S, eboss.fidcosmo, p=1.0, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=1.6: {eboss.compute_effective_quantities(data_S, eboss.fidcosmo, p=1.6, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=3.0: {eboss.compute_effective_quantities(data_S, eboss.fidcosmo, p=3.0, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=0.0: {eboss.compute_effective_quantities(data_S, eboss.fidcosmo, p=0.0, P0=3e4)}\n', file=f)
    
    
    
random_N = eboss.read_randoms('N', 'dr16', NN_weights=False)
random_S = eboss.read_randoms('S', 'dr16', NN_weights=False)

#numpy.random.seed(42*(random_N.comm.rank+1000+1))

#valid = select_subsample(random_N.size, int(1000000//random_N.comm.size))
#random_N = random_N[valid]
#
#valid = select_subsample(random_S.size, int(1440000//random_S.comm.size))
#random_S = random_S[valid]

eboss.finalize_data(random_N, eboss.fidcosmo, 'dr16', P0_FKP=3e4)
eboss.finalize_data(random_S, eboss.fidcosmo, 'dr16', P0_FKP=3e4)

with open(eboss.results_dir+'/effective_redshifts_dr16.txt', 'a') as f:
    print('RANDOMS\n', file=f)
    print('NGC\n', file=f)
    print(f'''z_min = {min(random_N['Z'].compute())}, z_max = {max(random_N['Z'].compute())}\n''', file=f)
    print(f'''(z_eff, nbar_eff, norm)\n''', file=f)
    print(f'no z cut NGC data FKP: {eboss.compute_effective_quantities(random_N, eboss.fidcosmo, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=1.0: {eboss.compute_effective_quantities(random_N, eboss.fidcosmo, p=1.0, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=1.6: {eboss.compute_effective_quantities(random_N, eboss.fidcosmo, p=1.6, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=3.0: {eboss.compute_effective_quantities(random_N, eboss.fidcosmo, p=3.0, P0=3e4)}\n', file=f)
    print(f'no z cut NGC data p=0.0: {eboss.compute_effective_quantities(random_N, eboss.fidcosmo, p=0.0, P0=3e4)}\n', file=f)
    
   
    
    print('SGC\n', file=f)
    print(f'''z_min = {min(random_S['Z'].compute())}, z_max = {max(random_S['Z'].compute())}\n''', file=f)
    print(f'''(z_eff, nbar_eff, norm)\n''', file=f)
    print(f'no z cut SGC data FKP: {eboss.compute_effective_quantities(random_S, eboss.fidcosmo, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=1.0: {eboss.compute_effective_quantities(random_S, eboss.fidcosmo, p=1.0, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=1.6: {eboss.compute_effective_quantities(random_S, eboss.fidcosmo, p=1.6, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=3.0: {eboss.compute_effective_quantities(random_S, eboss.fidcosmo, p=3.0, P0=3e4)}\n', file=f)
    print(f'no z cut SGC data p=0.0: {eboss.compute_effective_quantities(random_S, eboss.fidcosmo, p=0.0, P0=3e4)}\n', file=f)