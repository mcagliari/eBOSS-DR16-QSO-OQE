import numpy as np

def confidence_interval(chain, cl):
    '''Compute the confidence interval of a unimodal distribution corresponding to the confidence level cl.'''
    
    d = np.sort(np.copy(chain))
    n = len(chain)
    
    n_samples = np.floor(n * cl).astype(int)
    
    int_width = d[n_samples:] - d[:n-n_samples]
    min_int = np.argmin(int_width)
    
    return np.array([d[min_int], d[min_int+n_samples]])

def gelman_rubin_R(chains):
    '''Compute the Gelman-Rubin statistics R'''
    size = chains.shape[1]
    chain_mean = np.zeros(size)
    chain_var = np.zeros(size)
    
    for i in range(size):
        chain_mean[i] = np.mean(chains[:,i])
        chain_var[i] = np.var(chains[:,i], ddof=1)
        
    L = chains.shape[0]
    B = L * np.var(chain_mean, ddof=1)
    W = np.mean(chain_var)
    
    R = ((L - 1.) / L * W + 1. / L * B) / W
    
    return R