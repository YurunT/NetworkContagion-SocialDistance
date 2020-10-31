from scipy.special import gammaln
import numpy as np
def log_factorial(x):
    """Returns the logarithm of x!
    Also accepts lists and NumPy arrays in place of x."""
    return gammaln(np.array(x)+1)

def get_log_multinomial_coeffecient(N_vec):
    '''
    N_vec = (x1, x2, ..., xk)
    sum(N_vec) = n
    return n!/(x1!*x2!*...*xk!)
    '''
    n_range = sum(N_vec) # should be n_range - 1
    result = log_factorial(n_range) - sum(log_factorial(N_vec)) 
    return result

def div(nu, de):
    if de == 0:
        res = 0
    else:
        res = nu * 1.0 / de
    return res