import numpy as np

def woodbury_solve(sigma_diag, X, B):
    sig_inv = 1/sigma_diag
    if B.ndim == 2:
        sig_invB = sig_inv[:,np.newaxis] * B
    else:
        sig_invB = sig_inv * B
    sig_invX = sig_inv[:, np.newaxis] * X
    I = np.eye(X.shape[1])

    A = I + X.T @ sig_invX
    C = X.T @ sig_invB
    rhs = np.linalg.solve(A,C)

    return sig_invB - sig_invX@rhs
    

def sherman_morrison(A_inv, u, v):
    Au = A_inv@u
    den = 1 + v @ Au
    if den == 0:
        raise ValueError("Update is singular")
    
    num = np.outer(Au, v.T @ A_inv)
    

    return A_inv - (num/den)
    
    

def logdet_rank_one_update(logdet_A, A_inv, u, v):
    factor = 1 + v.T@A_inv@u
    if factor <=0:
        raise ValueError("Updated determinant is not positive")
    
    return float(logdet_A + np.log(factor) )