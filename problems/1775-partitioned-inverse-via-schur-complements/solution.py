import numpy as np

def partitioned_inverse(M, k):
    # Split into blocks 
    E = M[:k, :k]
    F = M[:k, k:]
    G = M[k:, :k]
    H = M[k:, k:]

    # H inverse
    H_inv = np.linalg.inv(H)

    # Schur Complement wrt H 
    S = E - F@H_inv@G

    # S Inverse
    S_inv = np.linalg.inv(S)

    #Individual blocks of M Inverse
    blockA = S_inv
    blockB = -S_inv@F@H_inv
    blockC = -(H_inv@G@S_inv)
    blockD = H_inv + H_inv@G@S_inv@F@H_inv

    M_inv = np.block([
        [blockA, blockB],
        [blockC, blockD]
    ])
    return M_inv

    

def block_slogdet(M, k):
    # Split into blocks 
    E = M[:k, :k]
    F = M[:k, k:]
    G = M[k:, :k]
    H = M[k:, k:]

    # H inverse
    H_inv = np.linalg.inv(H)

    # Schur Complement wrt H 
    S = E - F@H_inv@G

    sign_S, log_S = np.linalg.slogdet(S)
    sign_H, log_H = np.linalg.slogdet(H)

    return (sign_H * sign_S, log_H + log_S)

def bordered_inverse(A_inv, b, c):
    A_inv = np.asarray(A_inv, dtype= float)
    b = np.asarray(b, dtype = float).reshape(-1, 1)
    A_inv_b = A_inv @ b

    #Schur Complement wrt A
    S = float((c - (b.T@ A_inv_b)).item()) 

    if abs(S) < 1e-12:
        raise ValueError("S is too close to 0")

    #S Inverse 
    S_inv = 1.0/S
    
    bT_A_inv = b.T@A_inv

    blockA = A_inv + A_inv_b*S_inv@bT_A_inv
    blockB = - A_inv_b*S_inv
    blockC = - S_inv*bT_A_inv
    blockD = np.array([[S_inv]])

    M_inv = np.block([
        [blockA, blockB],
        [blockC, blockD]
    ])

    return M_inv
    
    
