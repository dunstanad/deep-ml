def inverse_2x2(matrix: list[list[float]]) -> list[list[float]] | None:
    """
    Calculate the inverse of a 2x2 matrix.
    
    Args:
        matrix: A 2x2 matrix represented as [[a, b], [c, d]]
    
    Returns:
        The inverse matrix as a 2x2 list, or None if the matrix is singular
        (i.e., determinant equals zero)
    """
    # Your code here
    ad = 1
    bc = 1

    ad = matrix[0][0] * matrix[1][1]
    bc = matrix[0][1] * matrix[1][0]
    det = ad - bc
    if det == 0:
        return None
    
    det_inv = 1/det
    mat_inv = [ [0 for _ in range(2)] for _ in range(2)  ]
    mat_inv[0][0] = matrix[1][1]
    mat_inv[1][1] = matrix[0][0]
    mat_inv[0][1] = -1 * matrix[0][1]
    mat_inv[1][0] = -1 * matrix[1][0]

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            mat_inv[i][j] = det_inv * mat_inv[i][j]

    return mat_inv


