import numpy as np
from scipy.sparse import csr_matrix, csc_matrix

def jaccard_similarity(A, B):
    """jaccard_similarity
    Calculates the jaccard similarity between a vector a and all of the rows in
    a matrix B.

    Args:
        a (:obj:`array-like`): A vector of length L.
        B (:obj:`array-like`): A matrix or 2d array with L columns.

    Returns:
        j (:obj:`np.array`): A vector containing the row-wise jaccard
            similarities between a and B.
    """
    intersection = A.multiply(B).sum(axis=1)
    union = ((A + B) > 0).sum(axis=1)
    j = np.divide(intersection, union)
#     if isinstance(B, (csr_matrix, csc_matrix)):
#         l = B.shape[0]
#         A = np.tile(a.todense(), [l, 1])
#     else:
#         l = len(B)
#         A = np.tile(a, [l, 1])
#  
#     intersection = np.sum(np.multiply(A, B), axis=1)
#     union = np.sum(((A + B) > 0), axis=1)
#     j = np.divide(intersection, union)
    return j

