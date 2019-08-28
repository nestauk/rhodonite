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
    return j

def jaccard_similarity_set(a, b):
    """jaccard_similarity"""
    a = set(a)
    b = set(b)
    intersection = len(a.intersection(b))
    union = len(a.union(b))
#     intersection = len(list(set(a).intersection(b)))
#     union = (len(a) + len(b)) - intersection
    return intersection / union

