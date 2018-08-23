import numpy as np


def jaccard_similarity(a, B):
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
    A = np.tile(a, [len(B), 1])
    intersection = np.sum(np.multiply(A, B), axis=1)
    union = np.sum(((A + B) > 0).astype(int))
    j = np.divide(intersection, union)
    return j

