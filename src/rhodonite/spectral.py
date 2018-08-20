import numpy as np


def association_strength(g):
    """assocation_strength_adj
    Calculates the symmetric association strength matrix from the
    edge coocurrence matrix and vertex occurrence vector.

    Args:
        n_cooccurrences (int): Total number of coocurrences.
        coocurs (:obj:`np.matrix`): Coocurrence matrix of the network.
        occurs (:obj:`np.array`): Array containing the total number of 
            occurrences for each vertex in the network.

    Returns:
        a_s (:obj:`np.matrix`): The association strength matrix for the
            coocurrence matrix.
    """
    

    a_s = np.divide((2 * g.n_coocurrences * g.coocurs),
            np.multiply(g.occurs, g.occurs.transpose()))
    return a_s
