import numpy as np

from rhodonite.utilities import edge_property_as_matrix


def association_strength(g):
    """assocation_strength_adj
    Calculates the symmetric association strength matrix from the
    edge coocurrence matrix and vertex occurrence vector.

    Args:
        g (:obj:`Graph`): Graph to use to calculate assocation strength.

    Returns:
        a_s (:obj:`PropertyMap`): Assocation strength property map.
    """
    coocurrences = edge_property_as_matrix(g, 'coocurrences')
    occurrences = g.vertex_properties['occurrences'].get_array()

    a_s_mat = np.divide((2 * g.n_coocurrences * coocurences),
            np.multiply(occurrences, occurrences.transpose()))
    
    a_s = g.new_edge_property('float')
    for s, t in g.edges():
        a_s[g.edge(s, t)] = a_s_mat[s][t]

    return a_s

def occurrences(g):
    """calculate_occurrences
    Calculates the number of times each element in the sequences occurs
    and puts them in to an array ordered by the elements' vertex IDs.

    Args:
        g (:obj:`Graph`):

    Returns:
        occurrences (:obj:`PropertyMap`):
    """
    
    occurrence_counts = Counter(flatten(g.sequences_ids))
    occurrences = g.new_vertex_property('int')
    for i, count in occurrence_counts.items()
        occurrences[i] = count
    return occurrences

def coocurrences(g):
    """coocurrence
    Creates a coocurrence matrix from a counter of edge coocurrences. This
    is equivalent to an adjacency matrix, weighted by coocurrence.

    Args:
        g (:obj:`Graph`):

    Returns:
        cooccurrences (:obj:`PropertyMap`):
    """
    sequences_flat = flatten(g.sequences_ids)
    cooccurrence_counts = Counter(sequences_flat)
    cooccurrences = g.new_edge_property('int')
    for (s, t), co in coocurrence_counts.items():
        coocurrences[self.edge(s, t)] = o
    return coocurrences
