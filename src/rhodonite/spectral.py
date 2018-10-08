import numpy as np

from collections import Counter

from rhodonite.utilities import edge_property_as_matrix, flatten


def association_strength(g):
    """assocation_strength_adj
    Calculates the symmetric association strength matrix from the
    edge cooccurrence matrix and vertex occurrence vector.

    The assocation strength is calculated as defined in van Eck 2009.

    a_s = 2 * total_cooccurrences * cooccurrences(a, b) /
        (occurrences(a) * occurrences(b))

    Args:
        g (:obj:`Graph`): Graph to use to calculate assocation strength.

    Returns:
        a_s (:obj:`PropertyMap`): Assocation strength property map.
    """
    cooccurrences = edge_property_as_matrix(g, 'cooccurrences')
    occurrences = np.array([g.vertex_properties['occurrences'].get_array()])

    a_s_mat = np.divide((2 * g.n_cooccurrences * cooccurrences),
            np.multiply(occurrences, occurrences.transpose()))
    a_s = g.new_edge_property('float')
    edge_vertices = [(int(e.source()), int(e.target()))
        for e in sorted(g.edges())]

    for s, t in edge_vertices:
        a_s[g.edge(s, t)] = a_s_mat[s][t]

    return a_s

def occurrences(g, sequences):
    """calculate_occurrences
    Calculates the number of times each element in the sequences occurs
    and puts them in to an array ordered by the elements' vertex IDs.

    Args:
        g (:obj:`Graph`):
        sequences (:obj:`iter` of :obj:`iter`):

    Returns
        o (:obj:`PropertyMap`): 
    """
    counts = Counter(flatten(sequences))
    o = g.new_vertex_property('int')
    for k, v in counts.items():
        o[k] = v
    return o

def cooccurrences(g):
    """cooccurrence
    Creates a cooccurrence matrix from a counter of edge cooccurrences. This
    is equivalent to an adjacency matrix, weighted by cooccurrence.

    Args:
        g (:obj:`Graph`):

    Returns:
        cooccurrences (:obj:`PropertyMap`):
    """
    co_sequences_flat = flatten(g.id_co_sequences)
    cooccurrence_counts = Counter(co_sequences_flat)
    cooccurrences = g.new_edge_property('int')
    for (s, t), co in cooccurrence_counts.items():
        s = g.tokenidx2vertex[s]
        t = g.tokenidx2vertex[t]
        cooccurrences[g.edge(s, t)] = co
    return cooccurrences
