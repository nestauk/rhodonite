import numpy as np

from collections import Counter
from graph_tool import edge_endpoint_property

from rhodonite.utilities import edge_property_as_matrix, flatten


def association_strength(g, occurrence_vprop, cooccurrence_eprop, norm=False, log=False):
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

    a_s = g.new_edge_property('float')

    o_source = edge_endpoint_property(g, occurrence_vprop, 'source')
    o_target = edge_endpoint_property(g, occurrence_vprop, 'target')
    n_cooccurrences = np.sum(cooccurrence_eprop.get_array())

    a_s.a = (
        (2 * n_cooccurrences * cooccurrence_eprop.a) / 
        (o_source.a * o_target.a)
            )

    if log:
        a_s.a = np.log10(a_s.get_array())
    return a_s

def conditional_probability(g, occurrence_vprop, cooccurrence_eprop, log=False):
    """conditional_probability
    Creates an edge property with the conditional probability between source
    and target. Only works for directed graphs. The conditional probability is
    calculated as the cooccurrence weight of the edge divided by the occurrence
    weight of the target.

    Args:
        g (:obj:`Graph`):

    Returns:
        conditional_probability_prop (:obj:`PropertyMap`): An edge property 
            mapping the conditional probability to each edge.
    """
    if g.is_directed():
        conditional_probability_prop = g.new_edge_property('float')
        o_target = edge_endpoint_property(g, occurrence_vprop, 'target')
        conditional_probability_prop.a = cooccurrence_eprop.a / o_target.a
        if log:
            conditional_probability_prop.a = np.log10(
                    conditional_probability_prop.a)
        return conditional_probability_prop

