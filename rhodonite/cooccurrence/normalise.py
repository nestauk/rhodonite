import numpy as np
from graph_tool import edge_endpoint_property


def association_strength(g, o_vprop, co_eprop, log=False):
    """assocation_strength
    Calculates the association strength, a symmetric and probabalistic cooccurrence
    edge weight normalisation.

    The assocation strength is calculated as defined in van Eck 2009.

    .. math::
        a = \frac{2 N c_{ij}}{o_{i} o{j}}
    
    if the graph is directed, or

    .. math::
        a = \frac{N c_{ij}}{o_{i} o{j}}
    
    if the graph is undirected, where N is the total number of cooccurrences, 
    :math:`c_{ij}` is the number of cooccurrences between vertices :math:`i` 
    and :math:`j`, and :math:`o_{i}` and :math:`o_{j}` are the respective 
    occurrence frequencies for those vertices.


    Args:
        g (:obj:`Graph`): Graph to use to calculate assocation strength.
        o_vprop (:obj:`graph_tool.PropertyMap`): A vertex property map containing
            vertex occurrence frequencies.
        co_eprop (:obj:`graph_tool.PropertyMap`): An edge property map containing
            edge cooccurrence frequencies.
        log (:obj:`bool`): If `True` association strength values are logged to
            base 10.


    Returns:
        a_s (:obj:`PropertyMap`): Assocation strength edge property map.
    """
    o_source = edge_endpoint_property(g, o_vprop, 'source')
    o_target = edge_endpoint_property(g, o_vprop, 'target')
    n_cooccurrences = np.sum(co_eprop.get_array())

    a_s = g.new_edge_property('float')
    if g.is_directed():
        a_s.a = (
            (n_cooccurrences * co_eprop.a) / 
            (o_source.a * o_target.a)
                )
    else:
        a_s.a = (
            (2 * n_cooccurrences * co_eprop.a) / 
            (o_source.a * o_target.a)
                )

    if log:
        a_s.a = np.log10(a_s.get_array())
    return a_s

def conditional_probability(g, occurrence_vprop, cooccurrence_eprop, log=False):
    """conditional_probability
    Creates an edge property with the conditional probability between source
    and target. Only works for directed graphs. The conditional probability, :math:`p`, of
    vertex :math:`i` cooccurring with vertex :math:`j`, given the occurrence frequency of vertex 
    :math:`j` is

    .. math::
        p_{i,j} = \frac{c_{ij}}{o_{j}}

    where :math:`c_{ij}` is the cooccurrence frequency between vertices :math:`i`
    and :math:`j`, and :math:`o_{j}` is the vertex occurrence frequency of 
    :math:`j`.

    Conditional probability is an assymetric normalisation. If a directed graph
    is supplied then a single edge property map is returned that contains 
    values for both edges :math:`E_{i, j}` and :math:`E_{j, i}` should they
    exist. If an undirected graph is supplied, then two edge property maps
    are returned; one containing values for :math:`p_{i,j}` and one containing
    values for every inverted edge, :math:`p_{j, i}`.

    Args:
        g (:obj:`Graph`): A graph.
        o_vprop (:obj:`graph_tool.PropertyMap`): A vertex property map containing
            vertex occurrence frequencies.
        co_eprop (:obj:`graph_tool.PropertyMap`): An edge property map containing
            edge cooccurrence frequencies.
        log (:obj:`bool`): If `True` association strength values are logged to
            base 10.

    Returns:
        c_p (:obj:`graph_tool.PropertyMap`): An edge property mapping the 
            conditional probability to each edge. Only returned if `g` is 
            directed.
        c_p_source (:obj:`graph_tool.PropertyMap`): An edge property mapping the 
            conditional probability of the source vertices given the target vertices
            for each edge. Only returned if `g` is undirected.
        c_p_target (:obj:`graph_tool.PropertyMap`): An edge property mapping the 
            conditional probability of the target vertices given the source vertices
            for each edge. Only returned if `g` is undirected.
    """
    if g.is_directed():
        c_p = g.new_edge_property('float')
        o_target = edge_endpoint_property(g, occurrence_vprop, 'target')
        c_p.a = cooccurrence_eprop.a / o_target.a
        if log:
            c_p.a = np.log10(c_p.a)
        return c_p
    else:
        c_p_source = g.new_edge_property('float')
        c_p_target = g.new_edge_property('float')
        o_target = edge_endpoint_property(g, occurrence_vprop, 'target')
        o_source = edge_endpoint_property(g, occurrence_vprop, 'source')
        c_p_source.a = cooccurrence_eprop.a / o_target.a
        c_p_target.a = cooccurrence_eprop.a / o_source.a

        if log:
            c_p_source.a = np.log10(c_p_source.a)
            c_p_target.a = np.log10(c_p_target.a)
        return c_p_source, c_p_target
