

def subgraph_eprop_values(g, vertices, eprop):
    '''get_subgraph_eprop
    Retrieves the edge property values of all edges in a subgraph.

    Args:
        g (:obj:`graph_tool.Graph`): The original graph.
        vertices (:obj:`iter`): Vertices of the subgraph.
        eprop (:obj:`graph_tool.PropertyMap`): An edge property for graph, g.

    Returns:
        vals (:obj:`array`): Edge property values from the subgraph.
    '''
    idx = g.new_vertex_property('bool')
    idx.a[vertices] = True
    gv = GraphView(g, vfilt=idx, skip_properties=True)
    vals = eprop.a[~gv.get_edge_filter()[0].ma.mask]
    return vals

def subgraph_eprop_agg(g, vertices, eprop, aggfuncs):
    '''subgraph_eprop_agg
    Retrieves the edge property values of all edges in a subgraph and returns 
    a list of user-provided aggregations.

    Args:
        g (:obj:`graph_tool.Graph`): The original graph.
        vertices (:obj:`iter`): Vertices of the subgraph.
        eprop (:obj:`graph_tool.PropertyMap`): An edge property of graph, g.
        aggfuncs (:obj:`iter`): A series of aggregation functions that can be
            applied to an array of values from the property map.
            e.g. aggfuncs=[np.mean, np.median] will return the mean and the 
            median, only retrieving the actual values from g once.

    Returns:
        agg (:obj:`list`): A list of aggregated values.
    '''
    vals = subgraph_eprop_agg(g, vertices, eprop)
    agg = [f(vals) for f in aggfuncs]
    return agg

