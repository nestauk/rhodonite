from operator import itemgetter

def subgraph_eprop_values(g, vertices, eprop):
    '''get_subgraph_eprop
    Retrieves the edge property values of all edges in a subgraph.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            The original graph.
        vertices : :obj:`iter` 
            Vertices of the subgraph.
        eprop : :obj:`graph_tool.EdgePropertyMap` 
            An edge property for graph, g.

    Returns
    -------
        vals : :obj:`array` Edge property values from the subgraph.
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

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            The original graph.
        vertices : :obj:`iter` 
            Vertices of the subgraph.
        eprop : :obj:`graph_tool.PropertyMap` 
            An edge property of graph, g.
        aggfuncs : :obj:`iter` 
            A series of aggregation functions that can be applied to an array of 
            values from the property map. e.g. aggfuncs=[np.mean, np.median] 
            will return the mean and the median, only retrieving the actual 
            values from g once.

    Returns
    -------
        agg : :obj:`list` A list of aggregated values.
    '''
    vals = subgraph_eprop_agg(g, vertices, eprop)
    agg = [f(vals) for f in aggfuncs]
    return agg

def dict_to_vertex_prop(g, prop_dict, val_type):
    '''dict_to_vertex_prop
    Maps dictionary values to a vertex property using keys.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
        prop_dict : :obj:`dict` 
            A dictionary where vertices or vertex indices, and values are the 
            value of a property or weight assigned to the vertex.
        vprop_kwargs: 
            See `graph_tool.PropertyMap`.
    
    Returns
        vprop : :obj:`graph_tool.PropertyMap`
    '''
    vprop = g.new_vertex_property(val_type)
    vprop.a = itemgetter(*range(g.num_vertices()))(prop_dict)
    return vprop

def dict_to_edge_prop(g, prop_dict, val_type, edge_index_dict=None):
    '''dict_to_edge_prop
    Maps dictionary values to an edge property using keys.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph`
            A graph.
        prop_dict : :obj:`dict`
            A dictionary where keys are tuples of source target vertices or 
            indices, and values are the value of a property or weight assigned 
            to the edge.
        val_type : :obj:`str` 
            See `graph_tool.PropertyMap`.
        eprop_kwargs: 
            Other edge property map kwargs. See `graph_tool.PropertyMap`
    
    Returns
        eprop : :obj:`graph_tool.PropertyMap` 
            An edge property map.
    '''
    if edge_index_dict is None:
        edge_index_dict = {(s, t): v for s, t, v in g.get_edges()}
    eprop = g.new_edge_property(val_type)
    indices = list(itemgetter(*prop_dict.keys())(edge_index_dict))
    eprop.a[indices] = list(prop_dict.values())
    return eprop

def add_dict_edge_list(g, prop_dict, prop_name, eval_type, **eprop_kwargs):
    '''add_dict_edge_list
    Converts a dictionary to graph edges with an edge property map.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
        prop_dict : :obj:`dict` 
            A dictionary where keys are tuples of source target vertices or 
            indices, and values are the value of a property or weight assigned 
            to the edge. vertices and values are the values to be inserted into 
            the property map.
        val_type : :obj:`str`
            The type of the property. See `graph_tool.PropertyMap`.
        eprop_kwargs: 
            See `graph_tool.PropertyMap`.
    '''
    eprop = g.new_edge_property(eval_type)
    g.ep[prop_name] = eprop
    edge_list = [(e[0], e[1], v) for e, v in prop_dict.items()]
    g.add_edge_list(edge_list, eprops=[g.ep[prop_name]])

def internalize_prop_dict(g, prop_dict, prop_tag='{}'):
    '''internalize_props_from_dict
    Internalises property maps contained in a dictionary and that correspond to
    graph, g, to that graph itself.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
        prop_dict : :obj:`dict` of :obj:`graph_tool.PropertyMap` 
            A dictionary of `graph_tool.PropertyMap` for a graph, g, its vertices 
            or edges.
        prop_tag : :obj:`str`, optional
            A common name to apply to the keys of internalized properties. 
            Curly braces should be included to indicate where the unique key for 
            each property map should be included. e.g. `my_graph_prop_{}`
            Defaults to key for each property map in prop_dict.
    '''
    for k, prop in prop_dict.items():
        kt = prop.key_type()
        if kt == 'v':
            g.vp[prop_tag.format(k)] = prop
        elif kt == 'e':
            g.ep[prop_tag.format(k)] = prop
        elif kt == 'g':
            g.gp[prop_tag.format(k)] = prop

