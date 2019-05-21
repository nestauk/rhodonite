import itertools
import numpy as np
import pathlib

from collections import defaultdict


def get_aggregate_vp(g, vp, vp_grouper, agg=None):
    """aggregate_property_map
    
    Args:
        g (:obj:`Graph`): A graph.
        vp (str): String representing an internal property map
            of graph, g.
        vp_grouper (str): String representing name of an internal
            property map that will be used to group by.
        agg (:obj:`function`): Function to aggregate by. For
            example, min, max, sum, numpy.mean, etc.
    Returns:
        (:obj:`iter` of float): Aggregated values from x. 
    """
    vp_vals = get_vp_values(g, vp)
    vp_agg = get_vp_values(g, vp_grouper)
    
    sid_x = vp_agg.argsort()
    # Get where the sorted version of base changes groups
    split_idx = np.flatnonzero(np.diff(vp_agg[sid_x]) > 0) + 1
    # OR np.unique(base[sidx],return_index=True)[1][1:]

    # Finally sort inp based on the sorted indices and split based on split_idx
    vp_vals_grouped = np.split(vp_vals[sid_x], split_idx)
    
    x = sorted(set(vp_agg))
    if agg: 
        y = [agg(vvg) for vvg in vp_vals_grouped]
    else:
        y = vp_vals_grouped

    return x, y

def get_vp_values(g, vertex_prop_name):
    """get_vp_values
    Retrieves a vertex property from a graph, taking into account any filter.
    
    Args:
        g (:obj:`Graph`): A graph.
        vertex_prop_name (str): The name of an internal vertex property.
        
    Returns:
        pm (:obj:`PropertyMapArray`): An array of the property map.
    """
    p_type = g.vp[vertex_prop_name].value_type()
    mask = g.get_vertex_filter()[0]
    if mask is not None:
        mask = np.where(mask.get_array())
        if p_type != 'string':
            pm = g.vp[vertex_prop_name].get_array()[mask]
        else:
            pm = [g.vp[vertex_prop_name][v]
                    for m, v in zip(mask, g.vertices()) if m == True]
    else:
        pm = g.vp[vertex_prop_name].get_array()
    
    return pm

def reverse_index_communities(nested_list):
    """reverse_mapping
    Creates a mapping from elements to container index.

    Args:
        nested_list (:obj:`list`): Iterable of iterables.

    Returns:
        d (:obj:`dict`): A dictionary mapping elements in the nested iterables
            to the index of the iterable that they appeared within.
    """
    d = defaultdict(set)
    for i, l in enumerate(nested_list):
        for element in l:
            d[element].add(i)
    return d

def clear_graph(g):
    """clear_graph
    Removes all edges and vertices from the graph.
    
    Args:
        g (:obj:`Graph`): A graph.
    """

    if len(list(g.edges())) > 0:
        for e in g.edges():
            g.remove_edge(e)
    if len(list(g.vertices())) > 0:
        for v in reversed(sorted(g.vertices())):
            g.remove_vertex(v)

def edge_property_as_matrix(g, prop_name):
    """edge_property_as_matrix
    Returns an edge property as a matrix. If the graph is undirected, a
    symmetric graph is returned.
    
    Args:
        g (:obj:`Graph`): A graph-tool type graph.
        prop_name (str): The name of an edge property belonging to Graph, g.

    Returns:
     mat (:obj:`np.array`): The edge property in matrix form, with each
     position (i, j) representing the property value between vertices i and j.
    """
    n_vertices = len([_ for _ in g.vertices()])
    mat = np.zeros((n_vertices, n_vertices))
    prop = g.ep[prop_name]
    edges = sorted(g.edges())
    if g.is_directed():
        for e in edges:
            s = int(e.source())
            t = int(e.target())            
            mat[s][t] = prop[e]
    else:
        for e in edges:
            s = int(e.source())
            t = int(e.target())            
            mat[s][t] = prop[e]
            mat[t][s] = prop[e]
    return mat

def get_edge_vertex_ids(g, sort=True):
    """get_edge_vertex_ids
    Returns a list of tuples containing the source and target for each edge in
    the graph.

    Args:
        g (:obj:`Graph`): A graph with edges.
        sorted (bool): Determines whether the edges will be returned sorted in
            ascending order or not. Defaults to True.

    Returns:
        A list of source target pairs for every edge in the graph, with the
        format, (source, target).
    """
    if sort:
        return sorted([(int(e.source()), int(e.target())) for e in g.edges()])
    else:
        return [(int(e.source()), int(e.target())) for e in g.edges()]

def save_edgelist(g, filepath, weight=None):
    """save_edgelist
    Saves the edges of a graph to a file in single a space separated format to
    a delineated file.

    Some example rows from a file with no weights:

    0 1
    0 2
    1 2

    Some example rows from a file with weights:

    0 1 0.125
    0 2 0.5
    1 2 0.25
    
    Args:
        g (:obj:`Graph`): 
        filepath (str): 
        weight (str) or (:obj:`PropertyMap`): 
    """
    
    edgelist = get_edge_vertex_ids(g)

    if weight:
        if type(weight) == str:
            edge_prop = g.edge_properties[weight].get_array()
        else:
            edge_prop = weight.get_array()
        with open(filepath, 'w') as f:
            for (s, t), w in zip(edgelist, weight):
                f.write('{} {} {}\n'.format(s, t, w))
    else:
        with open(filepath, 'w') as f:
            for s, t in edgelist:
                f.write('{} {}\n'.format(s, t))

def check_and_create_dir(dir):
    """check_and_create_dir
    Checks whether an output directory exists and creates it.
    """
    pathlib.Path(dir).mkdir(parents=True, exist_ok=True)

def flatten(list_of_iters):
    """flatten
    Flattens a list of iterables into a single flat list.
    
    Args:
    list_of_iters (:obj:`iter` of :obj:`iter`): A list of iterables.
    
    Returns:
    flat (:obj:`list`): A flat list.
    """
    flat = list(itertools.chain(*list_of_iters))
    return flat

def sequence_item_types(sequence):
    """sequence_item_types"""
    if all(isinstance(item, int) for item in sequence):
        item_types = 'int'
    elif all(isinstance(item, str) for item in sequence):
        item_types = 'string'
    elif all(isinstance(item, float) for item in sequence):
        item_types = 'float'
    else:
        item_types = 'object'
    return item_types

def window(seq, n=3):
    """window
    Yields a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...

    Args:
        seq (:obj:`iter`): Sequence to move the window across.

    Examples:
        >>> doc = ['a', 'b', 'c', 'd', 'e', 'f']

        >>> list(window(doc))
        [('a', 'b', 'c'), ('b', 'c', 'd'), ('c', 'd', 'e'), ('d', 'e', 'f')]
    """
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def label_isolated(g):
    """label_isolated
    Creates a vertex property map with True if a node has no neighbours,
    else False.

    Args:
        g (graph_tool.Graph): A graph.

    Returns:
        isolated_vp (graph_tool.PropertyMap): Property map labelling isolated
            vertices.
    """"
    isolated_vp = g.new_vertex_property('bool')
    for v in g.vertices():
        if v.out_degree() == 0:
            isolated_vp[v] = True
        else:
            isolated_vp[v] = False
    return isolated_vp
