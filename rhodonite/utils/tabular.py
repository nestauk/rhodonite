import pandas as pd
import numpy as np


def vertices_to_dataframe(g, vectors=False, keys=None):
    """vertices_to_dataframe
    Transforms a graph's vertices and their properties into a tabular format.

    Args:
        g (:obj:`Graph`): A graph.
        keys (:obj:`iter` of str): A list of property map keys to convert in
            to columns. If None, all property maps are converted. Default is
            None.
    Returns:
        vertex_df (:obj:`DataFrame`): A dataframe where each row represents a 
            vertex and the columns are properties of those edges. By default, 
            the dataframe will contain a column with the vertex id.
    """
    vertex_df = pd.DataFrame(list(g.vertices()), columns=['v'], dtype='int')
    filt = g.get_vertex_filter()
    if filt[0] is not None:
        filt = filt[0].a > 0
        filtered = True
    else:
        filtered=False
        filt = np.array([True] * g.num_vertices())
    if keys is None:
        keys = list(g.vp.keys())
    for k, vp in g.vp.items():
        if k in keys:
            vt = vp.value_type()
            if ('vector' not in vt) & ('string' not in vt) & ('object' not in vt):
                if ('int' in vt) | ('bool' in vt):
                    vertex_df[k] = vp.get_array()[filt]
                    vertex_df[k] = vertex_df[k].astype(int)
                elif 'double' in vt:
                    vertex_df[k] = vp.get_array()[filt]
                    vertex_df[k] = vertex_df[k].astype(float)
            elif ('vector' in vt) & ('string' not in vt) & (vectors == True):
                vertex_df[k] = [[i for i in vp[v]] for v, f in zip(g.vertices(), filt) if f]
    return vertex_df

def edges_to_dataframe(g, keys=None, drop_keys=[], sort=None, vectors=False):
    """edges_to_dataframe
    Transforms a graph's edges and their properties into a tabular format.

    Args:
        g (:obj:`Graph`): A graph.
        keys (:obj:`iter` of str): A list of property map keys to convert in
            to columns. If None, all property maps are converted. Default is
            None.
    Returns:
        edge_df (:obj:`DataFrame`): A dataframe where each row represents an 
            edge and the columns are properties of those edges. By default, the 
            dataframe will contain a column for the source vertices and another
            for the target vertices.
    """
    edge_df = pd.DataFrame(list(g.edges()), columns=['s', 't'], dtype='int')
    edge_df = pd.DataFrame(g.get_edges(), columns=['s', 't', 'e_index'], dtype='int')
    indices = edge_df['e_index'].values
    if keys is None:
        keys = list(g.ep.keys())
    for k, ep in g.ep.items():
        if k in keys:
            vt = ep.value_type()
            if ('vector' not in vt) & ('string' not in vt) & ('object' not in vt):
                if ('int' in vt) | ('bool' in vt):
                    edge_df[k] = ep.get_array()[indices]
                    edge_df[k] = edge_df[k].astype(int)
                elif 'double' in vt:
                    edge_df[k] = ep.get_array()[indices]
                    edge_df[k] = edge_df[k].astype(float)
            elif ('vector' in vt) & ('string' not in vt) & (vectors == True):
                edge_df[k] = [[i for i in vp[v]] for v in g.vertices()]
    if sort:
        edge_df.sort_values(sort, inplace=True)
    return edge_df
