import pandas


def vertices_to_dataframe(g, keys=None):
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
    for k, vp in g.vp.items():
        vt = vp.value_type()
        if ('vector' not in vt) & ('string' not in vt) & ('object' not in vt):
            if ('int' in vt) | ('bool' in vt):
                vertex_df[k] = vp.get_array()
                vertex_df[k] = vertex_df[k].astype(int)
            elif 'double' in vt:
                vertex_df[k] = vp.get_array()
                vertex_df[k] = vertex_df[k].astype(float)
    return vertex_df

def edges_to_dataframe(g, keys=None):
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
    for k, ep in g.ep.items():    
        vt = ep.value_type()
        if ('vector' not in vt) & ('string' not in vt) & ('object' not in vt):
            if ('int' in vt) | ('bool' in vt):
                edge_df[k] = ep.get_array()
                edge_df[k] = edge_df[k].astype(int)
            elif 'double' in vt:
                edge_df[k] = ep.get_array()
                edge_df[k] = edge_df[k].astype(float)
    return edge_df
