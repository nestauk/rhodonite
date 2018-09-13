import pandas


def vertices_to_dataframes(g, keys=None):
    """vertices_to_dataframe
    Transforms a graph's vertices and their properties into a tabular format.

    Args:
        g (:obj:`Graph`): A graph.
        keys (:obj:`iter` of str): A list of property map keys to convert in
            to columns. If None, all property maps are converted. Default is
            None.
    Returns:
        df (:obj:`DataFrame`): A dataframe where each row represents a vertex
            and the columns are properties of those edges. By default, the 
            dataframe will contain a column with the vertex id.
    """
    vps = {}

    vps['vertex'] = [v for v in g.vertices()]

    if keys is not None:
        for k in keys:
            vp = g.vp[k]
            if (vp.python_value_type != int) & (vp.python_value_type != float):
                vps[k] = [vp[v] for v in g.vertices()]
            else:
                vps[k] = vp.get_array()
    else:
        for k, vp in g.vp.items():
            if (vp.python_value_type != int) & (vp.python_value_type != float):
                vps[k] = [vp[v] for v in g.vertices()]
            else:
                vps[k] = vp.get_array()
    return pandas.DataFrame(vps)

def edges_to_dataframe(g, keys=None):
    """edges_to_dataframe
    Transforms a graph's edges and their properties into a tabular format.

    Args:
        g (:obj:`Graph`): A graph.
        keys (:obj:`iter` of str): A list of property map keys to convert in
            to columns. If None, all property maps are converted. Default is
            None.
    Returns:
        df (:obj:`DataFrame`): A dataframe where each row represents an edge
            and the columns are properties of those edges. By default, the 
            dataframe will contain a column for the source vertices and another
            for the target vertices.
    """
    eps = {}

    sources = []
    targets = []
    for e in g.edges():
        sources.append(int(e.source()))
        targets.append(int(e.target()))
    eps['source'] = sources
    eps['target'] = targets

    if keys is not None:
        for k in keys:
            eps[k] = [g.ep[k][e] for e in g.edges()]
    else:
        for k, ep in g.ep.items():
            eps[k] = [ep[e] for e in g.edges()]
    df = pandas.DataFrame(eps)
    return df
