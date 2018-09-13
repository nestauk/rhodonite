import pandas


def vertices_to_dataframes(g):
    vps = {}
    for k, vp in g.vp.items():
        if (vp.python_value_type != int) & (vp.python_value_type != float):
            vps[k] = [vp[v] for v in g.vertices()]
        else:
            vps[k] = vp.get_array()
    return pandas.DataFrame(vps)

def edges_to_dataframe(g):
    eps = {}
    for k, ep in g.ep.items():
        eps[k] = [g.ep[k][e] for e in g.edges()]
    sources = []
    targets = []
    for e in g.edges():
        sources.append(int(e.source()))
        targets.append(int(e.target()))
    eps['source'] = sources
    eps['target'] = targets
    return pandas.DataFrame(eps)
