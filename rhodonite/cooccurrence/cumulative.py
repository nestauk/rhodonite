from collections import Counter
from graph_tool import Graph

from rhodonite.utilities.misc import dict_to_vertex_prop
from rhodonite.cooccurrence.basic import cooccurrence_counts


def cumulative_cooccurrence_graph(steps, sequences, directed=False):
    '''cumulative_cooccurrence_graph

    Args:
        steps (:obj:`iter` of :obj:`int` or :obj:`str`): A series that contains
            sequential labels for the nested groups.
        sequences (:obj:`iter` of :obj:`iter` of :obj:`int`): Nested iterable
            of integers representing vertices in the graph. Number of nested
            iterables should be equal to `len(steps)`.
        directed (:obj:`bool`): Currently has no effect. In future this will
            determine whether to build a bi-directional cooccurrence graph.
    Returns:
        g (:obj:`graph_tool.Graph`): A graph. Vertices are elements. Edges
            link terms that have cooccurred at least once in the series.
        o_props (:obj:`dict`): Property maps with vertex occurrence values at
             each step.
        o_cumsum_props (:obj:`dict`): Property maps with cumulative vertex
            cooccurrence values at each step.
        co_props (:obj:`dict`): Property maps with edge cooccurrnce values
            at each step.
        co_cumsum_props (:obj:`dict`): Property maps with cumulative edge
            cooccurrence values at each step.
    '''

    g = Graph(directed=directed)

    o_total = Counter(chain(*sequences))
    n_vertices = len(o_total)
    g.add_vertex(n_vertices)
    o_max = dict_to_vertex_prop(g, o_total)

    co_total = cooccurrence_counts(chain(*sequences))
    edge_list = ((c[0], c[1], count) for c, count in co_total.items())
    co_max = g.new_edge_property('int')
    g.add_edge_list(edge_list, eprops=[co_max])

    edges = g.get_edges()
    edge_indices = dict(zip([(e[0], e[1]) for e in edges], edges[:, 2]))

    o_cumsum = g.new_vertex_property('int')
    co_cumsum = g.new_edge_property('int')
    o_props = {}
    co_props = {}
    o_cumsum_props = {}
    co_cumsum_props = {}
    for step in steps[:-1]:
        o_step = Counter(chain(*sequence_dict[step]))
        o_props[step] = dict_to_vertex_prop(g, o_step)
        o_cumsum.a = o_cumsum.a + o_props[step].a
        o_cumsum_props[step] = o_cumsum
        
        combos = (combinations(sorted(ids), 2) for ids in sequence_dict[step])
        co_step = Counter(chain(*combos))
        co_props[step] = dict_to_edge_prop(g, co_step, edge_indices)
        co_cumsum.a = co_cumsum.a + co_props[step].a
        co_cumsum_props[step] = co_cumsum
    
    # fill in the last step without needing to count occurrences
    # or cooccurrences
    step_max = steps[-1]
    o = g.new_vertex_property('int')
    co = g.new_edge_property('int')
    o.a = o_max.a - o_cumsum.a
    co.a = co_max.a - co_cumsum.a
    o_props[step_max] = o
    co_props[step_max] = co

    o_cumsum_props[step_max] = o_max
    co_cumsums_props[step_max] = co_max
    
    return g, o_props, o_cumsum_props, co_props, co_cumsum_props

def label_cumulative_cooccurrence(g, co_props):
    '''label_cumulative_cooccurrence
    Calculates the cumulative cooccurrence values for edges in a graph using
    the edge cooccurrence values.

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        co_props (:obj:`dict`): Contains `PropertyMaps` of edge cooccurrence
            values.
    Returns:
        co_cumsum_props (:obj:`dict`): Contains `PropertyMaps` of type `int` with
        edge cumulative cooccurrence values.
    '''
    co_cumsum_props = {}
    for i, (step, co_prop) in enumerate(co_props.items()):
        co_cumsum_prop = g.new_edge_property('int')
        if i == 0:
            co_cumsum_prop.a = co_prop.a
        else:
            co_cumsum_prop.a = co_cumsum_props[step-1].a + co_prop.a
        co_cumsum_props[step] = co_cumsum_prop
    return co_cumsum_props

def label_cumulative_occurrence(g, o_props):
    """label_cumulative_occurrence
    Calculates the cumulative occurrence values for vertices in a graph using
    the vertex occurrence values.

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        o_props (:obj:`dict`): Contains `PropertyMaps` of vertex occurrence
            values.
    Returns:
        o_cumsum_props (:obj:`dict`): Contains `PropertyMaps` of type 'int' with
        vertex cumulative occurrence values.
    """
    o_cumsum_props = {}
    for i, (step, o_prop) in enumerate(o_props.items()):
        o_cumsum_prop = g.new_vertex_property('int')
        if i == 0:
            o_cumsum_prop.a = o_prop.a
        else:
            o_cumsum_prop.a = (o_cumsum_props[step-1].a + o_prop.a)
        o_cumsum_props[step] = o_cumsum_prop
    return o_cumsum_props

def label_active_edge(g, co_props):
    '''label_active_edge
    Determines whether an edge is active (has at least one cooccurrence) at each
    step.

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        co_props (:obj:`dict`): Contains `PropertyMaps` of edge cooccurrence
            values.
    Returns:
        active_edge_props (:obj:`dict`): Contains `PropertyMaps` of type `bool`
            that is True if an edge is active.
    '''
    active_edge_props = {}
    for step, co_prop in co_props.items():
        active_edge_prop = g.new_edge_property('bool')
        active_edge_prop.a = co_prop.a > 0
        active_edge_props[step] = active_prop
    return active_edge_props

def label_new_edge(g, co_props, label_old=True, label_steps=False):
    '''label_new_edge
    Determines whether an edge has appeared for the first time at a given step.

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        co_props (:obj:`dict`): Contains `PropertyMaps` of edge cooccurrence
            values.
        label_old (:obj:`bool`): If True will also return `PropertyMaps` that
            indicate whether an edge has existed at a previous step.
        label_steps (:obj:`bool`): If True will also return a `PropertyMap` of
            type `int` that indicates the step at which an edge first appeared.

    Returns:
        new_edge_props (:obj:`dict`):
        old_edge_props (:obj:`dict`):
        edge_step_prop (:obj:`PropertyMap`):
        
        
    '''
    new_edge_props = {}
    if label_old:
        old_edge_props = {}

    _edge_tracker = g.new_edge_property('bool')

    for step, co_prop in co_props.items():
        new_edge_prop = g.new_edge_property('bool')
        new_edge_prop.a = (co_prop.a > 0) & (_edge_tracker.a == False)
        new_edge_props[step] = new_edge_prop
        if label_old:
            old_edge_prop = g.new_edge_property('bool')
            old_edge_prop.a = _edge_tracker.a
            old_edge_props[step] = old_edge_prop
        _edge_tracker.a = _edge_tracker.a + new_edge_prop.a
        
    if label_steps:
        steps = list(co_props.keys())
        edge_step_prop = g.new_edge_property('int')
        start = 0
        end = np.sum(new_edge_props[steps[0]].a)
        for i, step in enumerate(steps):
            if i > 0:
                start = int(start + np.sum(new_edge_props[steps[i-1]].a))
                end = int(end + np.sum(new_edge_props[step].a))
            edge_step_prop.a[start:end] = step

    if label_old & (not label_steps):
        return (new_edge_props, old_edge_props)
    elif label_steps & (not label_old):
        return (new_edge_props, edge_step_prop)
    elif label_old & label_steps:
        return (new_edge_props, old_edge_props, edge_step_prop)
    else:
        return new_edge_props

def label_edge_activity(g, co_props):
    '''label_edge_activity
    Determines whether an edge is new, reinforcing or inactive at any given 
    step.

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        co_props (:obj:`dict`): Contains `PropertyMaps` of edge cooccurrence
            values.

    Returns:
        new_edge_props (:obj:`dict`):
        reinforcing_edge_props (:obj:`dict`):
        inactive_edge_props (:obj:`dict`):
    '''
    reinforcing_edge_props = {}
    inactive_edge_props = {}

    active_edge_props = label_active_edge(g, co_props)
    new_edge_props, old_edge_props = label_new_edge(g, co_props)
    for step in active_edge_props.keys():
        reinforcing_eprop = g.new_vertex_property('bool')
        reinforcing_eprop.a = active_edge_props[step].a & old_edge_props[step].a
        reinforcing_edge_props[step] = reinforcing_eprop

        inactive_eprop = g.new_vertex_property('bool')
        inactive_eprop.a = ((old_edge_props[step].a > 0) 
                            & (active_edge_props[step].a == 0))
        inactive_edge_props[step] = inactive_eprop

    return new_edge_props, reinforcing_edge_props, inactive_edge_props

def label_new_vertex(g, o_props, label_steps=False):
    '''label_new_vertex

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        o_props (:obj:`dict`): A dictionary of `PropertyMaps` containing the 
            vertex occurrence values at each step.
        label_steps (:obj:`bool`): Returns a `PropertyMap` that indicates the
            first step that each vertex appeared.

    Returns:
        new_vertex_props (:obj:`dict`):
        vertex_step_prop (:obj:`graph_tool.PropertyMap`):
    '''
    new_vertex_props = {}
    steps = sorted(steps)
    _vertex_tracker = g.new_vertex_property('bool')
    for step, o_prop in o_props.items():
        new_vertex_prop = g.new_vertex_property('bool')
        new_vertex_prop.a = (o_prop.a > 0) & (_vertex_tracker.a == False)
        new_vertex_props[step] = new_vertex_prop
        _vertex_tracker.a = _vertex_tracker.a + new_vertex_prop.a
    if label_steps:
        vertex_step_prop = g.new_vertex_property('int')
        start = 0
        end = np.sum(new_vertex_props[steps[0]].a)
        for i, step in enumerate(steps):
            if i > 0:
                start = int(start + np.sum(new_vertex_props[step - 1].a))
                end = int(end + np.sum(new_vertex_props[step].a))
            vertex_step_prop.a[start:end] = step
        return (new_vertex_props, vertex_step_prop)
    else:
        return new_vertex_props

def label_linkage(g, new_eprops, new_vprops):
    """label_linkage
    Labels new edges as fresh, mixed or mature.
     - Fresh: edges that connect two vertices that have both appeared for the
        first time at a given step.
     - Mixed: edges that connect a new vertex with a vertex that was already
        in the graph.
     - Mature: edges that connect two pre-existing vertices.

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        new_eprops (:obj:`dict`):
        new_vprops (:obj:`dict`):

    Returns:
        fresh_eprops (:obj:`dict`):
        mixed_eprops (:obj:`dict`):
        mature_eprops (:obj:`dict`):
    """
    
    fresh_eprops = {}
    mixed_eprops = {}
    mature_eprops = {}
    
    for step, new_eprop in new_eprops.items():
        new_vprop = new_vprops[step]
        s_new_prop = g.new_edge_property('bool')
        t_new_prop = g.new_edge_property('bool')
        s_new_prop = edge_endpoint_property(g, new_vprop, 'source')
        t_new_prop = edge_endpoint_property(g, new_vprop, 'target')
    
        fresh_eprop = g.new_edge_property('bool')
        mixed_eprop = g.new_edge_property('bool')
        mature_eprop = g.new_edge_property('bool')
    
        fresh_eprop.a = (s_new_prop.a & t_new_prop.a)  & (new_eprop.a)
        mixed_eprop.a = (s_new_prop.a != t_new_prop.a)  & (new_eprop.a)
        mature_eprop.a = ((fresh_eprop.a == False) 
                          & (mixed_eprop.a == False) 
                          & (new_eprop.a))
        
        fresh_eprops[step] = fresh_eprop
        mixed_eprops[step] = mixed_eprop
        mature_eprops[step] = mature_eprop
        
    return fresh_eprops, mixed_eprops, mature_eprops

def label_k_score(g, co_cumsum_eprops, first=True):
    '''label_k_score
    Calculates the K score for each existing edge at each step:

    .. math::
        k_{ij, T} = \frac{1}{1 + C_{ij, T-1}}

    Args:
        g (:obj:`graph_tool.Graph`):
        co_cumsum_eprops (:obj:`dict`):

    Returns:
        k_score_eprops (:obj:`dict`):
    '''
    k_score_eprops = {}
    for i, step in enumerate(co_cumsum_eprops.keys()):
        if i == 0:
            if first:
                k_score_eprop = g.new_edge_property('float', val=1)
        else:
            k_score_eprop = g.new_edge_property('float')
            k_score_eprop.a = 1 / (co_cumsum_eprops[step_prev].a + 1)
        step_prev = step
        k_score_eprops[step] = k_score_eprop
    return k_score_eprops
            
def label_k_growth(g, co_eprops, co_cumsum_eprops):
    '''label_k_growth
    Calculates the K growth score for edges at each step.

    .. math::
        k_{ij, T} = \frac{c_{ij, T}}{1 + C_{ij, T-1}}

    Args:
        g (:obj:`graph_tool.Graph`): A graph.
        co_eprops (:obj:`dict`): Contains edge cooccurrence property maps.
        co_cumsum_eprops (:obj:`dict`): Contains edge cumulative cooccurrence
            property maps.
    Returns:
        k_growth_eprops (:obj:`dict`): Contains edge k growth property maps.
    '''
    k_growth_eprops = {}
    for i, (step, co_eprop) in enumerate(co_eprops.items()):
        if i == 0:
            continue
        else:
            k_growth_eprop = g.new_edge_property('float')
            k_growth_eprop.a = (co_eprops[step].a / 
                                    (co_cumsum_eprops[step].a 
                                    - co_eprops[step].a + 1))
            k_growth_eprops[step] = k_growth_eprop
    return k_growth_eprops
