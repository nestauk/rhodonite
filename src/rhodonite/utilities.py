import itertools
import numpy as np
import pathlib


def get_vp_values(g, vertex_prop_name):
    """get_vp_values
    Retrieves a vertex property from a graph, taking into account any filter.
    
    Args:
        g (:obj:`Graph`): A graph.
        vertex_prop_name (str): The name of an internal vertex property.
        
    Returns:
        pm (:obj:`PropertyMapArray`): An array of the property map.
    """
    mask = g.get_vertex_filter()[0]
    if mask is not None:
        mask = np.where(mask.get_array())
        pm = g.vp[vertex_prop_name].get_array()[mask]
    else:
        pm = g.vp[vertex_prop_name].get_array()
    
    return pm

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

    mat = np.zeros((g.n_vertices, g.n_vertices))
    prop = g.edge_properties[prop_name]
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
    flat = [item for iter in list_of_iters for item in iter]
    return flat

def seq2cooccurrence(seq):
    """seq2cooccurrence
    Converts an iterable sequence to a list of the set of tuples that 
    represent all the possible cooccurrences.

    Args:
        seq (:obj:`iter`): List of elements
        dedupe(

    Returns:
        cooccurrences (:obj:`list` of :obj:`tuple`): List of tuples. Each
            tuple is sorted.

    Examples:
        >>> doc = ['me', 'myself', 'irene']

        >>> sorted(seq2cooccurrence(doc))
        [('irene', 'me'), ('irene', 'myself'), ('me', 'myself')]
    """
    combos = list(itertools.combinations(set(seq), r=2))
    cooccurrences = list(set([tuple(sorted(c)) for c in combos]))
    return cooccurrences

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


#     def get_occurrence(self, label):
#         """get_occurrences
#         Returns the number of occurrences for a label from the occurrences
#         vector.
# 
#         Args:
#             label: A label from the original set of sequences.
# 
#         Returns:
#             The number of times the label occurred within the sequences.
#         """
# 
#         return self.occurences[self.label2vertex[label]]


#     def get_cooccurrence(self, label_0, label_1):
#         """get_edge_cooccurrence
#         Get the number of cooccurrences between two labels.
# 
#         Args:
#             label_0: A label from the sequence.
#             label_1: A label from the sequences
# 
#         Returns:
#             co: The number of times the two labels coocur.
#         """
#         i_0 = self.label2vertex[label_0]
#         i_1 = self.label2_vertex[label_1]
#         co = self.cooccurrences[i_0][i_1]
#         return co

# def edge_cooccurrence_counts(cooccurrences):
#     """cooccurrence_counts
#     Takes a corpus of document cooccurrence combinations and returns a
#     Counter object for them across the entire corpus.
#     
#     Args:
#         cooccurrences (:obj:`list` of :obj:`list` of :obj:`tuple`): 
#             Corpus of documents expressed as their cooccurrence pairs.
#             
#     Returns:
#         cooccurrence_counts (:obj:`Counter`): Counter with keys as edges
#             and values as number of cooccurrences between the two vertices.
#     """
#     cooccurrence_counts = Counter(flatten(cooccurrences))
#     return cooccurrence_counts
# 
# def vertex_degree_centrality(cooccurrence_counts):
#     """vertex_degree_centrality
#     Takes a Counter of edge cooccurrences and returns the degree centrality
#     for each vertex.
#     
#     Args:
#         cooccurrence_counts (:obj:`Counter`): Counter with keys as edges
#             and values as number of cooccurrences between the two vertices.
#             
#     Returns: 
#         vertex_degrees (:obj:`Counter`): Counter with keys as vertices
#             and values as degree centralities.
#     """
#     vertex_degrees = Counter()
#     for vertices, count in cooccurrence_counts.items():
#         v_0 = vertices[0]
#         v_1 = vertices[1]
#         if v_0 in vertex_degrees:
#             vertex_degrees[v_0] += 1
#         else:
#             vertex_degrees[v_0] = 1
#             
#         if v_1 in vertex_degrees:
#             vertex_degrees[v_1] += 1
#         else:
#             vertex_degrees[v_1] = 1
# 
#     return vertex_degrees
# 
# def vertex_cooccurrence_centrality(cooccurrence_counts):
#     """vertex_cooccurrence_centrality
#     Takes a Counter of edge cooccurences and returns the cooccurrence centrality
#     for each vertex. This is a summation of all cooccurrences for each vertex.
#     
#     Args:
#         cooccurrence_counts (:obj:`Counter`): Counter with keys as edges
#             and values as number of cooccurrences between the two vertices.
#             
#     Returns:
#         vertex_cooccurrences (:obj:`Counter`): Counter with keys as vertices
#             and values as number of cooccurrences.
#     """
#     
#     vertex_cooccurrences = Counter()
#     for vertices, count in cooccurrence_counts.items():
#         v_0 = vertices[0]
#         v_1 = vertices[1]
#         if v_0 in vertex_cooccurrences:
#             vertex_cooccurrences[v_0] += count
#         else:
#             vertex_cooccurrences[v_0] = count
#         if v_1 in vertex_cooccurrences:
#             vertex_cooccurrences[v_1] += count
#         else:
#             vertex_cooccurrences[v_1] = count
#         self.vertex_cooccurrences = vertex_cooccurrences
#     return vertex_cooccurrences
# 
# def edge_assocation_strength(self, n_coocurs, edge_coocurs, source_count,
#         target_count):
#     """assocation strength
#     Calculates the probabalistic assocation strength between two vertices
#     in a cooccurrence network.
# 
#     Args:
#         n_coocurs (int): Total number of cooccurrences in the network.
#         edge_coocurs (int): Number cooccurrences between the  source and
#             target vertices. 
#         source_count (int): Total number of occurences of the source
#             vertex.
#         target_count (int): Total number of occurrences of the target
#             vertex.
# 
#     Returns:
#         a_s (float): Calculated association strength between the source and
#             target vertices.
#     """
# 
#     a_s = (2 * n_coocurs * edge_coocurs) / (source_count * target_count)
#     return a_s

#     def edge_association_strengths(self):
#         association_strength_prop = self.new_edge_property("float")
#         for s, t in edges:
#             association_strength_prop[self.edge(s, t)] = assoc_strength(
#                     n_cooccurrences,
#                     edge_cooccurrences[tuple(sorted([s, t]))],
#                     vertex_cooccurrences[s],
#                     vertex_cooccurrences[t]
#                     )
#             self.edge_properties['association_strength'] = association_strength_prop
