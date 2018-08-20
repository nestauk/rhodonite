import itertools


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

def seq2coocurrences(seq):
    """seq2coocurrences
    Converts an iterable sequence to a list of the set of tuples that 
    represent all the possible coocurrences.

    Args:
        seq (:obj:`iter`): List of elements
        dedupe(

    Returns:
        coocurrences (:obj:`list` of :obj:`tuple`): List of tuples. Each
            tuple is sorted.

    Examples:
        >>> doc = ['me', 'myself', 'irene']

        >>> seq2coocurrence(doc)
        [('me', 'myself'), ('irene', 'myself'), ('irene', 'me')]
    """
    combos = list(itertools.combinations(set(seq), r=2))
    coocurrences = list(set([tuple(sorted(c)) for c in combos]))
    return coocurrences

def window(seq, n=3):
    """window
    Yields a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...

    Args:
        seq (:obj:`iter`): Sequence to move the window across.

    Examples:
        >>> doc = ['rockets', 'and', 'moonshots', 'blame', 'it', 'on',
                'the', 'have', 'nots']

        >>> list(window(doc))
        [('rockets', 'and', 'moonshots'),
         ('and', 'moonshots', 'blame'),
         ('moonshots', 'blame', 'it'),
         ('it, 'on', 'the'),
         ('the', 'have', 'nots),
         ]
    """
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

# def edge_coocurrence_counts(coocurrences):
#     """coocurrence_counts
#     Takes a corpus of document coocurrence combinations and returns a
#     Counter object for them across the entire corpus.
#     
#     Args:
#         coocurrences (:obj:`list` of :obj:`list` of :obj:`tuple`): 
#             Corpus of documents expressed as their coocurrence pairs.
#             
#     Returns:
#         coocurrence_counts (:obj:`Counter`): Counter with keys as edges
#             and values as number of coocurrences between the two vertices.
#     """
#     coocurrence_counts = Counter(flatten(coocurrences))
#     return coocurrence_counts
# 
# def vertex_degree_centrality(coocurrence_counts):
#     """vertex_degree_centrality
#     Takes a Counter of edge coocurrences and returns the degree centrality
#     for each vertex.
#     
#     Args:
#         coocurrence_counts (:obj:`Counter`): Counter with keys as edges
#             and values as number of coocurrences between the two vertices.
#             
#     Returns: 
#         vertex_degrees (:obj:`Counter`): Counter with keys as vertices
#             and values as degree centralities.
#     """
#     vertex_degrees = Counter()
#     for vertices, count in coocurrence_counts.items():
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
# def vertex_coocurrence_centrality(coocurrence_counts):
#     """vertex_coocurrence_centrality
#     Takes a Counter of edge cooccurences and returns the coocurrence centrality
#     for each vertex. This is a summation of all coocurrences for each vertex.
#     
#     Args:
#         coocurrence_counts (:obj:`Counter`): Counter with keys as edges
#             and values as number of coocurrences between the two vertices.
#             
#     Returns:
#         vertex_coocurrences (:obj:`Counter`): Counter with keys as vertices
#             and values as number of coocurrences.
#     """
#     
#     vertex_coocurrences = Counter()
#     for vertices, count in coocurrence_counts.items():
#         v_0 = vertices[0]
#         v_1 = vertices[1]
#         if v_0 in vertex_coocurrences:
#             vertex_coocurrences[v_0] += count
#         else:
#             vertex_coocurrences[v_0] = count
#         if v_1 in vertex_coocurrences:
#             vertex_coocurrences[v_1] += count
#         else:
#             vertex_coocurrences[v_1] = count
#         self.vertex_coocurrences = vertex_coocurrences
#     return vertex_coocurrences
# 
# def edge_assocation_strength(self, n_coocurs, edge_coocurs, source_count,
#         target_count):
#     """assocation strength
#     Calculates the probabalistic assocation strength between two vertices
#     in a coocurrence network.
# 
#     Args:
#         n_coocurs (int): Total number of coocurrences in the network.
#         edge_coocurs (int): Number coocurrences between the  source and
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
#                     n_coocurrences,
#                     edge_coocurrences[tuple(sorted([s, t]))],
#                     vertex_coocurrences[s],
#                     vertex_coocurrences[t]
#                     )
#             self.edge_properties['association_strength'] = association_strength_prop
