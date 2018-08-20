from rhodonite.base import CoocurrenceGraph
from rhodonite.utilities import window, flatten

class SlidingWindowGraph(CoocurrenceGraph):

    def __init__(self, directed=False, **kwargs):
        """SlidingWindowGraph
        A coocurrence graph built using a sliding window method.
        
        Args:
            sequences (:obj:`iter` of :obj:`iter`):
            dictionary (dict):
            window_size (int):
            **kwargs:

        Parameters:
            label2vertex (dict): A dictionary mapping the elements in sequences
                to their corresponding vertex id in the graph.
            n_edges (int): Total number of edges in the graph. Calculated when
                the graph is built.
            n_coocurrences (int): Total number of coocurrences in the graph.
                Calculated when the graph is built.
        """
        super().__init__(directed=directed, **kwargs)
        self.window_size = window_size
        self.n_edges = 0
        self.n_coocurrences = 0

    def sliding_window_coocurrences(seqs, n=n):
            """sliding_window_coocurrences
            Converts a sequments in a corpus into lists of coocurrence tuples 
            using a sliding window method. The window moves across each 
            sequment, and on each iteration, links the term at the centre of 
            each window with the terms on either side.
    
            Args:
                seqs (:obj:`iter` of :obj:`iter` of :obj:`str`): A corpus of
                    sequments.
                n (int): Size of the window.
    
            Returns:
                coocurrences:
            """
	s_w_c = []
	for seq in seqs:
	    seq_indices = range(len(seq))
	    seq_indices = window(seq_indices, n=n)
	    seq_indices = flatten([seq2coocurrences(t) for t in seq_indices])
	    seq = [sorted((seq[a], seq[b])) for a, b in list(set(seq_indices))]
	    seq = [tuple((a, b)) for a, b in seq if a !=b]
	    s_w_c.append(seq)
	return s_w_c
    
    def build(self, sequences, dictionary=None, window_size=3):
        self.sequences = sequences
        self.dictionary = dictionary
        self.window_size = window_size

        # convert the sequence elements to integer ids and get coocurrences
        labels = sorted(list(set(flatten(self.sequences))))
        self.n_vertices = len(labels)
        if self.dictionary is None:
            self.dictionary = {k: v for k, v in zip(labels, range(len(labels)))}
            id2vertex = {i: i for i in range(n_vertices)}
        else:
            id2vertex = {v: i for v, i in
                    zip(self.dictionary.values(), range(n_vertices))}
        
        # create a dict to access vertices by their labels
        self.label2vertex = {}
        for label, id in self.dictionary.items():
            self.label2vertex[label] = id2vertex[id]

        sequences_ids = [self.sequence2verticses(s) for s in self.sequences]
        self.ocurrences = self.calculate_occurrences()
        self.sequence_coocurrences = self.sliding_window_coocurrences(
                sequences_ids,
                self.window_size
                )
        sequence_coocurrences_flat = flatten(self.sequence_coocurrences)
        self.coocurrences = self.calculate_coocurrences(
                Counter(sequence_coocurrences_flat),
                self.n_vertices
                )
        self.n_coocurrences = len(sequence_coocurrences_flat)
        
        edges = list(set(flatten(self.sequence_coocurrences))
        # edge_n_coocurrences = self.edge_coocurrence_counts(coocurrences)
        # vertex_n_edges = self.vertex_degree_centralitiy(edge_coocurrences)
        self.n_edges = len(edges)

        # finally build the graph
        self.add_vertex(n_vertices)
        self.add_edge_list(edges)

