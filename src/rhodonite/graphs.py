from rhodonite.base import CooccurrenceGraph
from rhodonite.utilities import window, flatten, seq2cooccurrence
from rhodonite.spectral import cooccurrences, occurrences


class SlidingWindowGraph(CooccurrenceGraph):

    def __init__(self, window_size=3, **kwargs):
        """SlidingWindowGraph
        A cooccurrence graph built using a sliding window method.
        
        Args:
            **kwargs:

        Parameters:
            label2vertex (dict): A dictionary mapping the elements in sequences
                to their corresponding vertex id in the graph.
            n_edges (int): Total number of edges in the graph. Calculated when
                the graph is built.
            n_cooccurrences (int): Total number of cooccurrences in the graph.
                Calculated when the graph is built.
        """
        super().__init__(directed=False, **kwargs)
        self.window_size = window_size
        self.n_edges = 0
        self.n_cooccurrences = 0

    def sliding_window_cooccurrences(self, seqs, n):
        """sliding_window_cooccurrences
        Converts a sequments in a corpus into lists of cooccurrence tuples 
        using a sliding window method. The window moves across each 
        sequment, and on each iteration, links the term at the centre of 
        each window with the terms on either side.

        Args:
            seqs (:obj:`iter` of :obj:`iter` of :obj:`str`): A corpus of
                sequments.
            n (int): Size of the window.

        Returns:
            cooccurrences:
        """
        s_w_c = []
        for seq in seqs:
            seq_indices = range(len(seq))
            seq_indices = window(seq_indices, n=n)
            seq_indices = flatten([seq2cooccurrence(t) for t in seq_indices])
            seq = [sorted((seq[a], seq[b])) for a, b in list(set(seq_indices))]
            seq = sorted([tuple((a, b)) for a, b in seq if a !=b])
            s_w_c.append(seq)
        return s_w_c
    
    def build(self, sequences, dictionary=None):
        """build
        Builds the graph from a set of sequences and an optional dictionary.
        Overwrites parent class.

        Args:
            sequences (:obj:`iter` of :obj:`iter`):
            dictionary (dict): Dictionary mapping labels to IDs. Default is
                None.
            window_size (int):

        """
        self.sequences = sequences
        self.dictionary = dictionary

        # convert the sequence elements to integer ids and get cooccurrences
        labels = sorted(list(set(flatten(self.sequences))))
        self.n_vertices = len(labels)
        if self.dictionary is None:
            self.dictionary = {k: v for k, v in zip(labels, range(len(labels)))}
            id2vertex = {i: i for i in range(self.n_vertices)}
        else:
            id2vertex = {v: i for v, i in
                    zip(self.dictionary.values(), range(self.n_vertices))}
        
        # create a dict to access vertices by their labels
        self.label2vertex = {}
        for label, id in self.dictionary.items():
            self.label2vertex[label] = id2vertex[id]
        
        self.add_vertex(self.n_vertices)
        self.sequences_ids = [self.sequence2vertices(s) for s in self.sequences]
        o = occurrences(self)
        self.vertex_properties['occurrences'] = o

        self.sequence_cooccurrences = self.sliding_window_cooccurrences(
                self.sequences_ids,
                self.window_size
                )
        sequence_cooccurrences_flat = flatten(self.sequence_cooccurrences)
        edges = list(set(flatten(self.sequence_cooccurrences)))
        self.n_edges = len(edges)
        self.add_edge_list(edges)

        self.n_cooccurrences = len(sequence_cooccurrences_flat)
        co = cooccurrences(self)
        self.edge_properties['cooccurrences'] = co

