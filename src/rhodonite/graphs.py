from gensim.corpora import Dictionary

from rhodonite.base import CooccurrenceGraph
from rhodonite.utilities import window, flatten, seq2cooccurrence
from rhodonite.spectral import cooccurrences, occurrences


class SlidingWindowGraph(CooccurrenceGraph):

    def __init__(self, sequences, dictionary=None, window_size=3, **kwargs):
        """SlidingWindowGraph
        A cooccurrence graph built using a sliding window method.
        
        Args:
            **kwargs:

        Parameters:
            label2vertex (dict): A dictionary mapping the elements in sequences
                to their corresponding vertex idx in the graph.
            n_edges (int): Total number of edges in the graph. Calculated when
                the graph is built.
            n_cooccurrences (int): Total number of cooccurrences in the graph.
                Calculated when the graph is built.
        """
        super().__init__(directed=False, **kwargs)
        self.sequences = sequences
        self.dictionary = dictionary
        self.window_size = window_size
        self.n_edges = 0
        self.n_cooccurrences = 0

    def sliding_window_cooccurrences(self, seqs, n):
        """sliding_window_cooccurrences
        Converts a sequments in a corpus into lists of cooccurrence tuples 
        using a sliding window method. The window moves across each 
        sequment, and on each iteration, links the term at the centre of 
        each window with the terms on either sidxe.

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

    def prepare(self):
        """prepare
        Creates the graph vertices and prepares dictionaries for mapping terms
        and term ids to vertices.
        """
        labels = sorted(list(set(flatten(self.sequences))))
        self.n_vertices = len(labels)
        self.add_vertex(self.n_vertices)

        if self.dictionary is  None:
            self.dictionary = Dictionary(sequences)

        vertex_tokens = self.new_vertex_property('string')
        vertex_token_idxs = self.new_vertex_property('int')
        self.token2vertex = {}
        self.tokenidx2vertex = {}
        for i, token in enumerate(labels):
            vertex_tokens[i] = token
            vertex_token_idxs[i] = self.dictionary.token2id[token]
            self.token2vertex[token] = i
            self.tokenidx2vertex[self.dictionary.token2id[token]] = i

        self.idx_sequences = [self.dictionary.doc2idx(s)
                for s in self.sequences]
        return self

    def build(self):
        """build
        Builds the coocurrence graph from the sequences. Calculates the number
        of occurrences of each term and the number of cooccurrences between
        each cooccurring term pair.
        """
        o = occurrences(self)
        self.vertex_properties['occurrences'] = o

        self.id_co_sequences = self.sliding_window_cooccurrences(
                self.idx_sequences,
                self.window_size
                )
        co_sequences_flat = [(self.tokenidx2vertex[s], self.tokenidx2vertex[t])
                for s, t in flatten(self.id_co_sequences)]
        edges = list(set(co_sequences_flat))
        self.n_edges = len(edges)
        self.add_edge_list(edges)

        self.n_cooccurrences = len(co_sequences_flat)
        co = cooccurrences(self)
        self.edge_properties['cooccurrences'] = co

