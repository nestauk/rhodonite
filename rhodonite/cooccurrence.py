import itertools
import numpy
import warnings

from collections import Counter, defaultdict
from gensim.corpora import Dictionary
from graph_tool.all import Graph
from itertools import chain, combinations

from rhodonite.utilities import window, flatten, sequence_item_types


class CooccurrenceGraph(Graph):

    def __init__(self, directed=False, *args, **kwargs):
        """CooccurrenceGraph
        A graph class that has methods for contructing cooccurrence networks.

        The resulting network is one where vertices are a series of items, and
        the edges represent there being at least one cooccurence between a pair
        of items. An internalised ``PropertyMap`` provides the edge cooccurrence
        value.
        
        Args:
            *args
            **kwargs

        Attributes:
            vp.occurrence (int):
            vp.item:
            vp.item_value:
            ep.cooccurrence (int):
            ep.distance (vector<int>):
        """
        super().__init__(directed=directed, *args, **kwargs)

    def from_matrix(self, matrix, dictionary):
        """from matrix
        Constructs a coocurrence network from a cooccurrence matrix with N
        columns and rows, where N is the total number of unique items. The
        value of the cooccurrence matrix is zero at all points other than
        those where a cooccurrence occurs. The value at point (i, j) is the
        number of cooccurrences between the items i and j.

        Args:
            matrix:
            dictionary:

        Returns:
            self:
        """
        pass

    def from_sequences(self, sequences, dictionary, label_items=True, 
            item_type='string', distance_agg=None):
        """from_sequences
        Constructs a cooccurrence network from a series of sequences (an
        iterable of iterables), for example, a corpus of tokenized documents.
        Depending on the value of ``window_size``, either a sliding window is
        used to identify cooccurrences between neigbouring elements in the
        nested sequences, or cooccurrences are counted between all elements in
        each sequence.

        Args:
            sequences (:obj:`iter` of :obj:`iter` :obj:`iter`):
            dictionary (dict):
            window_size (int):
            distance_agg (function): 

        Returns:
            self
        """
        num_items = len(dictionary.keys())
        self.add_vertex(num_items)
        
        if label_items:
            items_vp = self.new_vertex_property(item_type)
            
            for i, item in dictionary.items():
                items_vp[i] = item
            
            self.vp['item'] = items_vp

        occurrences_vp = self.occurrences(sequences, dictionary)
        self.vertex_properties['occurrence'] = occurrences_vp

        cooccurrences = self.generate_cooccurrences(sequences)

        if self.is_directed():
            cooccurrences.update({k[::-1]: v for k, v in cooccurrences.items()})

        cooccurrences_ep = self.new_edge_property('int')
        self.ep['cooccurrence'] = cooccurrences_ep

        self.add_edge_list(
                ((k[0], k[1], v) for k, v in cooccurrences.items()),
                eprops=[cooccurrences_ep]
                )

        return self
    
    def generate_cooccurrences(self, sequences):
        co_pairs = (itertools.combinations(sorted(s), 2) for s in sequences)
        co = Counter(flatten(co_pairs))
        return co

    def seqs2cooccurrences(self, sequences, window_size):
        """seqs2cooccurrences
        From an iterable of sequences, this function generates cooccurrence
        pairs, and the distances between the pair elements, as they are found
        in the sequences.

        Args:
            sequences (:obj:`iter` of :obj:`iter`): An iterable containing
                sequences of elements, for which the cooccurrences are to be
                found.
            window_size (int): The window size that slides over each sequence,
                determining the maximum distance between elements that can be
                condsidered for cooccurrence.

        Yields:
            co_pair (:obj:`tuple`): A pair of cooccurring elements.
            dist (int): The distance between them in the original sequence.
        """
        distances = defaultdict(list)
        cooccurrences = defaultdict(int)
        for seq in sequences:
            if window_size < 2:
                n = len(seq)
            elif len(seq) < window_size:
                n = len(seq)
            else:
                n = window_size
            idx_windows = window(range(len(seq)), n=n)
            idx_pairs = flatten(
                    [self.cooccurrence_pairs(iw) for iw in idx_windows]
                    )
            seq_dists = [abs(a - b) for a, b in set(idx_pairs)]
            seq_pairs = [(seq[a], seq[b]) for a, b in set(idx_pairs)]
            co_pair_dists = [(sp, sd) for sp, sd in zip(seq_pairs, seq_dists)
                    if sp[0] != sp[1]]
            for pair, dist in co_pair_dists:
                pair = tuple(sorted(pair))
                distances[pair].append(dist)
                cooccurrences[pair] += 1
        return cooccurrences, distances
     
    def cooccurrence_pairs(self, seq):
        """coocurrence_pairs
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
        combos = itertools.combinations(set(seq), r=2)
        cooccurrences = set([tuple(sorted(c)) for c in combos])
        return cooccurrences

    def occurrences(self, sequences, dictionary):
        """calculate_occurrences
        Calculates the number of times each element in the sequences occurs
        and puts them in to an array ordered by the elements' vertex IDs.

        Args:
            sequences (:obj:`iter` of :obj:`iter`):

        Returns
            o (:obj:`PropertyMap`): 
        """
        counts = Counter(flatten(sequences))
        o = self.new_vertex_property('int')
        for k in dictionary.keys():
            o[k] = counts.get(k, 0)
        return o

def cooccurrence_counts(sequences):
    '''cooccurrence_counts
    Counts the cooccurring pairs from a nested sequence.

    Args:
        sequences (:obj:`iter` of :obj:`iter`):

    Returns:
        (:obj:`collections.Counter`):
        
    '''
    combos = (combinations(sorted(sequence), 2) for sequence in sequences)
    return Counter(chain(*combos))

def occurrence_counts(sequences):
    '''occurrence_counts
    Counts the occurrences of each element in a nested sequence.

    Args:
        sequences (:obj:`iter`):

    Returns
        (:obj:`collections.Counter`)
    '''
    return Counter(chain(*sequences))

