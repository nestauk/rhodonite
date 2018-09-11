import itertools
import numpy
import warnings

from collections import Counter, defaultdict
from gensim.corpora import Dictionary
from graph_tool.all import Graph

from rhodonite.base import CooccurrenceGraph
from rhodonite.utilities import window, flatten, sequence_item_types


class CooccurrenceGraph(Graph):

    def __init__(self, *args, **kwargs):
        """SlidingWindowGraph
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
        super().__init__(directed=False, *args, **kwargs)

    def from_matrix(self, matrix, dictionary=None):
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

    def from_sequences(self, sequences, window_size=2, dictionary=None):
        """from_sequences
        Constructs a cooccurrence network from a series of sequences (an
        iterable of iterables), for example, a corpus of tokenized documents.
        Depending on the value of ``window_size``, either a sliding window is
        used to identify cooccurrences between neigbouring elements in the
        nested sequences, or cooccurrences are counted between all elements in
        each sequence.
        """
        sequences_flat = flatten(sequences)
        sequence_items = set(sequences_flat)
        self.add_vertex(len(sequence_items))

        item_types = sequence_item_types(sequence_items) 
        items_vp = self.new_vertex_property(item_types)
        
        if dictionary is not None:
            dict_item_types = sequence_item_types(dictionary.keys())
            item_values_vp = self.new_vertex_property(dict_item_types)
        
        item_vertex_mapping = {}
        for i, item in enumerate(sequence_items):
            items_vp[i] = item
            item_vertex_mapping[item] = i
            if dictionary is not None:
                item_values_vp[i] = dictionary[item]
        
        self.vp['item'] = items_vp
        if dictionary is not None:
            self.vp['item_value'] = item_values_vp

        occurrences_vp = self.occurrences(sequences, item_vertex_mapping)
        self.vertex_properties['occurrence'] = occurrences_vp
        
        length_shortest_seq = min([len(s) for s in sequences])
        if window_size > length_shortest_seq:
            warnings.warn(('Parameter window_size is larger than the shortest'
                ' sequence length. Cooccurrences will not be counted from' 
                ' sequences with lengths smaller than the window size.'))
        
        distances = defaultdict(list)
        cooccurrences = defaultdict(int)

        for co_pair, dist in self.sequences2cooccurrences(sequences, window_size):
            source = item_vertex_mapping[co_pair[0]]
            target = item_vertex_mapping[co_pair[1]]
            distances[(source, target)].append(dist)
            cooccurrences[(source, target)] += 1

        self.add_edge_list(cooccurrences.keys())

        distances_ep = self.new_edge_property('vector<int>')
        cooccurrences_ep = self.new_edge_property('int')

        for co_pair, cooccurrence in cooccurrences.items():
            cooccurrences_ep[co_pair] = cooccurrences[co_pair]
            distances_ep[co_pair] = numpy.array(distances[co_pair])
        
        self.ep['cooccurrence'] = cooccurrences_ep
        self.ep['distance'] = distances_ep

        return self

    def sequences2cooccurrences(self, sequences, window_size):
        """sequences2cooccurrences
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
        for seq in sequences:
            if window_size < 2:
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
                yield pair, dist
     
    def cooccurrence_pairs(self, seq):
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
        combos = itertools.combinations(set(seq), r=2)
        cooccurrences = set([tuple(sorted(c)) for c in combos])
        return cooccurrences

    def occurrences(self, sequences, item_vertex_mapping):
        """calculate_occurrences
        Calculates the number of times each element in the sequences occurs
        and puts them in to an array ordered by the elements' vertex IDs.

        Args:
            g (:obj:`Graph`):
            sequences (:obj:`iter` of :obj:`iter`):

        Returns
            o (:obj:`PropertyMap`): 
        """
        counts = Counter(flatten(sequences))
        o = self.new_vertex_property('int')
        for k, v in counts.items():
            o[item_vertex_mapping[k]] = v
        return o

