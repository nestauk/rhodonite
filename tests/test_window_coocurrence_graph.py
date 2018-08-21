import pytest
import numpy as np
from numpy.testing import assert_array_equal
from rhodonite.graphs import SlidingWindowGraph


class TestSlidingWindowGraph():

    def test_init_window_size(self):
        g = SlidingWindowGraph(window_size=5)
        assert g.window_size == 5
    
    def test_init_n_edges(self):
        g = SlidingWindowGraph()
        assert g.n_edges == 0

    def test_init_n_cooccurrences(self):
        g = SlidingWindowGraph()
        assert g.n_cooccurrences == 0

    def test_init_is_directed(self):
        g = SlidingWindowGraph()
        assert g.is_directed() == False

    def test_sliding_window_cooccurrences_single_doc_corpus_n3(self):
        docs = [['dog', 'cat', 'fish', 'fish_food']]
        g = SlidingWindowGraph()
        co = g.sliding_window_cooccurrences(docs, 3)
        # window function should give
        # [('dog', 'cat', 'fish'), ('cat', 'fish', 'fish_food')]
        # therefore
        co_expected = [
                set([('cat', 'dog'), ('cat', 'fish'), ('cat', 'fish_food'),
                ('fish', 'fish_food'), ('dog', 'fish')])]
        co = [set(c) for c in co]
        assert co == co_expected

    def test_sliding_window_cooccurrences_multi_doc_corpus_n3(self):
        docs = [['dog', 'cat', 'fish', 'fish_food'],
                ['fox', 'rabbit', 'lettuce'],
                ['a', 'b', 'c', 'd', 'e']
                ]
        g = SlidingWindowGraph()
        co = g.sliding_window_cooccurrences(docs, 3)
        # window function should give
        # [('dog', 'cat', 'fish'), ('cat', 'fish', 'fish_food')]
        # therefore
        co_expected = [
                set([('cat', 'dog'), ('cat', 'fish'), ('cat', 'fish_food'),
                    ('dog', 'fish'), ('fish', 'fish_food')]),
                set([('fox', 'rabbit'), ('fox', 'lettuce'),
                    ('lettuce', 'rabbit')]),
                set([('a', 'b'), ('a', 'c'), ('b', 'c'), ('b', 'd'), ('c', 'd'),
                 ('c', 'e'), ('d', 'e')])
                ]
        co = [set(c) for c in co]
        assert co == co_expected
    
#     @pytest.mark.usefixtures('build_sliding_window_graph_single_doc_no_dict_n3')
    def test_build_single_doc_no_dict_n3(self, sliding_window_graph_single_doc_no_dict_n3):
        g = sliding_window_graph_single_doc_no_dict_n3
        dictionary_expected = {'a': 0, 'b': 1, 'c': 2, 'd': 3, 'e': 4}
        label2vertex_expected = {'a': 0, 'b': 1, 'c': 2, 'd': 3, 'e': 4}
        assert g.window_size == 3
        assert g.n_vertices == 5
        assert g.dictionary == dictionary_expected
        assert g.label2vertex == label2vertex_expected

        occurrences_expected = np.array([2, 1, 1, 1, 1])
        assert_array_equal(np.array(g.vp['occurrences'].get_array()),
                occurrences_expected)
        
        seq_co_expected = [[(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3),
            (2, 3), (2, 4), (3, 4)]]
        assert g.sequence_cooccurrences == seq_co_expected


