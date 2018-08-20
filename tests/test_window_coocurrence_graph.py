import pytest

from rhodonite.graphs import SlidingWindowGraph

class TestSlidingWindowGraph():

    def test_init_window_size(self):
        g = SlidingWindowGraph(window_size=5)
        assert g.window_size == 5
    
    def test_init_n_edges(self):
        g = SlidingWindowGraph()
        assert g.n_edges == 0

    def test_init_n_coocurrences(self):
        g = SlidingWindowGraph()
        assert g.n_coocurrences == 0

    def test_init_is_directed(self):
        g = SlidingWindowGraph()
        assert g.is_directed() == False

    def test_sliding_window_coocurrences_single_doc_corpus(self):
        docs = [['dog', 'cat', 'fish', 'fish_food']]
        g = SlidingWindowGraph()
        co = g.sliding_window_coocurrences(docs, 3)
        # window function should give
        # [('dog', 'cat', 'fish'), ('cat', 'fish', 'fish_food')]
        # therefore
        co_expected = [('dog', 'cat'), ('cat', 'fish'), ('fish', 'fish_food')]
        assert co == co_expected


