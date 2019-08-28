import pytest
import mock
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from rhodonite.cooccurrence.basic import *
from rhodonite.utils.misc import flatten


@pytest.fixture
def sequences():
    return [[0, 1, 2], [1, 2, 3], [4, 5]]

@pytest.fixture
def co_graph_undirected(sequences):
    sequences = sequences
    g, o, co = cooccurrence_graph(sequences)
    return g, o, co
    
@pytest.fixture
def co_graph_directed(sequences):
    sequences = sequences
    g, o, co = cooccurrence_graph(sequences, directed=True)
    return g, o, co

def mock_co_counts():
    co = {(0, 1): 1, (0, 2): 1, (1, 2): 2, (1, 3): 1, (2, 3): 1, (4, 5): 1}
    return co

def mock_flatten():
    flat = (0, 1, 2, 1, 2, 3, 4, 5)
    return flat

class TestBuildCooccurrenceGraphUndirected:
    
    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_num_vertices_co_graph_undirected(self, mocked_co_counts, 
            mocked_flatten, co_graph_undirected):
        g, o, co = co_graph_undirected
        assert g.num_vertices() == 6

    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_num_edges_co_graph_undirected(self, mocked_co_counts, 
            mocked_flatten, co_graph_undirected):
        g, o, co = co_graph_undirected
        assert g.num_edges() == 6

    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_o_co_graph_undirected(self, mocked_co_counts, 
            mocked_flatten, co_graph_undirected):
        g, o, co = co_graph_undirected
        o_expected = np.array([1, 2, 2, 1, 1, 1])
        assert_array_equal(o.a, o_expected)

    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_co_co_graph_undirected(self, mocked_co_counts, 
            mocked_flatten, co_graph_undirected):
        g, o, co = co_graph_undirected
        co_expected = np.array([1, 1, 2, 1, 1, 1])
        assert_array_equal(co.a, co_expected)


class TestBuildCooccurrenceGraphDirected:
    
    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_num_vertices_co_graph_directed(self, mocked_co_counts, 
            mocked_flatten, co_graph_directed):
        g, o, co = co_graph_directed
        assert g.num_vertices() == 6

    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_num_edges_co_graph_directed(self, mocked_co_counts, 
            mocked_flatten, co_graph_directed):
        g, o, co = co_graph_directed
        assert g.num_edges() == 12

    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_o_co_graph_directed(self, mocked_co_counts, 
            mocked_flatten, co_graph_directed):
        g, o, co = co_graph_directed
        o_expected = np.array([1, 2, 2, 1, 1, 1])
        assert_array_equal(o.a, o_expected)

    @mock.patch('test_basic.flatten', return_value=mock_flatten)
    @mock.patch('test_basic.cooccurrence_counts', return_value=mock_co_counts)
    def test_co_co_graph_directed(self, mocked_co_counts, 
            mocked_flatten, co_graph_directed):
        g, o, co = co_graph_directed
        co_expected = np.array([1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1 ])
        assert_array_equal(co.a, co_expected)
