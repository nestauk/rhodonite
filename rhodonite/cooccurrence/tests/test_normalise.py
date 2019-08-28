from graph_tool import Graph
import pytest
import mock
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from rhodonite.cooccurrence.normalise import *

@pytest.fixture
def co_graph_directed():
    '''co_graph_directed
    '''
    g = Graph(directed=True)
    g.add_vertex(2)
    edges = [(0, 1), (1, 0), (0, 2), (2, 0), (1, 2), (2, 1)]
    g.add_edge_list(edges)
    o = g.new_vertex_property('int')
    o.a = np.array([3, 4, 2])
    co = g.new_edge_property('int')
    co.a = np.array([2, 2, 1, 1, 2, 2])
    return g, o, co

@pytest.fixture
def co_graph_undirected():
    '''co_graph_undirected
    '''
    g = Graph(directed=False)
    g.add_vertex(2)
    edges = [(0, 1), (0, 2), (1, 2)]
    g.add_edge_list(edges)
    o = g.new_vertex_property('int')
    o.a = np.array([3, 4, 2])
    co = g.new_edge_property('int')
    co.a = np.array([2, 1, 2])
    return g, o, co

def test_association_strength_directed(co_graph_directed):
    g, o, co = co_graph_directed
    a_s = association_strength(g, o, co)
    a_s_expected = np.array([5/3, 5/3, 5/3, 5/3, 5/2, 5/2])
    assert_array_almost_equal(a_s.a, a_s_expected)

def test_association_strength_undirected(co_graph_undirected):
    g, o, co = co_graph_undirected
    a_s = association_strength(g, o, co)
    a_s_expected = np.array([5/3, 5/3, 5/2])
    assert_array_almost_equal(a_s.a, a_s_expected)

def test_conditional_probability_directed(co_graph_directed):
    g, o, co = co_graph_directed
    c_p = conditional_probability(g, o, co)
    c_p_expected = np.array([1/2, 2/3, 1/2, 1/3, 1, 1/2])
    assert_array_almost_equal(c_p.a, c_p_expected)

def test_conditional_probability_undirected(co_graph_undirected):
    g, o, co = co_graph_undirected
    c_p_source, c_p_target = conditional_probability(g, o, co)
    c_p_s_expected = np.array([1/2, 1/2, 1])
    c_p_t_expected = np.array([2/3, 1/3, 1/2])
    assert_array_almost_equal(c_p_source.a, c_p_s_expected)
    assert_array_almost_equal(c_p_target.a, c_p_t_expected)
