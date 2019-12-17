from graph_tool import Graph
import pytest
import mock
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from rhodonite.utils.graph import *

@pytest.fixture
def graph():
    g = Graph()
    v = g.add_vertex(5)
    edges = [(0, 1, 0.1),
             (0, 2, 0.2),
             (0, 4, 0.4),
             (1, 2, 0.2),
             (1, 3, 0.3),
             (1, 4, 0.4),
             (2, 4, 0.4)]
    eprop = g.new_ep('float')
    g.add_edge_list(edges, eprops=[eprop])
    return g, eprop

def test_subgraph_eprop_values(graph):
    g, eprop = graph
    vertices = [0, 1, 2]
    vals = subgraph_eprop_values(g, vertices, eprop)
    expected = np.array([0.1, 0.2, 0.2])
    assert_array_equal(vals, expected)

# def test_subgraph_eprop_agg(graph):
#     g, eprop = graph
#     vertices = [0, 1, 2]
#     vals = subgraph_eprop_agg(g, vertices, eprop, aggfuncs=[np.median, np.mean])
#     expected = np.array([0.1, 0.2, 0.2])
#     median = np.median(expected)
#     mean = np.mean(expected)
#     assert_array_equal(vals, [median, mean])

def test_dict_to_vertex_prop(graph):
    g, eprop = graph
    vprop_dict = {0: 10, 1: 20, 2: 30, 3: 40, 4: 50}
    vprop = dict_to_vertex_prop(g, vprop_dict, 'int')
    vprop_expected = np.array([10, 20, 30, 40, 50])
    assert_array_equal(vprop.a, vprop_expected)

def test_dict_to_edge_prop(graph):
    g, _ = graph
    eprop_dict = {(0, 1): 10,
                  (0, 2): 20,
                  (0, 4): 30,
                  (1, 2): 40,
                  (1, 3): 50,
                  (1, 4): 60,
                  (2, 4): 70,
                  }
    eindex_dict = {(0, 1): 0,
                   (0, 2): 1,
                   (0, 4): 2,
                   (1, 2): 3,
                   (1, 3): 4,
                   (1, 4): 5,
                   (2, 4): 6,
                   }
    eprop = dict_to_edge_prop(g, eprop_dict, 'int', eindex_dict)
    expected = np.array([10, 20, 30, 40, 50, 60, 70])
    assert_array_equal(eprop.a, expected)

