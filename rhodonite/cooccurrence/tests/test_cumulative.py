import pytest
import numpy as np
from numpy.testing import assert_array_equal

from rhodonite.cooccurrence.cumulative import *

@pytest.fixture
def simple_cumco_graph():
    '''vertex_steps_sequences
    
    Step 0:
        Existing vertices: []
        New vertices: [0, 1, 2, 3, 4]
        Inactive edges: []
        Reinforcing edges: []
        New edges:
            Mature: []
            Mixed: []
            Fresh: [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    Step 1:
        Existing vertices: [0, 1, 2, 3, 4]
        New vertices: [5, 6]
        Inactive edges: [(1, 3), (2, 3)]
        Reinforcing edges: [(0, 1), (0, 2), (1, 2)]
        New edges:
            Mature: [(3, 4)]
            Mixed: [(4, 5), (4, 6)]
            Fresh: [(5, 6)]

    Step 2:
        Existing vertices: [0, 1, 2, 3, 4, 5, 6]
        New vertices: [7, 8]
        Inactive edges: [(1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 6)]
        Reinforcing edges: [(0, 1), (0, 2), (1, 2)]
        New edges:
            Mature:
            Mixed: [(6, 7), (6, 8)]
            Fresh: [(7, 8)]
    '''

    steps = [0, 1, 2]
    sequences = [
            [
                [0, 1, 2],
                [1, 2, 3],
                [4]
                ],
            [
                [0, 1, 2],
                [3, 4],
                [4, 5, 6]
                ],
            [
                [0, 1, 2],
                [6, 7, 8]
                ]
            ]
    g, o, o_cumsum, co, co_cumsum = cumulative_cooccurrence_graph(
                steps, sequences)
    return g, o, o_cumsum, co, co_cumsum 

class TestBuildCumulativeCooccurrenceGraph:

    def test_num_vertices_cumulative_cooccurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        assert g.num_vertices() == 9

    def test_num_edges_cumulative_cooccurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        assert g.num_edges() == 12

    def test_o_cumulative_coocurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        o_expected = {
                0: np.array([1, 2, 2, 1, 1, 0, 0, 0, 0]), 
                1: np.array([1, 1, 1, 1, 2, 1, 1, 0, 0]), 
                2: np.array([1, 1, 1, 0, 0, 0, 1, 1, 1])}

        for step in g.gp['steps']:
            assert_array_equal(o[step].a, o_expected[step])

    def test_o_cumsum_cumulative_coocurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        o_cumsum_expected = {
                0: np.array([1, 2, 2, 1, 1, 0, 0, 0, 0]), 
                1: np.array([2, 3, 3, 2, 3, 1, 1, 0, 0]), 
                2: np.array([3, 4, 4, 2, 3, 1, 2, 1, 1])}

        for step in g.gp['steps']:
            assert_array_equal(o_cumsum[step].a, o_cumsum_expected[step])

