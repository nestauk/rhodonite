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

    def test_co_cumulative_coocurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph

        co_expected = {
                0: {(0, 1): 1,
                    (0, 2): 1,
                    (1, 2): 2,
                    (1, 3): 1,
                    (2, 3): 1,
                    (3, 4): 0,
                    (4, 5): 0,
                    (4, 6): 0,
                    (5, 6): 0,
                    (6, 7): 0,
                    (6, 8): 0,
                    (7, 8): 0,
                    }, 
                1: {(0, 1): 1,
                    (0, 2): 1,
                    (1, 2): 1,
                    (1, 3): 0,
                    (2, 3): 0,
                    (3, 4): 1,
                    (4, 5): 1,
                    (4, 6): 1,
                    (5, 6): 1,
                    (6, 7): 0,
                    (6, 8): 0,
                    (7, 8): 0,
                    }, 
                2: {(0, 1): 1,
                    (0, 2): 1,
                    (1, 2): 1,
                    (1, 3): 0,
                    (2, 3): 0,
                    (3, 4): 0,
                    (4, 5): 0,
                    (4, 6): 0,
                    (5, 6): 0,
                    (6, 7): 1,
                    (6, 8): 1,
                    (7, 8): 1,
                    }
                }
        edges = [(int(e.source()), int(e.target())) for e in g.edges()]
        for step in g.gp['steps']:
            co_dict = {e: co[step][e] for e in edges}
            assert_array_equal(co_dict, co_expected[step])

    def test_co_cumsum_cumulative_coocurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        co_cumsum_expected = {
                0: {(0, 1): 1,
                    (0, 2): 1,
                    (1, 2): 2,
                    (1, 3): 1,
                    (2, 3): 1,
                    (3, 4): 0,
                    (4, 5): 0,
                    (4, 6): 0,
                    (5, 6): 0,
                    (6, 7): 0,
                    (6, 8): 0,
                    (7, 8): 0,
                    }, 
                1: {(0, 1): 2,
                    (0, 2): 2,
                    (1, 2): 3,
                    (1, 3): 1,
                    (2, 3): 1,
                    (3, 4): 1,
                    (4, 5): 1,
                    (4, 6): 1,
                    (5, 6): 1,
                    (6, 7): 0,
                    (6, 8): 0,
                    (7, 8): 0,
                    }, 
                2: {(0, 1): 3,
                    (0, 2): 3,
                    (1, 2): 4,
                    (1, 3): 1,
                    (2, 3): 1,
                    (3, 4): 1,
                    (4, 5): 1,
                    (4, 6): 1,
                    (5, 6): 1,
                    (6, 7): 1,
                    (6, 8): 1,
                    (7, 8): 1,
                    }
                }
        edges = [(int(e.source()), int(e.target())) for e in g.edges()]
        for step in g.gp['steps']:
            co_cumsum_dict = {e: co_cumsum[step][e] for e in edges}
            assert_array_equal(co_cumsum_dict, co_cumsum_expected[step])

    def test_label_active_edges(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        active_expected = ({
            0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0]),
            2: np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1])
            })
        active = label_active_edge(g, co)
        for step in g.gp['steps']:
            assert_array_equal(active[step].a, active_expected[step])

    def test_label_new_edges_new(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new_expected = ({
            0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])
            })
        new = label_new_edge(g, co, label_old=False)
        for step in g.gp['steps']:
            assert_array_equal(new[step].a, new_expected[step])

    def test_label_new_edges_old(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        old_expected = ({
            0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            2: np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0])
            })
        new, old = label_new_edge(g, co)
        for step in g.gp['steps']:
            assert_array_equal(old[step].a, old_expected[step])

    def test_label_new_edges_steps(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        steps_expected = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2])
        new, step_prop = label_new_edge(g, co, label_old=False, label_steps=True)
        assert_array_equal(step_prop.a, steps_expected)

    def test_label_edge_activity_new(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new, reinforcing, inactive = label_edge_activity(g, co)
        new_expected = ({
            0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])
            })
        for step in g.gp['steps']:
            assert_array_equal(new[step].a, new_expected[step])

    def test_label_edge_activity_reinforcing(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new, reinforcing, inactive = label_edge_activity(g, co)
        reinforcing_expected = ({
            0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            2: np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            })
        for step in g.gp['steps']:
            assert_array_equal(reinforcing[step].a, reinforcing_expected[step])

    def test_label_edge_activity_inactive(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new, reinforcing, inactive = label_edge_activity(g, co)
        inactive_expected = ({
            0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            2: np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0])
            })
        for step in g.gp['steps']:
            assert_array_equal(inactive[step].a, inactive_expected[step])

    def test_label_new_vertex_new(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new = label_new_vertex(g, o)
        new_expected = {
                0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0]),
                1: np.array([0, 0, 0, 0, 0, 1, 1, 0, 0]),
                2: np.array([0, 0, 0, 0, 0, 0, 0, 1, 1]),
                }

        for step in g.gp['steps']:
            assert_array_equal(new[step].a, new_expected[step])

    def test_label_linkage_fresh(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new_vprops = {}
        new_vprops[0] = g.new_vertex_property('bool')
        new_vprops[0].a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0])
        new_vprops[1] = g.new_vertex_property('bool')
        new_vprops[1].a = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0])
        new_vprops[2] = g.new_vertex_property('bool')
        new_vprops[2].a = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1])
        
        new_eprops = {}
        new_eprops[0] = g.new_edge_property('bool')
        new_eprops[0].a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
        new_eprops[1] = g.new_edge_property('bool')
        new_eprops[1].a = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0])
        new_eprops[2] = g.new_edge_property('bool')
        new_eprops[2].a = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])

        additive, extending, joining = label_linkage(g, new_eprops, new_vprops)
        additive_expected = {
            0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
            }
        for step in g.gp['steps']:
            assert_array_equal(additive[step].a, additive_expected[step])

