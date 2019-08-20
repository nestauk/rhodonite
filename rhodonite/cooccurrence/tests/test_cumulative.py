import pytest
import mock
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

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
    Edge order: (0, 1), (0, 2), (1, 2). (2, 3), (3, 4), (4, 5), (4, 6), (5, 6), (6, 7), (6, 8), (7, 8)
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
    
@pytest.fixture
def new_vertex_props(simple_cumco_graph):
    g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
    new_vprops = {}
    new_vprops[0] = g.new_vertex_property('bool')
    new_vprops[0].a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0])
    new_vprops[1] = g.new_vertex_property('bool')
    new_vprops[1].a = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0])
    new_vprops[2] = g.new_vertex_property('bool')
    new_vprops[2].a = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1])
    return new_vprops

@pytest.fixture
def new_edge_props(simple_cumco_graph):
    g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
    new_eprops = {}
    new_eprops[0] = g.new_edge_property('bool')
    new_eprops[0].a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    new_eprops[1] = g.new_edge_property('bool')
    new_eprops[1].a = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0])
    new_eprops[2] = g.new_edge_property('bool')
    new_eprops[2].a = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])
    return new_eprops

def active_edge_props(simple_cumco_graph):
    g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
    active_eprops = {}
    active_eprops[0] = g.new_edge_property('bool')
    active_eprops[0].a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    active_eprops[1] = g.new_edge_property('bool')
    active_eprops[1].a = np.array([1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0])
    active_eprops[2] = g.new_edge_property('bool')
    active_eprops[2].a = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1])
    return active_eprops

def old_edge_props(simple_cumco_graph):
    g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
    old_eprops = {}
    old_eprops[0] = g.new_edge_property('bool')
    old_eprops[0].a = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    old_eprops[1] = g.new_edge_property('bool')
    old_eprops[1].a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    old_eprops[2] = g.new_edge_property('bool')
    old_eprops[2].a = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0])
    return old_eprops

class TestBuildCumulativeCooccurrenceGraph:

    def test_num_vertices_cumulative_cooccurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        assert g.num_vertices() == 9

    def test_num_edges_cumulative_cooccurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        assert g.num_edges() == 12

    def test_steps_cumulative_cooccurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        steps_expected = np.array([0, 1, 2], dtype=np.int32)
        assert_array_equal(g.gp['steps'], steps_expected)

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
                0: np.array([1, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
                1: np.array([1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0]),
                2: np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1])
                }
        for step in g.gp['steps']:
            assert_array_equal(co[step].a, co_expected[step])

    def test_co_cumsum_cumulative_coocurrence_graph(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        co_cumsum_expected = {
                0: np.array([1, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
                1: np.array([2, 2, 3, 1, 1, 1, 1, 1, 1, 0, 0, 0]),
                2: np.array([3, 3, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1])
                }
        for step in g.gp['steps']:
            assert_array_equal(co_cumsum[step].a, co_cumsum_expected[step])

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
    
    @mock.patch('test_cumulative.label_active_edge', return_value=active_edge_props)
    @mock.patch('test_cumulative.label_new_edge', return_value=new_edge_props)
    def test_label_edge_activity_new(self, mocked_new_edge, mocked_active_edge,
            simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new, reinforcing, inactive = label_edge_activity(g, co)
        new_expected = ({
            0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])
            })
        for step in g.gp['steps']:
            assert_array_equal(new[step].a, new_expected[step])

    @mock.patch('test_cumulative.label_active_edge', return_value=active_edge_props)
    @mock.patch('test_cumulative.label_new_edge', return_value=new_edge_props)
    def test_label_edge_activity_reinforcing(self, mocked_new_edge, mocked_active_edge,
            simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        new, reinforcing, inactive = label_edge_activity(g, co)
        reinforcing_expected = ({
            0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            2: np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            })
        for step in g.gp['steps']:
            assert_array_equal(reinforcing[step].a, reinforcing_expected[step])

    @mock.patch('test_cumulative.label_active_edge', return_value=active_edge_props)
    @mock.patch('test_cumulative.label_new_edge', return_value=new_edge_props)
    def test_label_edge_activity_inactive(self, mocked_new_edge, mocked_active_edge, 
            simple_cumco_graph):
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

    def test_label_linkage_additive(self, simple_cumco_graph, new_vertex_props, 
            new_edge_props):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        additive, extending, joining = label_linkage(g, new_edge_props, 
                new_vertex_props)
        additive_expected = {
            0: np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
            }
        for step in g.gp['steps']:
            assert_array_equal(additive[step].a, additive_expected[step])

    def test_label_linkage_extending(self, simple_cumco_graph, new_vertex_props, 
            new_edge_props):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        additive, extending, joining = label_linkage(g, new_edge_props, 
                new_vertex_props)
        extending_expected = {
            0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0])
            }
        for step in g.gp['steps']:
            assert_array_equal(extending[step].a, extending_expected[step])

    def test_label_linkage_joining(self, simple_cumco_graph, new_vertex_props, 
            new_edge_props):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        additive, extending, joining = label_linkage(g, new_edge_props, 
                new_vertex_props)
        joining_expected = {
            0: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            1: np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
            2: np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            }
        for step in g.gp['steps']:
            assert_array_equal(joining[step].a, joining_expected[step])

    def test_label_k_score(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        k_score = label_k_score(g, co_cumsum)
        k_score_expected = {
            0: np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
            1: np.array([1/2, 1/2, 1/3, 1/2, 1/2, 1, 1, 1, 1, 1, 1, 1]),
            2: np.array([1/3, 1/3, 1/4, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1, 1, 1])
            }
        for step in g.gp['steps']:
            assert_array_almost_equal(k_score[step].a, k_score_expected[step])

    def test_label_k_score_no_first(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        k_score = label_k_score(g, co_cumsum, first=False)
        k_score_expected = {
            1: np.array([1/2, 1/2, 1/3, 1/2, 1/2, 1, 1, 1, 1, 1, 1, 1]),
            2: np.array([1/3, 1/3, 1/4, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1, 1, 1])
            }
        for step in g.gp['steps'][1:]:
            assert_array_almost_equal(k_score[step].a, k_score_expected[step])

    def test_label_k_growth(self, simple_cumco_graph):
        g, o, o_cumsum, co, co_cumsum = simple_cumco_graph
        k_score = label_k_growth(g, co, co_cumsum)
        k_score_expected = {
            1: np.array([1/2, 1/2, 1/3, 0, 0, 1, 1, 1, 1, 0, 0, 0]),
            2: np.array([1/3, 1/3, 1/4, 0, 0, 0, 0, 0, 0, 1, 1, 1])
            }
        for step in g.gp['steps'][1:]:
            assert_array_almost_equal(k_score[step].a, k_score_expected[step])
