from graph_tool import Graph
from pandas.testing import assert_frame_equal
import pytest

from rhodonite.utils.tabular import *


@pytest.fixture
def graph():
    g = Graph()
    g.add_vertex(3)
    g.add_edge_list([(0, 1), (1, 2)])

    vprop_i = g.new_vp('int')
    vprop_i.a = [10, 20, 30]
    vprop_v = g.new_vp('vector<int>')
    vprop_v_ = np.array([
        [1, 2, 3], 
        [4, 5, 6],
        [7, 8, 9]
        ])
    vprop_v.set_2d_array(vprop_v_)
    vprop_s = g.new_vp('string')
    vprop_s_ = ['a', 'b', 'c']
    for s, v in zip(vprop_s_, g.vertices()):
        vprop_s[v] = s
    g.vp['v_int_prop'] = vprop_i
    g.vp['v_vec_prop'] = vprop_v
    g.vp['v_str_prop'] = vprop_s

    eprop_i = g.new_ep('int')
    eprop_i.a = np.array([100, 200])
    eprop_v = g.new_ep('vector<int>')
    eprop_v_ = np.array([
        [1, 2],
        [3, 4],
        [5, 6]
        ])
    eprop_v.set_2d_array(eprop_v_)
    eprop_s = g.new_ep('string')
    eprop_s_ = ['x', 'y']
    for s, e in zip(eprop_s_, g.edges()):
        eprop_s[e] = s
    g.ep['e_int_prop'] = eprop_i
    g.ep['e_vec_prop'] = eprop_v
    g.ep['e_str_prop'] = eprop_s

    return g

def test_vertices_to_dataframe(graph):
    g = graph
    df = vertices_to_dataframe(g, vectors=True)
    expected = pd.DataFrame({
        'v': [0, 1, 2],
        'v_int_prop': [10, 20, 30],
        'v_vec_prop': [
            np.array([1, 4, 7]),
            np.array([2, 5, 8]),
            np.array([3, 6, 9])
            ],
        # 'v_str_prop': ['a', 'b', 'c'] feature not added yet
        })
    assert_frame_equal(df, expected)

def test_edges_to_dataframe(graph):
    g = graph
    df = edges_to_dataframe(g, vectors=True)
    expected = pd.DataFrame({
        's': [0, 1],
        't': [1, 2],
        'e_index': [0, 1],
        'e_int_prop': [100, 200],
        'e_vec_prop': [
            np.array([1, 3, 5]),
            np.array([2, 4, 6]),
            ],
        # 'e_str_prop': ['x', 'y'] feature not added yet
        })
    assert_frame_equal(df, expected)

