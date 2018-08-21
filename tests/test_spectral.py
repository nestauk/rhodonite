import pytest
import numpy as np

from numpy.testing import assert_array_equal, assert_allclose

from rhodonite.spectral import occurrences, cooccurrences, association_strength

def test_occurrence_multi_doc(sliding_window_graph_multi_doc_no_dict_n3):
    g = sliding_window_graph_multi_doc_no_dict_n3
    o = np.array(occurrences(g).get_array())
    o_expected = np.array([3, 2, 2, 3, 2])
    assert_array_equal(o, o_expected)

def test_cooccurrence_multi_doc(sliding_window_graph_multi_doc_no_dict_n3):
    g = sliding_window_graph_multi_doc_no_dict_n3
    edges = sorted(g.edges())
    co = {}
    for e, c in zip(edges, np.array(cooccurrences(g).get_array())):
        s = int(e.source())
        t = int(e.target())
        co[(s, t)] = c

    co_expected = {
            (0, 1): 1,
            (0, 2): 2,
            (0, 3): 2,
            (0, 4): 1,
            (1, 2): 1,
            (1, 3): 2,
            (1, 4): 1,
            (2, 3): 2,
            (2, 4): 1,
            (3, 4): 2,
            }
    assert co == co_expected

def test_assoc_str_multi_doc(sliding_window_graph_multi_doc_no_dict_n3):
    g = sliding_window_graph_multi_doc_no_dict_n3
    edges = sorted(g.edges())
    a_s = np.array(sorted(association_strength(g).get_array()))

    n_cooccurrences = 15
    k = 2 * n_cooccurrences
    co_o = [1/6, 2/6, 2/9, 1/6, 1/4, 2/6, 1/4, 2/6, 1/4, 2/6,] 
    a_s_expected = np.array(sorted([k * c for c in co_o]))
   
    assert len(a_s) == len(a_s_expected)
    assert_allclose(a_s, a_s_expected)


