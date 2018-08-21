import pytest

from rhodonite.graphs import SlidingWindowGraph


@pytest.fixture()
def sliding_window_graph_single_doc_no_dict_n3():
    g = SlidingWindowGraph(window_size=3)
    seqs = [['a', 'b', 'c', 'd', 'e', 'a']]
    g.build(seqs)
    return g

@pytest.fixture()
def sliding_window_graph_multi_doc_no_dict_n3():
    g = SlidingWindowGraph(window_size=3)
    seqs = [['a', 'b', 'c', 'd'],
            ['a', 'c', 'd', 'e'],
            ['b', 'e', 'd', 'a']
            ]
    g.build(seqs)
    return g
