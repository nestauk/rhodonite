import pytest

from graph_tool import Graph

@pytest.fixture
def vertex_sequences_5():
    '''cooccurrence_sequences
    '''
    s = [
            [0, 1, 2],
            [1, 2, 3],
            [4]
            ]
    return s

