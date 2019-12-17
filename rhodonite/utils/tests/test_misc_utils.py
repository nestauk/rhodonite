import pytest
from numpy.testing import assert_array_equal

from rhodonite.utils.misc import *

@pytest.fixture
def sequences():
    s = [
            [1, 2, 3],
            [0, 4, 5, 6], 
            [1, 4, 7], 
            [5, 2, 3, 8],
            [7, 3, 9, 1]
            ]
    return s

@pytest.fixture
def sequence():
    s = [1, 2, 3]
    return s

def test_reverse_index(sequences):
    expected = {
            0: [1],
            1: [0, 2, 4],
            2: [0, 3],
            3: [0, 3, 4],
            4: [1, 2],
            5: [1, 3],
            6: [1],
            7: [2, 4],
            8: [3],
            9: [4]
            }
    mapping = reverse_index(sequences)
    for k, v in mapping.items():
        assert_array_equal(v, expected[k])

def test_flatten(sequences):
    flat = list(flatten(sequences))
    expected = [1, 2, 3, 0, 4, 5, 6, 1, 4, 7, 5, 2, 3, 8, 7, 3, 9, 1]
    assert_array_equal(flat, expected)

def test_clique_unions_3(sequence):
    expected = [
            (1,), (2,), (3,),
            (1, 2), (1, 3), (2, 3),
            (1, 2, 3)
            ]
    results = clique_unions(sequence, 3)
    assert_array_equal(expected, results)

def test_recursive_combinations_3(sequence):
    expected = [
            (1,), (2,), (3,),
            (1, 2), (1, 3), (2, 3),
            (1, 2, 3)
            ]
    results = list(recursive_combinations(sequence, 3))
    assert_array_equal(expected, results)
