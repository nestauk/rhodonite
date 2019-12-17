from numpy.testing import assert_almost_equal

from rhodonite.utils.math import *

def test_jaccard_similarity_set():
    a = {1, 2, 3, 4, 5}
    b = {2, 3, 4, 5}
    c = {2, 3, 4, 5, 6, 7, 8}
    expected_a_b = 4 / 5
    expected_a_c = 4 / 8
    expected_b_c = 4 / 7
    a_b = jaccard_similarity_set(a, b)
    a_c = jaccard_similarity_set(a, c)
    b_c = jaccard_similarity_set(b, c)
    assert_almost_equal(a_b, expected_a_b)
    assert_almost_equal(a_c, expected_a_c)
    assert_almost_equal(b_c, expected_b_c)
