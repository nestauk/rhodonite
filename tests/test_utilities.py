import pytest

from rhodonite.utilities import (edge_property_as_matrix, flatten, window,
        seq2cooccurrence)

def test_flatten_list_of_lists():
    nested = [
	    ['a', 'b', 'c'],
	    ['d','e'],
	    ['f']
	    ]
    flat = ['a', 'b', 'c', 'd', 'e', 'f']
    assert flatten(nested) == flat

def test_seq2cooccurrence_no_duplicates():
    doc = ['me', 'myself', 'irene']
    co_expected = sorted([
        ('irene', 'me'),
        ('irene', 'myself'),
        ('me', 'myself')
        ])
    co = sorted(seq2cooccurrence(doc))
    assert co == co_expected

def test_seq2cooccurrence_duplicates():
    doc = ['me', 'me', 'me', 'myself', 'irene']
    co_expected = sorted([
        ('irene', 'me'),
        ('irene', 'myself'),
        ('me', 'myself')
        ])
    co = sorted(seq2cooccurrence(doc))
    assert co == co_expected

def test_window_n3_no_duplicates():
    doc = ['a', 'b', 'c', 'd', 'e', 'f']
    groups_expected = [
            ('a', 'b', 'c'),
            ('b', 'c', 'd'),
            ('c', 'd', 'e'),
            ('d', 'e', 'f')
            ]
    groups = list(window(doc, n=3))
    assert groups == groups_expected

def test_window_n5_no_duplicates():
    doc = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    groups_expected = [
            ('a', 'b', 'c', 'd', 'e'),
            ('b', 'c', 'd', 'e', 'f'),
            ('c', 'd', 'e', 'f', 'g'),
            ]
    groups = list(window(doc, n=5))
    assert groups == groups_expected
