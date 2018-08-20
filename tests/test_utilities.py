import pytest

from rhodonite.utilities import edge_property_as_matrix, flatten, window

def test_flatten_list_of_lists():
    nested = [
	    ['a', 'b', 'c'],
	    ['d','e'],
	    ['f']
	    ]
    flat = ['a', 'b', 'c', 'd', 'e', 'f']
    assert flatten(nested) == flat

