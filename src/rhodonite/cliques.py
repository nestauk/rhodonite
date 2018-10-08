import itertools
import os
import shutil

from collections import defaultdict
from subprocess import call


def load_cliques_cfinder(file_path):
    """load_cliques
    Loads cliques from a CFinder output file into a list of tuples.
    Args:
        file_path (str): The path to the CFinder output file. This is normally
            in a directory of outputs and named "cliques".
    Returns:
        cliques (:obj:`list` of :obj:`tuple`): A list of all of the cliques
            found by CFinder. Each clique is represented as a tuple of
            vertices.
    """
    with open(file_path, 'r') as f:
        clique_data = f.read().splitlines()
    cliques = []
    for cd in clique_data:
        if len(cd) > 0:
            if cd[0].isdigit():
                clique = cd.split(' ')[1:-1]
                clique = tuple(sorted([int(i) for i in clique]))
                cliques.append(clique)
    return cliques

def reverse_index_cliques(clique_set):
    """reverse_index_cliques
    Takes a set of network cliques and return all possible combinations of
    cliques where all cliques in a combination contain at least one common
    value.

    Args:
        clique_set (:obj:`iter` of :obj:`iter`): A set of cliques where 
            each element in the nested iterable contain vertices in the
            network.

    Returns:
        clique_union_indices (:obj:`list` of :obj:`tuple`): A list of the
            combinations of clique indices.
        clique_union_vertices (:obj:`list` of :obj:`tuple`): A list of the
            sets of vertices that comprise the clique combinations.
    """
    mapping = defaultdict(list)
    for i, cs in enumerate(clique_set):
        for vertex in cs:
            mapping[vertex].append(i)
    mapping = {k: tuple(v) for k, v in mapping.items()}
    return mapping

def clique_unions(clique_indices, limit):
    """clique_unions
    Create combinations of cliques up to limit.

    Args:
        clique_indices (:obj:`iter` of int): List of indices of cliques.
        limit (int): The maximum number of cliques in each union.

    Returns:
        combos (:obj:`iter` of :obj:`iter` of int): Tuples of clique
            combinations.
    """
    combos = []
    for l in range(1, limit):
        for combo in itertools.combinations(clique_indices, l):
            combos.append(tuple(combo))
    return combos

def is_subset(needle, haystack):
    """is_subset
    Check if needle is ordered subset of haystack.
    """

    if len(haystack) < len(needle): return False

    index = 0
    for element in needle:
        try:
            index = haystack.index(element, index) + 1
        except ValueError:
            return False
    else:
        return True

def filter_subsets(lists):
    """filter_subsets
    Given list of lists, return new list of lists without subsets.
    """

    for needle in lists:
        if not any(is_subset(needle, haystack) for haystack in lists
            if needle is not haystack):
            yield needle

