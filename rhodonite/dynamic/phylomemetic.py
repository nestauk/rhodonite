import logging
import numpy as np
import pandas
import math
import os

from collections import defaultdict
from graph_tool.all import Graph, GraphView, edge_endpoint_property
from graph_tool.topology import all_predecessors, shortest_distance
from itertools import repeat, combinations
from operator import itemgetter
from sklearn.preprocessing import MultiLabelBinarizer
from multiprocessing import Pool, cpu_count

from rhodonite.utils.misc import (window, flatten, clear_graph, 
        get_aggregate_vp, reverse_index_communities, clique_unions, 
        reverse_index, recursive_combinations)
from rhodonite.utils.math import jaccard_similarity, jaccard_similarity_set
from rhodonite.utils.tabular import vertices_to_dataframe

import time

logger = logging.getLogger(__name__)

def label_ages(g):
    """label_ages
    Get the ages of each vertex in a graph and put them into
    a property map. The age is defined as the number of steps between
    a vertex and its most distant antecedent.
    
    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
        
    Returns
    -------
        age_vp : :obj:`graph_tool.VertexPropertyMap` 
            A property map containing the age of each vertex.
    """
    age_vp = g.new_vertex_property('int')
    for v in g.vertices():
        if v.out_degree() == 0:
            age_vp[v] = 0
        else:
            parent_ages = [age_vp[p] for p in v.out_neighbours()]
            max_age = np.max(parent_ages)
            age_vp[v] = max_age + 1
    return age_vp

def community_density(community, co_graph, density_prop, fill=0):
    """community_density
    Calculate the density of a clique based on the number of occurrences
    and coocurrences of the terms within it. Based on Callon et al. 1991.

    Parameters
    ----------
        community : :obj:`iter` of :obj:`int`: 
            A set of terms that comprise a single clique.
        g : :obj:`graph_tool.Graph` 
            The coocurrence graph from which the clique originated
        fill :obj:`float`: 
            A number to fill in for the cooccurrence value if none exists.
    Returns
    -------
        density :obj:`float`: The density of the clique.
    """
    card = len(community)
    densities = density_prop.a[community]
    density = 1 / card * np.sum(densities)
    return density

def label_density(g, cooccurrence_graphs, norm=None):
    """label_density
    Creates a property map with the density of each vertex based on the items
    contained within the community it represents.
    Requires a cooccurrence graphs.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
        cooccurrence_graphs : :obj:`iter` of :obj:`Graph` 
            A list of cooccurrence graphs for each time period.
        norm :obj:`function`: 
            A normalisation function.

    Returns
    -------
        community_densities : :obj:`PropertyMap` A property map containing the
            densities of the phylomemetic graph vertices.
    """
    community_densities = g.new_vertex_property('float')

    time_steps = sorted(np.unique(g.vp['label'].a))
    
    co_time_step = None
    for v in g.vertices():
        ts = g.vp['label'][v]
        if ts != co_time_step:
            co_time_step = ts
            co = cooccurrence_graphs[ts]
            o_source = edge_endpoint_property(co, co.vp['occurrence'], 'source') 
            o_target = edge_endpoint_property(co, co.vp['occurrence'], 'target')

            density_prop = co.new_edge_property('float')
            density_prop.a = (
                    np.square(co.ep['cooccurrence'].a) / 
                    (o_source.a * o_target.a)
                    )
        community_densities[v] = community_density(
                g.vp['item'][v], 
                co,
                density_prop
                )
    return community_densities

def label_emergence(g):
    """label_emergence
    Creates a property map that identifies whether a vertex is categorised as
    ephemeral, emerging, steady, or declining.

    The categories are represented as integers with the following mapping:
        {0: 'ephemeral',
         1: 'emerging',
         2: 'steady',
         3: 'declining'}
         
    Parameters
    ----------
        g : :obj:`Graph`
            A graph.

    Returns
    -------
        emergence : :obj:`graph_tool.VertexPropertyMap` 
            An integer property map that represents the emergence category of 
            each vertex.
    """
    emergence = g.new_vertex_property('int')
    for v in g.vertices():
        if (v.in_degree() == 0) & (v.out_degree() == 0):
            # ephemeral
            emergence[v] = 0
        elif (v.out_degree() == 0) & (v.in_degree() > 0):
            # emerging
            emergence[v] = 1
        elif (v.out_degree() > 0) & (v.in_degree() == 0):
            # declining
            emergence[v] = 3
        else:
            # steady
            emergence[v] = 2
    return emergence

def label_special_events(g):
    """label_special_events
    Creates property maps that identify whether a vertex belongs to the
    possible categories of "special events": branching or merging.
    
    Parameters
    ----------
        g : :obj:`Graph` 
            A graph.

    Returns
    -------
        branching : :obj:`graph_tool.VertexPropertyMap` 
            A boolean property map that is True where a vertex's out degree is 
            greater than 1.
        merging : :obj:`graph_tool.VertexPropertyMap` 
            A boolean property map that is True where a vertex's in degree is 
            greater than 1.
    """
    branching = g.new_vertex_property('bool')
    merging = g.new_vertex_property('bool')
    for v in g.vertices():
        if v.out_degree() >= 2:
            merging[v] = True
        if v.in_degree() >= 2:
            branching[v] = True
    return branching, merging

def label_cross_pollination(g, merging_prop, agg=np.mean):
    """label_cross_pollination
    Cross-pollination is defined as the average pairwise Jaccard similarity
    between the parents of a merging node.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
        merging_prop : :obj:`graph_tool.VertexPropertyMap` 
            A vertex property map that is True where a vertex is a merging node.
        agg: :obj:`function` 
            An aggregation function. Typically mean or median.

    Returns
    -------
        cross_pol_prop : :obj:`graph_tool.VertexPropertyMap`
            Contains cross pollination values of each vertex.
    """
    cross_poll_prop = g.new_vertex_property('float')
    g_merging = GraphView(g, vfilt=merging_prop, skip_efilt=True)
    for v in g_merging.vertices():
        parents = [g.vp['item'][p] for p in g.vertex(v).out_neighbors()]
        jaccard = agg(
            [1 - jaccard_similarity_set(list(c[0]), list(c[1]))
            for c in combinations(parents, 2)]
        )
        cross_poll_prop[v] = jaccard
    return cross_poll_prop

def label_diversification(g, branching_prop, agg=np.mean):
    """label_diversification
    Diversification is defined as the average pairwise Jaccard similarity
    between the children of a branching node.

    Parameters
    ----------
        g : :obj:`Graph` 
            A graph.
        branching_prop : :obj:`graph_tool.VertexPropertyMap` 
            A verted property map that is True where a vertex is a branching 
            node.
        agg: :obj:`function` 
            An aggregation function. Typically mean or median.
    """
    diversification_prop = g.new_vertex_property('float')
    g_branching = GraphView(g, vfilt=branching_prop, skip_efilt=True)
    for v in g_branching.vertices():
        children = [g.vp['item'][c] for c in g.vertex(v).in_neighbors()]
        jaccard = agg(
            [1 - jaccard_similarity_set(list(c[0]), list(c[1]))
            for c in combinations(children, 2)]
        )
        diversification_prop[v] = jaccard
    return diversification_prop

def label_item_emergence(g):
    """label_item_emergence
    Creates property maps that label whether a vertex contains an item that is
    emerging. An item is defined as emerging if it appears for the first time in
    the network.

    Parameters
    ----------
        g : :obj:`graph_tool.Graph` 
            A graph.
    Returns
    -------
        item_emergence_count : :obj:`graph_tool.VertexPropertyMap` 
            Property map with values representing the number of emerging items 
            in a vertex.
        item_emergence_item : :obj:`graph_tool.VertexPropertyMap` 
            Property map with a list of emerging items for each vertex.
    """
    item_emergence_count = g.new_vertex_property('int')
    item_emergence_item = g.new_vertex_property('vector<int>')

    for term, vertices in reverse_mapping(g).items():
        # find first time item appears
        v_0 = sorted(vertices)[0]
        label = g.vp['label'][v_0]
        for v in vertices:
            if g.vp['label'][v] == label:
                item_emergence_count[v] += 1
                item_emergence_item[v].append(term)
    return item_emergence_count, item_emergence_item

def label_item_inheritance(g, item_emergence_item_prop):
    """label_item_inheritance
    Creates property maps that labels vertices as containing recombining or
    reconducting items. An item is recombining for a community if it exists
    previously in the phylomemetic network, but not in a direct ancestor of 
    the community. An item is reconducting for a community if it exists in a
    direct parent of the community.
    
    Parameters
    ----------
        g : :obj:`graph_tool.Graph` A graph.
        item_emergence_item_prop : :obj:`graph_tool.VertexPropertyMap` 
            Items that are emerging as calculated by label_item_emergence.
    Returns
    -------
        item_recombination_count : :obj:`graph_tool.VertexPropertyMap` 
            Property map with values representing the number of recombining 
            items in a vertex.
        item_recombination_item : :obj:`graph_tool.VertexPropertyMap` 
            Property map with a list of recombining items for each vertex.
        item_reconduction_count : :obj:`graph_tool.VertexPropertyMap` 
            Property map with values representing the number of reconducting 
            items in a vertex.
            
        item_reconduction_item : :obj:`graph_tool.PropertyMap` Property map
            with a list of reconducting items for each vertex.
    """
    item_reconduction_count = g.new_vertex_property('int')
    item_reconduction_item = g.new_vertex_property('vector<int>')
    
    item_recombination_count = g.new_vertex_property('int')
    item_recombination_item = g.new_vertex_property('vector<int>')
    
    vertex_parents = {}
    vertex_predecessors = {}
    
    for term, vertices in reverse_mapping(g).items():
        for v in vertices:
            if term not in item_emergence_item_prop[v]:
                if v not in vertex_parents:
                    dist_map, pred_map, pred_visited = shortest_distance(
                        g,
                        source=v,
                        pred_map=True,
                        return_reached=True,
                    )
                    parents = set([g.vertex(p) for p in pred_visited])
                    vertex_parents[v] = parents
                else:
                    parents = vertex_parents[v]

                if v not in vertex_predecessors:
                    label = g.vp['label'][v]
                    predecessors = set([p for p in vertices if g.vp['label'][p] < label])
                    vertex_predecessors[v] = predecessors
                else:
                    predecessors = vertex_predecessors[v]

                overlap = predecessors.intersection(parents)
                n_overlap = len(overlap)

                if n_overlap > 0:
                    item_reconduction_item[v].append(term)
                    item_reconduction_count[v] += 1
                elif n_overlap == 0:
                    item_recombination_item[v].append(term)
                    item_recombination_count[v] += 1
    return (
        item_recombination_count, 
        item_recombination_item, 
        item_reconduction_count, 
        item_reconduction_item
    )

def reverse_mapping(g):
    """_reverse_mapping
    Returns the reverse mapping of items to vertices for the graph. If it
    does not exist, it is created.

    Returns
    -------
        self.reverse_mapping (dict): A mapping of community items to sets
            of vertices.
    """
    if not hasattr(g, 'reverse_mapping'):
        reverse_mapping = defaultdict(set)

        for v in g.vertices():
            items = g.vp['item'][v]
            for item in items:
                reverse_mapping[item].add(v)
        g.reverse_mapping = reverse_mapping

    return g.reverse_mapping

def find_links(args):
    """find_links
    Finds the inter-temporal links between a clique at time period and cliques
    occurring at previous time periods.

    Parameters
    ----------
        args : :obj:`tuple`, Contains:
        cf : :obj:`iter` of :obj:`int`: 
            The clique for which parents are being identified.
        cfi :obj:`int`: 
            The clique index relative to all cliques in the phylomemy.
        cps : :obj:`iter` of :obj:`iter` of :obj:`int`: 
            The cliques from the immediately previous time period to cf.
        pos : :obj:`tuple` 
            The start and end vertices in the previous time period.
        cpm : :obj:`scipy.sparse` 
            Binarized sparse matrix of communties from previous time period.
        mapping : :obj:`dict` 
            Maps community elements to communities.
        binarizer : :obj:`MultiLabelBinarizer` 
            A multi label binarizer fitted to the entire vocabulary across all 
            cliques.
        parent_limit :obj:`int`: 
            The limit for the size of possible combinations of potential parent 
            cliques that will be considered.
        threshold :obj:`float`: 
            The minimum percentile for the Jaccard similartiy between parents 
            and children, below which parent candidates will be discarded.
    
    Returns
    -------
        links : :obj:`list` 
            A list that contains the inter-temporal edges that have been found 
            as well as the corresponding Jaccard indexes. Each element is of the
            format: 
            `((source, target), group_jaccard_index, single_jaccard_index)`.
            `group_jaccard_index` is the Jaccard index of any union of parents
            that has been identified as the most likely antecedants of the and
            single_jaccard_index is the individual Jaccard similarity between
            each potential parent and the child.
    """
    cf, cfi, cps, pos, mapping, parent_limit = args

    start_p, end_p = pos
    start_f = end_p
    links = []
    possible_parents = []
    cp_indexes = []
    for element in cf:
        element_matches = list(mapping[element])
        if len(element_matches) == 0:
            continue
        else:
            cp_matches = [cps[i] for i in element_matches]
            j = [jaccard_similarity_set(cp, cf) for cp in cp_matches]
            if math.isclose(np.max(j), 1):
                direct_parent = np.nonzero([math.isclose(ji, 1) for ji in j])[0][0]
                parent_cp = element_matches[direct_parent]
                links.append(
                       (cfi + start_f, parent_cp + start_p, 1, 1)
                        )
                return links

            close_matches = np.nonzero(j == np.max(j))[0]
            close_match_indexes = [element_matches[i] for i in close_matches]
            possible_parents += close_match_indexes
    
    cp_indexes = set(possible_parents)
     
    if len(cp_indexes) > 0:
        cp_union_indices = clique_unions(cp_indexes, parent_limit)

        cp_union_vertices = []
        for cui in cp_union_indices:
            cuv = set(flatten([cps[i] for i in cui]))
            cp_union_vertices.append(cuv)
        j = [jaccard_similarity_set(cpu, cf) for cpu in cp_union_vertices]
        cp_union_j_mapping = {k[0]: v for k, v in zip(cp_union_indices, j) if len(k) == 1}
        j_max = np.max(j)
        parent_clique_indices = np.nonzero(j == j_max)[0]
        if len(parent_clique_indices) > 0:
            parent_cliques = itemgetter(*parent_clique_indices)(cp_union_indices)
            j_parents = [np.max(j_p) for j_p in j if np.max(j_p) > 0]
            if any(isinstance(i, tuple) for i in parent_cliques):
                parent_cliques = set(flatten(parent_cliques))

            for pc in parent_cliques:
                link = (
                        cfi + start_f, 
                        pc + start_p,
                        j_max,
                        cp_union_j_mapping[pc]
                    )
                links.append(link)
    return links

def filter_by_size(sequences, min_size=None, max_size=None):
    """filter_communities
    Returns the communities in a series that are larger than the minimum
    community size.

    Parameters
    ----------
        cliques_set : :obj:`iter` of :obj:`iter` 
            Series of iterables containing the vertices that make up cliques.
        min_size :obj:`int`: 
            The threshold size.

    Returns
    -------
        :obj:`filter`
    """
    if (min_size is not None) & (max_size is not None):
        f = lambda x: (len(x) >= min_size) & (len(x) <= max_size)
    elif (min_size is not None) & (max_size is None):
        f = lambda x: len(x) >= min_size
    elif (min_size is None) & (max_size is not None):
        f = lambda x: len(x) <= max_size
    else:
        f = lambda x: x is not None
    return filter(f, sequences)

def find_links_fast(communities_grouped, max_parent_combinations, offsets,
        min_forward_containment=0, min_backward_containment=0, 
        reverse_mappings=None):
    """find_links_fast
    Under development.
    """
    
    n_largest_combos = 0
    edge_list = []

    for i, communities in enumerate(communities_grouped):
        if i == 0:
            past_community_lengths = np.array([len(c) for c in communities])
            if reverse_mappings is None:
                past_reverse_mapping = reverse_index(communities)
            continue
        for community_id, community in zip(range(*offsets[i]), communities):
            matches = np.array(
                    list((flatten(itemgetter(*community)(past_reverse_mapping))
                        )))
            if len(matches) == 0:
                continue
            match_ids, match_counts = np.unique(matches, return_counts=True)
            match_lengths = past_community_lengths[match_ids]
            community_lengths = np.array([len(community)] * match_ids.shape[0])

            backward_containment = match_counts / community_lengths
            forward_containment = match_counts / match_lengths
            # check for perfect links
            mask_perfect = (backward_containment == 1) & (forward_containment ==1)
            if np.sum(mask_perfect) > 0:
                for match in match_ids[mask_perfect]:
                    edge_list.append((community_id, match, 1, 1))

            mask = ((backward_containment >= min_backward_containment)
                    & (forward_containment >= min_forward_containment))
            if np.sum(mask) > 0:
            # check for perfect links
            
                if max_parent_combinations == 1:
                    jaccard = (match_counts / 
                            (match_lengths + community_lengths - match_counts))
                    for match, j in zip(match_ids[mask], jaccard[mask]):
                        edge_list.append((community_id, match, j, j))
                else:
                    match_ids = match_ids[mask]
                    combos = recursive_combinations(match_ids, max_parent_combinations)
                    match_combos = []
                    for combo in combos:
                        combo_communities = itemgetter(*combo)(communities[i-1])
                        if type(combo_communities) != tuple:
                            combo_communities = (combo_communities,)
                            print(combo_communities)
                        match_combos.append(set(flatten(combo_communities)))
                    # match_combos = [set(flatten(itemgetter(*combo)(communities[i-1]))) 
                    #       for combo in combos]
                    n_match_combos = len(match_combos)
                    if n_match_combos > n_largest_combos:
                        n_largest_combos = n_match_combos
                        logging.info(('Finding links from {n_match_combos} possible '
                                    'antecedant combinations.'))
                    community = set(community)

                    jaccard = [jaccard_similarity_set(community, set(match)) 
                           for match in match_combos]
                    j_max = np.max(jaccard)
                    match_combo_ids = np.nonzero(list(map(lambda x: x == j_max)))[0]
                    for mci, j in zip(itemgetter(*match_combo_ids)(combos), itemgetter(*match_combo_ids)(jaccard)):
                        for match in mci:
                            j = jaccard_set_similarity(community, set(communities[i-1][match]))
                            edge_list.append(community_id, match, j_max)
    return edge_list


def phylomemetic_graph(steps, communities, min_size=3, max_size=50, 
        parent_limit=2, workers='auto', chunksize='auto', method='fast', 
        min_backwards_containment=0, min_forward_containment=0):
    '''phylomemetic_graph

    Parameters
    ----------
        steps : :obj:`iter` of :obj:`int`
        communities : :obj:`iter` of :obj:`iter` of :obj:`int`
        min_size : :obj:`int`
        max_size : :obj:`int`
        parent_limit : :obj:`int`
        workers : :obj:`int`
        chunksize : :obj:`int`
        method : :obj:`str`
        min_backwards_containment : :obj:`float`
        min_forward_containment : :obj:`float`

    Returns
    -------
        g : :obj:`graph_too.Graph`
        group_link_strength : :obj:`graph_tool.EdgePropertyMap`
        single_link_strength : :obj:`graph_tool.EdgePropertyMap`
        vertex_steps : :obj:`graph_tool.VertexPropertyMap`
        element_vertex_map : :obj:`dict`
    '''
    if workers == 'auto':
        workers = cpu_count() - 1

    communities_filt = []
    communities_lengths = []
    element_community_mappings = []
    for sequences in communities:
        s_filt = list(filter_by_size(sequences, min_size, max_size))
        communities_filt.append(s_filt)
        communities_lengths.append(len(s_filt))
        element_community_mappings.append(reverse_index(s_filt))

    community_vertex_maps = []
    communities_offsets = []
    cumsum_lengths = np.cumsum(communities_lengths)

    for length, count in zip(communities_lengths, cumsum_lengths):
        start = count - length
        end = count
        communities_offsets.append((start, end))
        community_vertex_maps.append(
                {c: v for c, v in zip(range(length), range(start, end))})

    n_communities = np.sum(communities_lengths)

    phylomemetic_links = []

    for i, (cps, cfs) in enumerate(window(communities_filt, 2)):
        n_cf = len(cfs)
        logger.info(f'Processing {i+1} of {len(communities)-1} periods')
        if chunksize == 'auto':
            chunksize_i = int(np.ceil((1 / workers) * n_cf))
        else:
            chunksize_i = chunksize
        
        with Pool(workers) as pool:
            phylomemetic_links.append(
                    pool.map(
                        find_links,
                        zip(
                            cfs,
                            range(0, len(cfs)),
                            repeat(cps, n_cf),
                            repeat(communities_offsets[i], n_cf),
                            repeat(element_community_mappings[i], n_cf),
                            repeat(parent_limit, n_cf),
                            ),
                        chunksize=chunksize_i,
                        )
                    )
            pool.close()
            pool.join()

    g = Graph(directed=True)
    g.add_vertex(n_communities)

    group_link_strength = g.new_edge_property('float')
    single_link_strength = g.new_edge_property('float')

    phylomemetic_links = flatten(flatten(phylomemetic_links))
    g.add_edge_list(
            phylomemetic_links,
            eprops=[group_link_strength, single_link_strength]
            )

    element_vertex_map = reverse_index_communities(flatten(communities_filt))

    vertex_steps = g.new_vertex_property('int')
    for (start, end), step in zip(communities_offsets, steps):
        vertex_steps.a[start:end] = step

    return (g, group_link_strength, single_link_strength, 
            vertex_steps, element_vertex_map)

# class PhylomemeticGraph(Graph):
#     
#     def __init__(self, *args, **kwargs):
#         """PhylomemeticGraph
# 
#         A graph that links longtitudinal sets of communities using a
#         "phylomemetic" algorithm (Chavalarias and Cointet 2013).
#         """
#         super().__init__(*args, **kwargs)
# 
#     def from_communities(self, community_sets, labels=None, min_clique_size=None,
#             parent_limit=2, workers='auto', chunksize='auto'):
#         """from_communities
#         Builds a phylomemetic graph from a set of communities.
# 
#         Parameters
#     ----------
#             community_sets : :obj:`iter` of :obj:`iter` of :obj:`iter` Sets of
#                 elements grouped into communities.
#             labels : :obj:`iter` Labels that identify the period at which each
#                 set of communities occurs.
#             min_clique_size (int): The minimum community size to be consider for
#                 for a node in the network.
#             workers (int): The number of CPUs to use.
#             parent_limit (int): The maximum combination size for unions of 
#                 potential parental candidates.
#             chunksize (int or str): The number of communities for each CPU to
#                 process with in each process. Default is 'auto'.
# 
#         Returns
#     -------
#             self
#         """
# 
#         if workers == 'auto':
#             workers = cpu_count() - 1
# 
#         community_sets_filt = []
#         community_sets_lengths = []
#         element_community_mappings = []
#         for communities in community_sets:
#             cfsilt = self.filter_communities(
#                     communities, min_clique_size
#                     )
#             community_sets_filt.append(cfsilt)
#             community_sets_lengths.append(len(cfsilt))
#             element_community_mappings.append(
#                 reverse_index_communities(cfsilt)
#                 )
# 
#         community_vertex_maps = []
#         community_sets_pos = []
#         cumsum_lengths = np.cumsum(community_sets_lengths)
#         for length, count in zip(community_sets_lengths, cumsum_lengths):
#             start = count - length
#             end = count
#             community_sets_pos.append((start, end))
#             community_vertex_maps.append(
#                     {c: v for c, v in zip(range(length), range(start, end))}
#                         )
# 
#         n_communities = sum(community_sets_lengths)
# 
#         vocab_all = list(set(flatten([flatten(c) for c in community_sets])))
#         len_vocab = len(vocab_all)
#  
#         phylomemetic_links = []
#         for i, (cps, cfs) in enumerate(window(community_sets, 2)):
#             n_cf = len(cfs)
#             logger.info(f'Processing {i+1} of {len(community_sets)-1} periods')
#             if chunksize == 'auto':
#                 chunksize_i = int(np.ceil((1 / workers) * n_cf))
#             else:
#                 chunksize_i = chunksize
#             
#             with Pool(workers) as pool:
#                 phylomemetic_links.append(
#                         pool.map(
#                             find_links,
#                             zip(
#                                 cfs,
#                                 range(0, len(cfs)),
#                                 repeat(cps, n_cf),
#                                 repeat(community_sets_pos[i], n_cf),
#                                 repeat(element_community_mappings[i], n_cf),
#                                 repeat(parent_limit, n_cf),
#                                 ),
#                             chunksize=chunksize_i,
#                             )
#                         )
#                 pool.close()
#                 pool.join()
#         if len(list(self.vertices())) == 0:
#             self.add_vertex(n_communities)
# 
#         self.ep['group_link_strength'] = self.new_edge_property('float')
#         self.ep['single_link_strength'] = self.new_edge_property('float')
# 
#         phylomemetic_links = flatten(flatten(phylomemetic_links))
#         self.add_edge_list(
#                 phylomemetic_links,
#                 eprops=[self.ep['group_link_strength'], self.ep['single_link_strength']]
#                 )
# 
#         community_labels = self.new_vertex_property('int')
#         community_items = self.new_vertex_property('vector<int>')
#         for i, ((start, end), communities) in enumerate(
#                 zip(community_sets_pos, community_sets)):
#             vertices = range(start, end)
#             for vertex, c in zip(vertices, communities):
#                 community_items[vertex] = np.array(list(c))
#                 community_labels[vertex] = labels[i]
# 
#         self.vp['item'] = community_items
#         self.vp['label'] = community_labels
# 
#         return self
# 
#     def filter_communities(self, cliques_set, min_clique_size):
#         """filter_communities
#         Returns the communities in a series that are larger than the minimum
#         community size.
# 
#         Parameters
#     ----------
#             cliques_set : :obj:`iter` of :obj:`iter` Series of iterables
#                 containing the vertices that make up cliques.
#             min_clique_size (int): The threshold size.
# 
#         Returns
#     -------
#             csf : :obj:`iter` of :obj:`iter` Series of filtered cliques that
#                 have length greater than or equal to min_clique_size.
#         """
#         csf = [c for c in cliques_set if len(c) >= min_clique_size]
#         return csf
# 
