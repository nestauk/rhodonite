import logging
import numpy as np
import pandas
import os

from collections import defaultdict
from graph_tool.all import Graph, GraphView
from itertools import repeat, combinations
from operator import itemgetter
from sklearn.preprocessing import MultiLabelBinarizer
from multiprocessing import Pool, cpu_count

from rhodonite.utilities import (window, flatten, clear_graph, 
        get_aggregate_vp, reverse_mapping)
from rhodonite.cliques import (filter_subsets, clique_unions,
        reverse_index_cliques, load_cliques_cfinder)
from rhodonite.similarity import jaccard_similarity, jaccard_similarity_set
from rhodonite.tabular import vertices_to_dataframe

import time


def label_ages(g):
    """label_ages
    Get the ages of each vertex in a graph and put them into
    a property map. The age is defined as the number of steps between
    a vertex and its most distant antecedent.
    
    Args:
        g (:obj:`Graph`): A graph.
        
    Returns:
        age_vp (:obj:`PropertyMap(`): A property map containing
            the age of each vertex.
    """

    # get all vertices with age 0 in dictionary
    ages = {}
    for v in g.vertices():
        if v.in_degree() == 0:
            ages[v] = 0

    # find youngest in-neighbour of each node
    # if it is in the age dict then get the new age by adding the difference
    # else append back on to the list to try again
    vertices = list(g.vertices())
    for v in vertices:
        if v in ages:
            continue
        else:
            year_v = g.vp['label'][v]

            predecessors = list(v.in_neighbors())
            years = [g.vp['label'][p] for p in predecessors]
            min_i = np.argmin(years)
            min_neighbor = predecessors[min_i]
            year_neighbor = years[min_i]
        if min_neighbor in ages:
            year_parent = ages[min_neighbor]
            ages[v] = ages[min_neighbor] + (year_v - year_neighbor)
        else:
            vertices.append(v)

    age_vp = g.new_vertex_property('int')

    for v in g.vertices():
        age_vp[v] = ages[v]
    
    return age_vp

def community_density(community, g):
    """community_density
    Calculate the density of a clique based on the number of occurrences
    and coocurrences of the terms within it. Based on Callon et al. 1991.

    Args:
        community (:obj:`iter` of int): A set of terms that comprise
            a single clique.
        g (:obj:`Graph`): The coocurrence graph from which the clique
            originated
    Returns:
        density (float): The density of the clique.
    """
    card = len(community)
    co = []
    o = []
    for i, j in combinations(community, 2):
        o_i = g.vp['occurrence'][i]
        o_j = g.vp['occurrence'][j]
        o.append(o_i * o_j)
        co.append(g.ep['cooccurrence'][(i, j)])
    density = 1 / card * np.sum(np.divide(np.square(co), o))
    return density

def label_density(g, cooccurrence_graphs, norm=None):
    """label_density
    Creates a property map with the density of each vertex based on the items
    contained within the community it represents.
    Requires a cooccurrence graphs.

    Args:
        g (:obj:`PhylomemeticGraph`): A phylomemetic graph.
        cooccurrence_graphs (:obj:`iter` of :obj:`CooccurrenceGraph`): A list
            of cooccurrence graphs for each time period.
        norm (function): A normalisation function.

    Returns:
        community_densities (:obj:`PropertyMap`): A property map containing the
            densities of the phylomemetic graph vertices.
    """
    community_densities = g.new_vertex_property('float')

    g_df = vertices_to_dataframe(g)
    time_steps = sorted(g_df['label'].unique())
    label_groups = g_df.groupby('label')

    for (_, group), co_graph in zip(label_groups, cooccurrence_graphs):
        densities = [community_density(c, co_graph) for c in group['item']]
        if norm is not None:
            densities = np.array(densities) / norm(densities)
        for v, d in zip(group['vertex'], densities):
            community_densities[v] = d

    return community_densities

def label_emergence(g):
    """label_emergence
    Creates a property map that identifies whether a vertex is categorised as
    ephemeral, emerging, steady, or declining.

    The categories are represented as integers with the following mapping:
        {'ephemeral': 0,
         'emerging': 1,
         'steady': 2,
         'declining': 3}
         
    Args:
        g (:obj:`Graph`): A graph.

    Returns:
        emergence (:obj:`PropertyMap`): An integer property map that represents
            the emergence category of each vertex.
    """
    emergence = g.new_vertex_property('int')
    for v in g.vertices():
        if (v.in_degree() == 0) & (v.out_degree() == 0):
            # ephemeral
            emergence[v] = 0
        elif (v.in_degree() == 0) & (v.out_degree() > 0):
            # emerging
            emergence[v] = 1
        elif (v.in_degree() > 0) & (v.out_degree() == 0):
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
    
    Args:
        g (:obj:`Graph`): A graph.

    Returns:
        branching (:obj:`PropertyMap`): A boolean property map that is True
            where a vertex's out degree is greater than 1.
        merging (:obj:`PropertyMap`): A boolean property map that is True where
            a vertex's in degree is greater than 1.
    """
    branching = g.new_vertex_property('bool')
    merging = g.new_vertex_property('bool')
    for v in g.vertices():
        if (v.in_degree() < 2) & (v.out_degree() >= 2):
            # branching
            branching[v] = True
            merging[v] = False
        elif (v.in_degree() >= 2) & (v.out_degree() < 2):
            # merging
            branching[v] = False
            merging[v] = True
        elif (v.in_degree() >= 2) & (v.out_degree() >= 2):
            # branching and merging
            branching[v] = True
            merging[v] = True
    return branching, merging

def label_cross_pollination(g, merging_prop, agg=np.mean):
    """label_cross_pollination
    Cross-pollination is defined as the average pairwise Jaccard similarity
    between the parents of a merging node.

    Args:
        g (:obj:`Graph`): A phylomemetic graph.
        merging_prop (:obj:`PropertyMap`): A verted property map that is True
            where a vertex is a merging node.
        agg(:obj:`function`): An aggregation function. Typically mean or
            median.
    """
    cross_poll_prop = g.new_vertex_property('float')
    g_merging = GraphView(g, vfilt=lambda v: merging_prop[v] ==  1)
    parents = []
    for v in g_merging.vertices():
        parents = [g.vp['item'][p] for p in g.vertex(v).in_neighbors()]
        jaccard = agg(
            [jaccard_similarity(list(c[0]), list(c[1]))
            for c in combinations(parents, 2)]
        )
        cross_poll_prop[v] = jaccard
    return cross_poll_prop

def label_diversification(g, branching_prop, agg=np.mean):
    """label_diversification
    Diversification is defined as the average pairwise Jaccard similarity
    between the children of a branching node.

    Args:
        g (:obj:`Graph`): A phylomemetic graph.
        branching_prop (:obj:`PropertyMap`): A verted property map that is True
            where a vertex is a branching node.
        agg(:obj:`function`): An aggregation function. Typically mean or
            median.
    """
    diversification_prop = g.new_vertex_property('float')
    g_branching = GraphView(g, vfilt=lambda v: branching_prop[v] ==  1)
    parents = []
    for v in g_branching.vertices():
        children = [g.vp['item'][c] for c in g.vertex(v).out_neighbors()]
        jaccard = agg(
            [jaccard_similarity(list(c[0]), list(c[1]))
            for c in combinations(children, 2)]
        )
        diversification_prop[v] = jaccard
    return diversification_prop

def find_links(args):
    """find_links
    Finds the inter-temporal links between a clique at time period and cliques
    occurring at previous time periods.

    Args:
        cf (:obj:`iter` of int): The clique for which parents are being
            identified.
        cfi (int): The clique index relative to all cliques in the phylomemy.
        cps (:obj:`iter` of :obj:`iter` of int): The cliques from the
            immediately previous time period to cf.
        pos (:obj:`tuple`): The start and end vertices in the previous time
            period.
        cpm (:obj:`scipy.sparse`): Binarized sparse matrix of communties from
            previous time period.
        mapping (:obj:`dict`): Maps community elements to communities.
        binarizer (:obj:`MultiLabelBinarizer`): A multi label binarizer fitted
            to the entire vocabulary across all cliques.
        parent_limit (int): The limit for the size of possible combinations
            of potential parent cliques that will be considered.
        threshold (float): The minimum percentile for the Jaccard similartiy
            between parents and children, below which parent candidates will
            be discarded.
    
    Returns:
        links: (:obj: `list`): A list that contains the inter-temporal edges
            that have been found as well as the corresponding Jaccard indexes.
            Each element is of the format:
            ((source, target), group_jaccard_index, single_jaccard_index).
            group_jaccard_index is the Jaccard index of any union of parents
            that has been identified as the most likely antecedants of the and
            single_jaccard_index is the individual Jaccard similarity between
            each potential parent and the child.
    """
    cf, cfi, cps, pos, cpm, mapping, binarizer, parent_limit, threshold = args

    start_p, end_p = pos
    start_f = end_p
    links = []
    possible_parents = []
    cp_indexes = []
    for element in cf:
        element_matches = mapping[element]
        if len(element_matches) == 0:
            continue
        else:
            cp_matches = [em for em in element_matches]
            match_matrix = cpm[cp_matches, :]
            cf_matrix = binarizer.transform(
                    [cf for _ in range(match_matrix.shape[0])]
                    )
            j = jaccard_similarity(cf_matrix, match_matrix)

            if np.max(j) == 1:
                direct_parent = np.nonzero(j == 1)[0][0]
                parent_cp = cp_matches[direct_parent]
                links.append(
                       ((cfi + start_f, parent_cp + start_p), 1, 1)
                        )
                return links

            thresh = np.percentile(j, threshold)
            high_j = np.nonzero(j >= thresh)[0]
            close_matches = [cp_matches[i] for i in high_j]
            possible_parents.append(close_matches)
    
    cp_indexes = sorted(set(flatten(possible_parents)))
     
    if len(cp_indexes) > 0:
        cp_thresh = [cps[i] for i in cp_indexes]
        cp_union_indices = clique_unions(cp_indexes, parent_limit)

        cp_union_vertices = []
        for cui in cp_union_indices:
            cuv = set(flatten([cps[i] for i in cui]))
            cp_union_vertices.append(cuv)

        cp_union_matrix = binarizer.transform(cp_union_vertices) 
        cf_matrix = binarizer.transform(
                [cf for _ in range(cp_union_matrix.shape[0])]
                )
        j = jaccard_similarity(cf_matrix, cp_union_matrix)
        cp_union_j_mapping = {k[0]: v[0, 0] for k, v in zip(cp_union_indices, j) if len(k) == 1}
        j_max = np.max(j)
        parent_clique_indices = np.nonzero(j == j_max)[0]
        if len(parent_clique_indices) > 0:
            parent_cliques = itemgetter(*parent_clique_indices)(cp_union_indices)
            j_parents = [np.max(j_p) for j_p in j if np.max(j_p) > 0]
            if any(isinstance(i, tuple) for i in parent_cliques):
                parent_cliques = set(flatten(parent_cliques))

            for pc in parent_cliques:
                link = (
                        (cfi + start_f, pc + start_p),
                        j_max,
                        cp_union_j_mapping[pc]
                    )
                if len(link) == 2:
                    print(link)
                links.append(link)
            return links


class PhylomemeticGraph(Graph):
    
    def __init__(self, *args, **kwargs):
        """PhylomemeticGraph

        A graph that links longtitudinal sets of communities using a
        "phylomemetic" algorithm (Chavalarias and Cointet 2013).
        """
        super().__init__(*args, **kwargs)

    def from_communities(self, community_sets, labels=None,
            min_clique_size=None, parent_limit=2, match_threshold=95,
            workers='auto', chunksize=0.2):
        """from_communities
        Builds a phylomemetic graph from a set of communities.

        Args:
            community_sets (:obj:`iter` of :obj:`iter` of :obj:`iter`): Sets of
                elements grouped into communities.
            labels (:obj:`iter`): Labels that identify the period at which each
                set of communities occurs.
            min_clique_size (int): The minimum community size to be consider for
                for a node in the network.
            workers (int): The number of CPUs to use.
            parent_limit (int): The maximum combination size for unions of 
                potential parental candidates.
            chunksize (int or float): The number of communities for each CPU to
                process with in each process.

        Returns:
            self
        """

        if workers == 'auto':
            workers = cpu_count() - 1

        community_sets_filt = []
        community_sets_lengths = []
        element_community_mappings = []
        for communities in community_sets:
            cfsilt = self.filter_communities(
                    communities, min_clique_size
                    )
            community_sets_filt.append(cfsilt)
            community_sets_lengths.append(len(cfsilt))
            element_community_mappings.append(
                reverse_mapping(cfsilt)
                )

        community_vertex_maps = []
        community_sets_pos = []
        cumsum_lengths = np.cumsum(community_sets_lengths)
        for length, count in zip(community_sets_lengths, cumsum_lengths):
            start = count - length
            end = count
            community_sets_pos.append((start, end))
            community_vertex_maps.append(
                    {c: v for c, v in zip(range(length), range(start, end))}
                        )

        n_communities = sum(community_sets_lengths)

        vocab_all = list(set(flatten([flatten(c) for c in community_sets])))
        len_vocab = len(vocab_all)
        binarizer_all = MultiLabelBinarizer(
                classes=vocab_all,
                sparse_output=True)
        binarizer_all.fit(range(0, len_vocab))

        community_matrices = [binarizer_all.transform(c)
                for c in community_sets]
 
        phylomemetic_links = []
        # find direct parents
        for i, (cps, cfs) in enumerate(window(community_sets, 2)):
            if type(chunksize) == float:
                n_processes = len(cfs)
                chunksize_i = int(np.ceil(n_processes * chunksize))
            else:
                chunksize_i = chunksize

            n_cf = len(cfs)
            cp_matrix = binarizer_all.transform(cps)
            
            with Pool(workers) as pool:
                phylomemetic_links.append(
                        pool.map(
                            find_links,
                            zip(
                                cfs,
                                range(0, len(cfs)),
                                repeat(cps, n_cf),
                                repeat(community_sets_pos[i], n_cf),
                                repeat(community_matrices[i], n_cf),
                                repeat(element_community_mappings[i], n_cf),
                                repeat(binarizer_all, n_cf),
#                                 repeat(delta_0, n_cf),
                                repeat(parent_limit, n_cf),
                                repeat(match_threshold, n_cf),
                                ),
                            chunksize=chunksize_i,
                            )
                        )
                pool.close()
                pool.join()
        if len(list(self.vertices())) == 0:
            self.add_vertex(n_communities)

        group_link_strengths = self.new_edge_property('float')
        single_link_strengths = self.new_edge_property('float')
        for pl in phylomemetic_links:
            pl = flatten([p for p in pl if p is not None])
            pl_edges = [p[0] for p in pl]
            pl_group_weights = [p[1] for p in pl]
            pl_single_weights = [p[2] for p in pl]
            self.add_edge_list(set(pl_edges))
            for e, gw, sw in zip(pl_edges, pl_group_weights, pl_single_weights):
                group_link_strengths[self.edge(e[0], e[1])] = gw
                single_link_strengths[self.edge(e[0], e[1])] = sw

        self.ep['group_link_strength'] = group_link_strengths
        self.ep['single_link_strength'] = single_link_strengths

        community_labels = self.new_vertex_property('int')
        community_items = self.new_vertex_property('vector<int>')

        for i, ((start, end), communities) in enumerate(
                zip(community_sets_pos, community_sets)):
            vertices = range(start, end)
            for vertex, c in zip(vertices, communities):
                community_items[vertex] = np.array(c)
                community_labels[vertex] = labels[i]

        self.vp['item'] = community_items
        self.vp['label'] = community_labels

        return self

    def filter_communities(self, cliques_set, min_clique_size):
        """filter_communities
        Returns the communities in a series that are larger than the minimum
        community size.

        Args:
            cliques_set (:obj:`iter` of :obj:`iter`): Series of iterables
                containing the vertices that make up cliques.
            min_clique_size (int): The threshold size.

        Returns:
            csf (:obj:`iter` of :obj:`iter`): Series of filtered cliques that
                have length greater than or equal to min_clique_size.
        """
        csf = [c for c in cliques_set if len(c) >= min_clique_size]
        return csf

