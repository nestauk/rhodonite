import logging
import numpy as np
import pandas
import os

from collections import defaultdict
from graph_tool.all import Graph, GraphView
from itertools import repeat, combinations
from operator import itemgetter
from sklearn.preprocessing import MultiLabelBinarizer
from multiprocessing import Pool

from rhodonite.cliques import (find_cliques_cfinder, filter_subsets,
        clique_unions, reverse_index_cliques, load_cliques_cfinder)
from rhodonite.similarity import jaccard_similarity
from rhodonite.utilities import window, flatten, clear_graph, get_aggregate_vp
from rhodonite.tabular import vertices_to_dataframe

import time


logger = logging.getLogger(__name__)

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
            year_v = g.vp['times'][v]

            predecessors = list(v.in_neighbors())
            years = [g.vp['times'][p] for p in predecessors]
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
        pp_matrices (:obj:`iter` of :obj:`matrix`): Multi label binarized
            cliques from all time periods previous to cf.
        pos (:obj:`iter` of :obj:`iter`): The start and ending positions of all
            each set of cliques in the previous time periods, relative to all
            cliques.
        binarizer (:obj:`MultiLabelBinarizer`): A multi label binarizer fitted
            to the entire vocabulary across all cliques.
        delta_0 (int): The threshold Jaccard index value for potential parent
            cliques.
        parent_limit (int): The limit for the size of possible combinations
            of potential parent cliques that will be considered.
    
    Returns:
        links: (:obj: `list`): A list that contains the inter-temporal edges
            that have been found as well as the corresponding Jaccard indexes.
            Each element is of the format ((source, target) jaccard_index).
    """
    cf, cfi, cps, pp_matrices, pos, binarizer, delta_0, parent_limit = args
    
    pos_tmp = pos.copy()
    start_f, end_f = pos_tmp.pop()
    pos_tmp = list(reversed(pos_tmp))
    
    links = []
    for i, pp_matrix in enumerate(pp_matrices):
        start_p, end_p = pos_tmp[i]
        cf_matrix = binarizer.transform(
                [cf for _ in range(pp_matrix.shape[0])]
                )
        j = jaccard_similarity(cf_matrix, pp_matrix)
        if np.max(j) == 1:
            direct_parent = np.nonzero(j == 1)[0][0]
            links.append(
                   ((direct_parent + start_p, cfi + start_f), 1)
                    )
            return links

        # keep the first set of jaccard similarities in case no exact match
        if i == 0:
            j_immediate = j.copy()
    
    start_p, end_p = pos_tmp[0]
    cp_indexes = np.nonzero(j_immediate > delta_0)[0]
    
    if len(cp_indexes) > 0:
        cp_thresh = [cps[i] for i in cp_indexes]
        cp_union_indices = clique_unions(cp_indexes, parent_limit)

        cp_union_vertices = []
        for cui in cp_union_indices:
            cuv = set(flatten([cps[i] for i in cui]))
            cp_union_vertices.append(cuv)

        cp_matrix_thresh = binarizer.transform(cp_union_vertices) 
        cf_matrix = binarizer.transform(
                [cf for _ in range(cp_matrix_thresh.shape[0])]
                )
        j_thresh = jaccard_similarity(cf_matrix, cp_matrix_thresh)
        j_max = np.max(j_thresh)
        parent_clique_indices = np.nonzero(j_thresh == j_max)[0]
        if len(parent_clique_indices) > 0:
            parent_cliques = itemgetter(*parent_clique_indices)(cp_union_indices)
            if any(isinstance(i, tuple) for i in parent_cliques):
                parent_cliques = flatten(parent_cliques)
            j_parents = [np.max(j) for j in j_thresh if np.max(j) > 0]

            for pc, j in zip(parent_cliques, j_parents):
                links.append(((pc + start_p, cfi + start_f), j))
            return links

class PhylomemeticGraph(Graph):
    
    def __init__(self, *args, **kwargs):
        """PhylomemeticGraph

        A graph that links longtitudinal sets of communities using a
        "phylomemetic" algorithm (Chavalarias and Cointet 2013).
        """
        super().__init__(*args, **kwargs)

    def from_communities(self, community_sets, labels=None,
            min_clique_size=None, workers=4, delta_0=0.3, parent_limit=2,
            color=False):
        """from_communities
        Builds a phylomemetic graph from a set of communities.

        Args:
            community_sets (:obj:`iter` of :obj:`iter` of :obj:`iter`):
            labels (:obj:`iter`):
            min_clique_size (int):
            workers (int):
            delta_0 (float):
            parent_limit (int):

        Returns:
            self
        """

        community_sets_filt = []
        community_sets_lengths = []
        for communities in community_sets:
            communities_filt = self.filter_communities(
                    communities, min_clique_size
                    )
            community_sets_filt.append(communities_filt)
            community_sets_lengths.append(len(communities))

        community_sets_pos = []
        for length, count in zip(
                community_sets_lengths, np.cumsum(community_sets_lengths)):
            community_sets_pos.append((count - length, count))
        n_communities = sum(community_sets_lengths)

        vocab_all = list(set(flatten([flatten(c) for c in community_sets])))
        len_vocab = len(vocab_all)
        binarizer_all = MultiLabelBinarizer(
                classes=vocab_all,
                sparse_output=True)
        binarizer_all.fit(range(0, len_vocab))

        binarized_community_sets = [binarizer_all.transform(c)
                for c in community_sets]
        
        phylomemetic_links = []
        # find direct parents
        for i, (communities_p, communities_f) in enumerate(window(community_sets, 2)):
            n_cf = len(communities_f)
            cp_matrix = binarizer_all.transform(communities_p)

            possible_parent_matrices = list(
                    reversed(binarized_community_sets[:i+1])
                    )
            # include positions of the current
            positions = community_sets_pos[:i+2]

            with Pool(workers) as pool:
                phylomemetic_links.append(
                        pool.map(
                            find_links,
                            zip(
                                communities_f,
                                range(0, len(communities_f)),
                                repeat(communities_p, n_cf),
                                repeat(possible_parent_matrices, n_cf),
                                repeat(positions, n_cf),
                                repeat(binarizer_all, n_cf),
                                repeat(delta_0, n_cf),
                                repeat(parent_limit, n_cf),
                                )
                            )
                        )
                    
                pool.close()
                pool.join()
                                            
        if len(list(self.vertices())) == 0:
            self.add_vertex(n_communities)

        link_strengths = self.new_edge_property('float')
        for pl in phylomemetic_links:
            pl = flatten([p for p in pl if p is not None])
            pl_edges = [p[0] for p in pl]
            pl_weights = [p[1] for p in pl]
            self.add_edge_list(set(pl_edges))
            for e, w in zip(pl_edges, pl_weights):
                link_strengths[self.edge(e[0], e[1])] = w

        self.ep['link_strength'] = link_strengths
        
        if color:
            colors = [i / len(labels) for i in range(0, len(labels))]
            community_color = self.new_vertex_property('float')

        community_labels = self.new_vertex_property('int')
        community_items = self.new_vertex_property('vector<int>')

        for i, ((start, end), communities) in enumerate(
                zip(community_sets_pos, community_sets)):
            vertices = range(start, end)
            for vertex, c in zip(vertices, communities):
                community_items[vertex] = np.array(c)
                community_labels[vertex] = labels[i]
                if color:
                    community_color[vertex] = colors[i]

        self.vp['item'] = community_items
        self.vp['label'] = community_labels
        if color:
            self.vp['color'] = community_color

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
    

