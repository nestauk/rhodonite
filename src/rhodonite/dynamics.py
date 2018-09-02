import logging
import numpy as np
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
from rhodonite.utilities import window, flatten, clear_graph

import time


logger = logging.getLogger(__name__)

def find_links(args):
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
            print('aye')
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
    
    def __init__(self, graphs, weights, dictionary, times,
            max_weight=None, min_weight=None, **kwargs):
        """PhylomemeticGraph
        """
        super().__init__(**kwargs)
        self.graphs = graphs
        if type(weights) == str:
            self.weights = [g.edge_properties[weights] for g in self.graphs]
        else:
            self.weights = weights
        self.dictionary = dictionary
        self.times = times
        self.len_vocab = len(dictionary.token2id)
        
        self.max_weight = max_weight
        self.min_weight = min_weight

        self.n_cooccurrence_vertices = len(dictionary.keys())
        n_periods = len(times)
        self.colors = [i / n_periods for i in range(n_periods)]

    def prepare(self, cliques_dir, cfinder_path=None,):
        """prepare
        """
        self.clique_sets = []
        for graph, weight, time in zip(self.graphs, self.weights, self.times):
            if cfinder_path is not None:
                # create a subgraph based on the min and max weight thresholds
                if (self.min_weight is not None) & (self.max_weight is None):
                    thresh_graph = GraphView(graph,
                        efilt=lambda e: self.min_weight <= weight[e])
                elif (self.max_weight is not None) & (self.min_weight is None):
                    thresh_graph = GraphView(graph,
                        efilt=lambda e: weight[e] <= self.max_weight)
                elif (self.max_weight is not None) & (self.min_weight is not None):
                    thresh_graph = GraphView(graph,
                        efilt=lambda e: self.min_weight <= weight[e] <= self.max_weight)
                else:
                    thresh_graph = graph
                
                # only vertices with more than one edge can be part of a clique
                thresh_graph = GraphView(thresh_graph, vfilt=lambda v: v.out_degree() > 1)

                # find all cliques in the subgraph
                output_dir = os.path.join(cliques_dir, str(time))
                cliques = find_cliques_cfinder(thresh_graph, cfinder_path,
                        output_dir, delete_outputs=False)
            else:
                cliques = load_cliques_cfinder(os.path.join(cliques_dir, str(time), 'cliques'))

            cliques_harmonised = []
            for c in cliques:
                ch = [graph.vp.vertex_token_idxs[v] for v in c]
                cliques_harmonised.append(ch)

            self.clique_sets.append(cliques_harmonised)

        return self
    
    def filter_cliques(self, cliques_set, min_clique_size):
        """filter_cliques
        Returns the cliques in a series that are larger than the minimum clique
        size.

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

    def build(self, workers=4, min_clique_size=3, delta_0=0.5, parent_limit=3):
        """build
        Creates the links in the phylomemetic network between cliques at time T
        and time T'.

        Args:
            workers (int):
            min_clique_size (int):
            log_every
        """
        
        clique_sets = []
        clique_sets_lengths = []
        for clique_set in self.clique_sets:
            clique_set_filt = self.filter_cliques(clique_set, min_clique_size)
            clique_sets.append(clique_set_filt)
            clique_sets_lengths.append(len(clique_set_filt))
        
        clique_sets_pos =[]
        for length, count in zip(
                clique_sets_lengths, np.cumsum(clique_sets_lengths)):
            clique_sets_pos.append((count - length, count))
        
        total_n_cliques = sum(clique_sets_lengths)
        
        vocab_all = list(set(flatten([flatten(c) for c in clique_sets])))
        len_vocab = len(vocab_all)
        binarizer_all = MultiLabelBinarizer(
                classes=vocab_all,
                sparse_output=True)
        binarizer_all.fit(range(0, len_vocab))

        binarized_clique_sets = [binarizer_all.transform(c)
                for c in clique_sets]
        
        phylomemetic_links = []
        # find direct parents
        for i, (cliques_p, cliques_f) in enumerate(window(clique_sets, 2)):
            n_cf = len(cliques_f)
            cp_matrix = binarizer_all.transform(cliques_p)

            possible_parent_matrices = list(
                    reversed(binarized_clique_sets[:i+1])
                    )
            # include positions of the current
            positions = clique_sets_pos[:i+2]

            with Pool(workers) as pool:
                phylomemetic_links.append(
                        pool.map(
                            find_links,
                            zip(
                                cliques_f,
                                range(0, len(cliques_f)),
                                repeat(cliques_p, n_cf),
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
            self.add_vertex(total_n_cliques)

        jaccard_weights = self.new_edge_property('float')
        for pl in phylomemetic_links:
            pl = flatten([p for p in pl if p is not None])
            pl_edges = [p[0] for p in pl]
            pl_weights = [p[1] for p in pl]
            self.add_edge_list(set(pl_edges))
            for e, w in zip(pl_edges, pl_weights):
                jaccard_weights[self.edge(e[0], e[1])] = w

        self.ep['jaccard_weights'] = jaccard_weights

        clique_times = self.new_vertex_property('int')
        clique_terms = self.new_vertex_property('vector<int>')
        clique_color = self.new_vertex_property('float')
        clique_densities = self.new_vertex_property('float')
       
        for i, ((start, end), cliques) in enumerate(zip(clique_sets_pos, clique_sets)):
            vertices = range(start, end)
            for vertex, c in zip(vertices, cliques):
                clique_densities[vertex] = self.calculate_clique_density(
                        c, self.graphs[i])
                clique_terms[vertex] = np.array(c)
                clique_times[vertex] = self.times[i]
                clique_color[vertex] = self.colors[i]
        
        self.vp['density'] = clique_densities
        self.vp['terms'] = clique_terms
        self.vp['times'] = clique_times
        self.vp['color'] = clique_color
        
        return self

    def calculate_clique_density(self, clique_terms, g):
        """calculate_clique_density
        Calculate the density of a clique based on the number of occurrences
        and coocurrences of the terms within it. Based on Callon et al. 1991.

        Args:
            clique_terms (:obj:`iter` of int): A set of terms that comprise
                a single clique.
            g (:obj:`Graph`): The coocurrence graph from which the clique
                originated
        Returns:
            density (float): The density of the clique.
        """
        clique_terms = [g.tokenidx2vertex[i] for i in clique_terms]
        card = len(clique_terms)
        co = []
        o = []
        for i, j in combinations(clique_terms, 2):
            o_i = g.vp['occurrences'][i]
            o_j = g.vp['occurrences'][j]
            o.append(o_i * o_j)
            co.append(g.ep['cooccurrences'][(i, j)])
        density = np.sum(np.divide(np.square(co), o))
        return density

