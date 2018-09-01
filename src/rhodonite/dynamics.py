import logging
import numpy as np
import os

from collections import defaultdict
from graph_tool.all import Graph, GraphView
from itertools import repeat
from operator import itemgetter
from sklearn.preprocessing import MultiLabelBinarizer
from multiprocessing import Pool

from rhodonite.cliques import (find_cliques_cfinder, filter_subsets,
        clique_unions, reverse_index_cliques, load_cliques_cfinder)
from rhodonite.similarity import jaccard_similarity
from rhodonite.utilities import window, flatten, clear_graph

import time

logger = logging.getLogger(__name__)

# def find_exact_parents(args):
#     exact_phylomemetic_links = []
#     child_indices = []
#     binarize all clique sets
#     for clique_set in clique_sets[1:]:
#         for i, clique in enumerate(clique_set):
#             for past_cliques in reversed(clique_sets[:i]):
#                 j = jaccard_similarity
#                 if argmax(j) == 1:
#                     exact_phylomemetic_links.append()
#                     child_indices.append()
#                     break
#     return exact_phylomemetic_links, child_indices

def find_links(args):
    cf, cfi, cps, cp_matrix, ppcs, pos, binarizer, delta_0, parent_limit = args
    
    pos_tmp = pos.copy()
    start_f, end_f = pos_tmp.pop()
    pos_tmp = list(reversed(pos_tmp))
    
    links = []
    for i, ppc in enumerate(ppcs):
        print(ppc)
        pp_matrix = binarizer.transform(ppc)
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
        
        parent_clique_indices = np.nonzero(j_thresh)[0]
        parent_cliques = itemgetter(*parent_clique_indices)(cp_union_indices)
        if any(isinstance(i, tuple) for i in parent_cliques):
            parent_cliques = flatten(parent_cliques)
        j_parents = [np.max(j) for j in j_thresh if np.max(j) > 0]

        for pc, j in zip(parent_cliques, j_parents):
            links.append(((pc + start_p, cfi + start_f), j))
        return links

# def find_links(args):
#     """find_links
#     Finds phylomemetic links between a set of cliques at a time T and a
#     set of cliques in the future at time T'.
# 
#     Args:
#         clique_f (:obj:`tuple`): A clique from time T'.
#         clique_f_i
#         cliques_p (:obj:`list` of :obj:`tuple`): The set of cliques
#             from time T.
#         vertex_cliques_p_mapping (:obj:`dict`): A reverse index mapping
#             of clique values to the list of cliques that they appear
#             within. e.g. {42: [1, 4, 5, 19]}
# 
#     Returns:
#         
#     """
#     cf, cfi, cpc, cp_matrix, cps, transformer, parent_limit, delta_0, delta_1 = args
# 
#     future_clique_vertex = cfi
#     
#     cf_matrix = transformer.transform(
#             [cf for _ in range(len(cps))]
#             )
#     jaccard_similarities = jaccard_similarity(cf_matrix, cp_matrix)
#     cp_indexes = np.nonzero(
#         jaccard_similarities > delta_0)[0]
# 
#     if len(cp_indexes) > 0:
#         cp_thresh = [cps[i] for i in cp_indexes]
#         map_back = {i_c: i_p for i_c, i_p in zip(range(len(cp_indexes)), cp_indexes)}
#         cp_union_indices = clique_unions(cp_indexes, parent_limit)
#         cp_union_vertices = []
#         for cui in cp_union_indices:
#             cuv = set(flatten([cps[i] for i in cui]))
#             cp_union_vertices.append(cuv)
#         cp_matrix_thresh = transformer.transform(cp_union_vertices)
#         cf_matrix = transformer.transform(
#                 [cf for _ in range(len(cp_union_vertices))]
#             )
# 
#         jaccard_similarities_1 = jaccard_similarity(cf_matrix, cp_matrix_thresh)
#         links = []
#         try:
#             if np.max(jaccard_similarities_1) >= delta_1:
#                 parent_clique_indices = flatten(np.argwhere(
#                     jaccard_similarities_1 == np.max(jaccard_similarities_1)))
#                 parent_cliques = itemgetter(*parent_clique_indices)(cp_union_indices)
#                 if any(isinstance(i, tuple) for i in parent_cliques):
#                     parent_cliques = flatten(parent_cliques)
#                 for pc in parent_cliques:
#                     links.append((pc + cpc, future_clique_vertex))
#                 return links
#         except ValueError:
#             pass

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
            vocab = list(set(flatten(cliques_p) + flatten(cliques_f)))
#             binarizer = MultiLabelBinarizer(classes=vocab, sparse_output=True)
            # pre-fit to full vocab for speed
#            binarizer.fit(range(0, len(vocab)))
            cp_matrix = binarizer_all.transform(cliques_p)

            possible_parent_matrices = list(
                    reversed(binarized_clique_sets[:i+1])
                    )
            possible_parent_clique_sets = list(reversed(
                clique_sets[:i+1]))
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
                                repeat(cp_matrix, n_cf),
                                repeat(possible_parent_clique_sets, n_cf),
                                repeat(positions, n_cf),
                                repeat(binarizer_all, n_cf),
                                repeat(delta_0, n_cf),
                                repeat(parent_limit, n_cf),
#                                     repeat(delta_1, 2),
                                )
                            )
                        )
                    
                pool.close()
                pool.join()
                                        
            
#         phylomemetic_links = []
#         cliques_past_count = 0 
#         for t, (cliques_past, cliques_future) in enumerate(window(self.clique_sets, 2)):
#             cliques_past = self.filter_cliques(cliques_past, min_clique_size)
#             cliques_future = self.filter_cliques(cliques_future, min_clique_size)
# #             vertex_cp_mapping = reverse_index_cliques(cliques_past)
# 
#             n_cf = len(cliques_future)
# 
#             vocab = list(set(flatten(cliques_past) + flatten(cliques_future)))
#             binarizer = MultiLabelBinarizer(classes=vocab, sparse_output=True)
#             # pre-fit to full vocab for speed
#             # means will only need to transform on each iteration
#             binarizer.fit(range(0, len(vocab)))
#             cf_matrix = binarizer.transform(cliques_past)
# 
#             cf_vertex_start = len(cliques_past) + cliques_past_count
#             cf_vertex_end = cf_vertex_start + n_cf
#             chunksize = int(n_cf / 10 * workers)
#             indices = range(cf_vertex_start, cf_vertex_end)
# 
#             with Pool(workers) as pool:
#                 phylomemetic_links.append(pool.map(
#                         find_links,
#                             zip(cliques_future, 
#                                 indices,
#                                 repeat(cliques_past_count, n_cf),
#                                 repeat(cf_matrix, n_cf),
#                                 repeat(cliques_past, n_cf),
#                                 repeat(binarizer, n_cf),
#                                 repeat(self.parent_limit, n_cf),
#                                 repeat(self.delta_0, n_cf),
#                                 repeat(self.delta_1, n_cf),
#                                 ),
#                             chunksize=chunksize
#                             ))
#                 pool.close()
#                 pool.join()
#             cliques_past_count += len(cliques_past)
# 
#         filtered_cliques = []
#         for clique_set in self.clique_sets:
#             filtered_cliques.append(
#                     self.filter_cliques(clique_set, min_clique_size))
# 
#         total_cliques = sum([len(c) for c in cliques_set])
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
       
        for i, ((start, end), cliques) in enumerate(zip(clique_sets_pos, clique_sets)):
            vertices = range(start, end)
            for vertex, c in zip(vertices, cliques):
                clique_terms[vertex] = np.array(c)
                clique_times[vertex] = self.times[i]
                clique_color[vertex] = self.colors[i]

#         vertex = 0
#         for i, cliques in enumerate(clique_sets):
#             for c in cliques:
#                 clique_terms[vertex] = np.array(c)
#                 clique_times[vertex] = self.times[i]
#                 clique_color[vertex] = self.colors[i]
#                 vertex += 1
        self.vp['terms'] = clique_terms
        self.vp['times'] = clique_times
        self.vp['color'] = clique_color
        
        return self
