from collections import defaultdict
from graph_tool.all import Graph, GraphView
import logging
import numpy as np
from operator import itemgetter
from sklearn.preprocessing import MultiLabelBinarizer
from scipy.sparse import csc_matrix

from rhodonite.cliques import find_cliques_cfinder, filter_subsets, clique_unions, reverse_index_cliques
from rhodonite.similarity import jaccard_similarity
from rhodonite.utilities import window, flatten

logger = logging.getLogger(__name__)

class PhylomemeticGraph(Graph):
    
    def __init__(self, graphs, weights, dictionary, times, 
            delta=0.5, parent_limit=3, min_clique_size=3, max_weight=None,
            min_weight=None, **kwargs):
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
        self.delta = delta
        self.parent_limit = parent_limit
        self.min_clique_size = min_clique_size
        
        self.max_weight = max_weight
        self.min_weight = min_weight

        self.n_cooccurrence_vertices = len(dictionary.keys())

    def prepare(self, cfinder_path, clique_output_dir):
        """prepare
        """

        self.clique_sets = []
        for graph, weight in zip(self.graphs, self.weights):
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
            cliques = find_cliques_cfinder(thresh_graph, cfinder_path)
            cliques = [c for c in cliques if len(c) >= self.min_clique_size]
            
            cliques_harmonised = []
            for c in cliques:
                ch = [graph.vp.vertex_token_idxs[v] for v in c]
                cliques_harmonised.append(ch)

            self.clique_sets.append(cliques_harmonised)

        return self

    def build(self, log_every=500):
        """build
        Creates the links in the phylomemetic network between cliques at time T
        and time T'.
        """
        if len(list(self.edges())) > 0:
            for e in self.edges():
                self.remove_edge(e)
        if len(list(self.vertices())) > 0:
            for v in list(reversed(sorted(self.vertices()))):
                self.remove_vertex(v)

        phylomemetic_links = []
        
        for t, (cliques_past, cliques_future) in enumerate(window(self.clique_sets, 2)):
            vertex_cp_mapping = reverse_index_cliques(cliques_past)
            num_cf = len(cliques_future)

            future_vertices_start = len(cliques_past)
#             cv = CountVectorizer(analyzer='char', dtype=bool,
#                     vocabulary=range(0, self.len_vocab))
#             cp_vertices_str = [(' ').join([str(c) for c in cuv])
#                     for cuv in cliques_past]
#             cv.fit(cp_vertices_str)
            vocab = list(set(flatten(cliques_past) + flatten(cliques_future)))
            cv = MultiLabelBinarizer(classes=vocab)
            # pre-fit for speed
            cv.fit(range(0, len(vocab)))

            for i, cf in enumerate(cliques_future):
                if i % log_every == 0:
                    logger.info("{} cliques processed out of for time {}".format(
                        i, len(cliques_future), t + 1))

                # find all groups of overlapping cliques
                cp_indexes = list(set([vertex_cp_mapping[v] for v in cf if v in vertex_cp_mapping]))
                # remove any groups that are subsets
                cp_indexes = list(filter_subsets(cp_indexes))
                # find maximal cliques
                cp_union_indices, cp_union_vertices = clique_unions(
                        cp_indexes, cliques_past, self.parent_limit)

                cp_union_vertices.append(cf)
            
#                     clique_vertices_str = [(' ').join([str(c) for c in cuv]) 
#                             for cuv in cp_union_vertices]
                clique_matrix = cv.transform(cp_union_vertices)
#                 clique_matrix = csc_matrix(clique_matrix)

                cf_vector = clique_matrix[-1]
                cp_matrix = clique_matrix[:-1]
                future_clique_vertex = i + future_vertices_start

                jaccard_similarities = jaccard_similarity(cf_vector, cp_matrix)
                try:
                    if np.max(jaccard_similarities) >= self.delta:
                        parent_clique_indices = flatten(np.argwhere(
                            jaccard_similarities == np.max(jaccard_similarities)))
                        parent_cliques = itemgetter(*parent_clique_indices)(cp_union_indices)
                        if any(isinstance(i, tuple) for i in parent_cliques):
                            parent_cliques = flatten(parent_cliques)
                        for pc in parent_cliques:
                            phylomemetic_links.append((pc, future_clique_vertex))
                except ValueError:
                    continue
        
        self.add_vertex(len(flatten(self.clique_sets)))
        self.add_edge_list(phylomemetic_links)

        clique_times = self.new_vertex_property('int')
        clique_terms = self.new_vertex_property('vector<int>')
        
        vertex = 0
        for i, cliques in enumerate(self.clique_sets):
            for c in cliques:
                clique_terms[vertex] = np.array(c)
                clique_times[vertex] = self.times[i]
                vertex += 1
        self.vp['terms'] = clique_terms
        self.vp['times'] = clique_times
        
        return self 
