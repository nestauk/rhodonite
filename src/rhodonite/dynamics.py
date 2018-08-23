from collections import defaultdict
from graph_tool.all import Graph, GraphView
from operator import itemgetter

from rhodonite.cliques import find_cliques_cfinder, clique_unions
from rhodonite.similarity import jaccard


class PhylomemticGraph(Graph):
    
    def __init__(self, graphs, weights, dictionary, delta=0.5,
            min_clique_size=3, max_weight=None, min_weight=None):
        """PhylomemeticGraph
        """
        self.graphs = graphs
        if type(normalisation_weight) == str:
            self.weights = [g.edge_properties[weights] for g in self.graphs]
        else:
            self.weights = weights
        self.dictionary = dictionary
        self.names = names
        self.normalisation_weight = normalisation_weight
        self.delta = delta
        self.min_clique_size = min_clique_size
        self.max_weight = max_weight

        self.new_vertex_property('vector<int>')

        self.n_cooccurrence_vertices = len(dictionary.keys())

    def prepare(self, cfinder_path, clique_output_dir):
        """prepare

        """

        self.clique_sets = []
        for graph, weight, name in zip(self.graphs, self.weights, self.names):
            # create a subgraph based on the min and max weight thresholds
            thresh_graph = GraphView(graph,
                    efilt=lambda e: min_weight <= weight <= max_weight)
            # find all cliques in the subgraph
            cliques = find_cliques_cfinder(thresh_graph, cfinder_path)
            cliques = [c for c in cliques if len(c) >= self.min_clique_size]
            self.clique_sets.append(clique_sets)

    def build(self):
        """build
        Creates the links in the phylomemetic network between cliques at time T
        and time T'.
        """
        
        for cliques_past, cliques_future in window(self.clique_sets, 2):
            phylomemetic_links = []

            self.add_vertex(len(cliques_past))
            self.add_edge_list(phylomemetic_links)
            
            future_vertices_start = len(list(self.vertices()))
            

            cp_union_indices, cp_union_vertices = clique_unions(cliques_past)
            cp_matrix = np.zeros(
                    (len(cp_union_vertices), self.n_cooccurences_vertices))
            for i, cuv in enumerate(cp_union_vertices):
                rows = [i] * len(cuv)
                cp_matrix[rows, cuv]

            for i, cf in enumerate(cliques_future):
                future_clique_vertex = i + future_vertices_start

                cf_vector = np.zeros(len(self.n_cooccurrence_vertices))
                cf_vector = np.put(cf_vector, cf, 1)
                jaccard_similarities = jaccard_similarity(cf_vector, cp_matrix)

                if np.max(jaccard_similarities) >= self.delta:
                    parent_clique_indices = flatten(np.argwhere(
                        jaccard_similarities == np.max(jaccard_similarities)))

                        parent_cliques = flatten(
                                itemgetter(*parent_clique_indices)(cliques_past))
                        for pc in parent_cliques:
                            phylomemetic_links.append((pc, future_clique_vertex))


