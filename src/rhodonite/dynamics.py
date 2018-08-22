from collections import defaultdict
from graph_tool.all import Graph, GraphView

from rhodonite.cliques import find_cliques_cfinder
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
        self.names = names
        self.normalisation_weight = normalisation_weight
        self.delta = delta
        self.min_clique_size = min_clique_size
        self.max_weight = max_weight

        self.new_vertex_property('vector<int>')

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
        
        phylomemetic_links = []

        for cliques_past, cliques_future in window(self.clique_sets, 2):
            self.add_vertex(len(cliques_past))
            self.add_edge_list(phylomemetic_links)

            overlapping_past_cliques = []

            cliques_past_unions = 

            for cf in cliques_future:
                jaccard_matches = []
                for cp in cliques_past:
                    if not set(cf).isdisjoint(cp):
                        j = jaccard(cf, cp)
                        if j >= self.delta:
                            jaccard_matches.append(j)
                        else:
                            jaccard_matches.append(0)

