from graph_tool.all import Graph


class CoocurrenceGraph(Graph):

    def __init__(self, sequences, dictionary=None, **kwargs):
        super().__init__(**kwargs)
        self.sequences = sequences
        self.dictionary = dictionary
    
    def sequence2vertices(self, sequence):
        """sequence2vertices
        Converts a sequence of labels to their vertex IDs.

        Args:
            sequence (:obj:`iter`): A sequence of labels or items to be
                converted.

        Returns:
            sequence_ids (:obj:`iter`): The original sequence of elements
                converted to their corresponding vertex IDs.
        """
        sequence_ids = [self.label2vertex[e] for e in sequence]
        return sequence_ids

    def calculate_occurrences(self, sequences_ids):
            """calculate_occurrences
            Calculates the number of times each element in the sequences occurs
            and puts them in to an array ordered by the elements' vertex IDs.

            Args:
                sequences (:obj:`iter` of :obj:`iter` of :obj:`int`): Sequences
                    converted to the ID representations.

            Returns:
                o (:obj:`np.array`): Array of occurrence counts sorted by
                    vertex order.
            """
        
        occurrence_counts = Counter(flatten(sequences))
        o = np.zeros(len(occurrence_counts))
        for i, count in occurrence_counts.items():
            o[i] = count
        self.occurrences = o

    def get_occurrence(self, label):
        """get_occurrences
        Returns the number of occurrences for a label from the occurrences
        vector.

        Args:
            label: A label from the original set of sequences.

        Returns:
            The number of times the label occurred within the sequences.
        """

        return self.occurences[self.label2vertex[label]]
 
    def calculate_coocurrences(self):
        """coocurrence
        Creates a coocurrence matrix from a counter of edge coocurrences. This
        is equivalent to an adjacency matrix, weighted by coocurrence.

        Args:
            coocurrences (:obj:`Counter`): Counter object mapping edges to
                the number of coocurrences between the source and target.
            n_vertices (int): The total number of vertices in the network.

        Returns:
            c (:obj:`np.array`): Coocurrence matrix.
        """
        
        c = np.zeros((self.n_vertices, self.n_vertices))
        for e, cs in coocurrences:
            s = e[0]
            t = e[1]
            c[s][t] = c[t][s] = cs
        return c

    def get_coocurrence(self, label_0, label_1):
        """get_edge_coocurrence
        Get the number of coocurrences between two labels.

        Args:
            label_0: A label from the sequence.
            label_1: A label from the sequences

        Returns:
            co: The number of times the two labels coocur.
        """
        i_0 = self.label2vertex[label_0]
        i_1 = self.label2_vertex[label_1]
        co = self.coocurrences[i_0][i_1]
        return co

    def association_strength(self):
        """assocation_strength_adj
        Calculates the symmetric association strength matrix from the
        edge coocurrence matrix and vertex occurrence vector.

        Args:
            n_cooccurrences (int): Total number of coocurrences.
            coocurs (:obj:`np.matrix`): Coocurrence matrix of the network.
            occurs (:obj:`np.array`): Array containing the total number of 
                occurrences for each vertex in the network.

        Returns:
            a_s (:obj:`np.matrix`): The association strength matrix for the
                coocurrence matrix.
        """
        

        a_s = np.divide((2 * self.n_coocurrences * self.coocurences),
                np.multiply(self.occurrences, self.occurrences.transpose()))
        return a_s

    def get_association_strength(self, label_0, label_1):
        """get_association_strength
        Returns the assocation strength between two labels from the assocation
        strength matrix.
 
        Args:
            label_0: A label from the sequence.
            label_1: A label from the sequences

        Returns:
            a_s: The association strength between the two labels.
        """
        i_0 = self.label2vertex[label_0]
        i_1 = self.label2_vertex[label_1]
        a_s = self.association_strengths[i_0][i_1]
        return a_s
        
