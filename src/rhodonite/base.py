from graph_tool.all import Graph


class CoocurrenceGraph(Graph):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    
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

    def build(self, sequences, dictionary=None):
        pass
