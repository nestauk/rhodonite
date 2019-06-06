import click
import pickle
import os
from gensim.corpora import Dictionary

from rhodonite.cooccurrence import CooccurrenceGraph

@click.command()
@click.option('--input', '-i', help='Input file of communities.')
@click.option('--output', '-o', help='Output file location. End with .gt suffix.')
@click.option(
        '--label_items', '-l', type=bool,
        help='Maximum number of possible parents to consider for each community.'
        )

@click.option(
        '--item_type', '-i', type=str,
        help='Maximum number of possible parents to consider for each community.'
        )

@click.option(
        '--distance_agg', '-i', type=function,
        help='Maximum number of possible parents to consider for each community.'
        )

def build(input, output, label_items, item_type, distance_agg):
        # input 1, input2
    """from_sequences
    Constructs a cooccurrence network from a series of sequences (an
    iterable of iterables), for example, a corpus of tokenized documents.
    Depending on the value of ``window_size``, either a sliding window is
    used to identify cooccurrences between neigbouring elements in the
    nested sequences, or cooccurrences are counted between all elements in
    each sequence.
    Args:
        sequences (:obj:`iter` of :obj:`iter` :obj:`iter`):
        dictionary (dict):
        label_items (bool):
        item_type (str):
        distance_agg (function):
    Returns:
        self
    """

    with open(input, 'rb') as f:
        sequences = pickle.load(f)

    dictionary = Dictionary(sequences)

    if callable(distance_agg) == True:
        distance_agg=distance_agg


    save_dir = os.path.join(*input.split(os.sep)[:-1])
    if os.path.isdir(save_dir):
        co_graphs = []
        for corpus_ids in corpora_ids:
            co = CooccurrenceGraph()
            co.from_sequences(
                sequences=corpus_ids,
                dictionary=dictionary,
                label_items=label_items,
                item_type=item_type,
                distance_agg=distance_agg
            )
            co.ep['association_strength'] = association_strength(co)
            co_graphs.append(co)

        co.save(output)

    else:
        click.echo('Output directory does not exist.')

if __name__ == '__main__':
    build()
