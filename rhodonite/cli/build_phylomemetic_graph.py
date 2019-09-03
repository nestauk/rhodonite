import click
import pickle
import os

from rhodonite.dynamic.phylomemetic import phylomemetic_graph

@click.command()
@click.option('--input', '-i', help='Input file of communities.')
@click.option('--output', '-o', help='Output file location. End with .gt suffix.')
@click.option(
        '--parent_limit', '-p', type=int,
        help='Maximum number of possible parents to consider for each community.'
        )
@click.option(
        '--min_clique_size', '-m', type=int, 
        help='Minimum number of elements for a community to be included.'
        )
@click.option(
        '--workers', '-w', default='auto', type=str,
        help=(
            'Number of processes to use. More workers means faster calculation, '
            'but higher memory overhead.'
            )
        )
@click.option(
        '--chunksize', '-c', default='auto', type=str,
        help='Maximum number of communities for each worker to work on.'
        )
def build(input, output, min_clique_size, parent_limit, 
        workers, chunksize):
    '''from_communities
    Creates and saves a phylomemetic graph from an input of temporal communities.

    Args:
        input (:obj:`str`): Path to input pickled dictionary of communities.
        output (:obj:`str`): Output directory for results (.gt format).
        min_clique_size (:obj:`int`): Minimum community size to consider.
        parent_limit (:obj:`int`): Maximum number of parents to consider.
        workers (:obj:`str` or :obj:`int`): Number of processes to use. Either 
            provide integer, or "auto". Default is "auto".
        chunksize (:obj:`str` or :obj:`int`): Number of communities for each 
            worker to process at a time. Either provide integer or "auto".
            Default is "auto".
    '''

    with open(input, 'rb') as f:
        communities = pickle.load(f)
    
    if chunksize.isnumeric():
        chunksize = int(chunksize)
    if workers.isnumeric():
        workers = int(workers)

    save_dir = os.path.join(*input.split(os.sep)[:-1])
    if os.path.isdir(save_dir):
        pg = phylomemetic_graph(
                community_sets=list(communities.values()),
                labels=list(communities.keys()),
                min_clique_size=min_clique_size, 
                parent_limit=parent_limit,
                workers=workers,
                chunksize=chunksize,
                )
        pg.save(output)
    else:
        click.echo('Output directory does not exist.')

if __name__ == '__main__':
    build()
