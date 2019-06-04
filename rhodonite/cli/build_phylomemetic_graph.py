import click
import pickle

from rhodonite.phylomemetic import PhylomemeticGraph

@click.command()
@click.option('--input', '-i', help='Input file of communities.')
@click.option('--output', '-o', help='Output file location. End with .gt suffix.')
@click.option(
        '--parent_limit', '-p',
        help='Maximum number of possible parents to consider for each community.'
        )
@click.option(
        '--min_clique_size', '-m', 
        help='Minimum number of elements for a community to be included.'
        )
@click.option(
        '--workers', '-w', default='auto',
        help=(
            'Number of processes to use. More workers means faster calculation, '
            'but higher memory overhead.'
            )
        )
@click.option(
        '--chunksize', '-c', default='auto',
        help='Maximum number of communities for each worker to work on.'
        )
def from_communities(input, output, min_clique_size, parent_limit, 
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
    
    save_dir = os.path(*input.split(os.sep)[:-1])
    if os.path.isdir(save_dir):
        pg = PhylomemeticGraph()
        pg.from_communities(
                communities, 
                min_clique_size=min_clique_size, 
                parent_limit=parent_limit,
                workers=workers,
                chunksize=chunksize,
                )
        pg.save(output)
    else:
        click.echo('Output directory does not exist.')

if __name__ == '__main__':
    from_communities()
