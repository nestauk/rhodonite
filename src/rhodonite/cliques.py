import os
import shutil

from subprocess import call
from rhodonite.utilities import save_edgelist, check_and_create_dir


def find_cliques_cfinder(g, cfinder_path, output_dir=None, delete_outputs=True,
        weight=None, **opts):
    """find_cliques_cfinder
    Finds the cliques in a graph using the CFinder tool.

    Args:
        g (:obj:`Graph`):
        cfinder_path:
        output_dir:
        delete_outputs:
        weight:
        **opts: Dictionary of flag-value pairs representing CFinder input
            options. From the CFinder prompt, these are:
            -i  specify input file.                       (Mandatory)
            -l  specify licence file with full path.      (Optional)
            -o  specify output directory.                 (Optional)
            -w  specify lower link weight threshold.      (Optional)
            -W  specify upper link weight threshold.      (Optional)
            -d  specify number of digits when creating
                the name of the default output directory
                of the link weight thresholded input.     (Optional)
            -t  specify maximal time allowed for
                clique search per node.                   (Optional)
            -D  search with directed method.              (Optional)
            -U  search with un-directed method.           (Default)
                (Declare explicitly the input and the
                modules to be un-directed.)
            -I  search with intensity method and specify
                the lower link weight intensity threshold
                for the k-cliques.                        (Optional)
            -k  specify the k-clique size.                (Optional)
                (Advised to use it only when a
                link weight intensity threshold is set.)
 
    Retuns:
        cliques (:obj:`list` of :obj:`tuple`): A list of all of the cliques
            found by CFinder. Each clique is represented as a tuple of
            vertices.
    """
    opts = dict(**opts)

    if output_dir is None:
        output_dir = os.path.abspath(os.path.join(cfinder_path, os.pardir))
        output_dir = os.path.join(output_dir, 'output')
        opts['-o'] = output_dir
    else:
        opts['-o'] = output_dir	
    input_path = os.path.abspath(os.path.join(cfinder_path, os.pardir))
    input_path = os.path.join(input_path, 'graph_edges.txt')
    opts['-i'] = input_path

    check_and_create_dir(output_dir)
    if weight is not None:
        save_edgelist(g, input_path, weight=weight)
    else:
        save_edgelist(g, input_path)
    run_cfinder(cfinder_path, opts)
    cliques = load_cliques_cfinder(os.path.join(output_dir, 'cliques'))
    if delete_outputs:
        shutil.rmtree(output_dir)
    return cliques

def load_cliques_cfinder(file_path):
    """load_cliques
    Loads cliques from a CFinder output file into a list of tuples.

    Args:
        file_path (str): The path to the CFinder output file. This is normally
            in a directory of outputs and named "cliques".

    Returns:
        cliques (:obj:`list` of :obj:`tuple`): A list of all of the cliques
            found by CFinder. Each clique is represented as a tuple of
            vertices.
    """
    with open(file_path, 'r') as f:
        clique_data = f.read().splitlines()
    cliques = []
    for cd in clique_data:
        if len(cd) > 0:
            if cd[0].isdigit():
                clique = cd.split(' ')[1:-1]
                clique = tuple([int(i) for i in clique])
                cliques.append(clique)
    return cliques
            

def run_cfinder(cfinder_path, opts):
    opts_list = [cfinder_path]
    for flag, value in opts.items():
        opts_list.append(flag)
        opts_list.append(value)
    call(opts_list)

