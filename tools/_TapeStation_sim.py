#!/usr/bin/env python3
__description__ =\
"""
Purpose: Perform double restriction digest on a given genome.
"""
__author__ = "Erick Samera; Michael Ke"
__version__ = "5.1.0"
__comments__ = "stable; multi-processing"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
import matplotlib.pyplot as plt
import pickle
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        'input',
        type=Path,
        help="path of input positions file (.pos)")

    args = parser.parse_args()
    return args
# --------------------------------------------------
def _get_fragment_size(start_pos: int, end_pos: int): return end_pos - start_pos
def _parse_fragments(_input_path: Path):
    """
    """
    
    fragment_sizes = []
    fragments_dict: dict = pickle.load(open(_input_path, mode='rb'))
    for chromosome in [i for i in fragments_dict if i !='metadata']:
        fragment_sizes += [_get_fragment_size(position[0], position[1]) for position in fragments_dict[chromosome]['fragment_positions']]
    
    x_vals = [i for i in range(max(fragment_sizes)+1)]
    y_vals = [0]*(max(fragment_sizes)+1)
    for fragment in fragment_sizes:
        y_vals[fragment] += 1
    
    return x_vals, y_vals
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    x_vals, y_vals = _parse_fragments(args.input)

    plt.plot(x_vals, y_vals)
    plt.xlim(10, 1000)
    plt.xscale("symlog")
    plt.show()

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()
