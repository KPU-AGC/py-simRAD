#!/usr/bin/env python3
__description__ =\
"""
Purpose: After performing double restriction enzyme digest, visualize the histogram of fragments.
"""
__author__ = "Erick Samera"
__version__ = "2.0.0"
__comments__ = "stable concept"
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
        'input_path',
        type=Path,
        help="path of input fragment length data (.pickle)")
    parser.add_argument(
        '-o',
        '--out',
        dest='output_path',
        metavar='FILE',
        type=Path,
        required=True,
        help="path of output file (.png/.jpg/.svg)")

    group_size_filter = parser.add_argument_group('fragment size filtering and binning options')
    group_size_filter.add_argument(
        '--filter',
        dest='filter',
        action='store_true',
        help="if --filter flag is active, filter the fragments to the binning argument (default=False)")
    group_size_filter.add_argument(
        '-m',
        '--min',
        dest='size_min',
        metavar='INT',
        type=int,
        default=100,
        help="min size of fragments to include in histogram (inclusive) (default=100)")
    group_size_filter.add_argument(
        '-M',
        '--max',
        dest='size_max',
        metavar='INT',
        type=int,
        default=600,
        help="max size of fragments to include in histogram (inclusive) (default=600)")
    group_bin_options = parser.add_mutually_exclusive_group(required=True)
    group_bin_options.add_argument(
        '-b',
        '--bins',
        dest='bin_num',
        metavar='INT',
        type=int,
        default=5,
        help="number of histogram bins (default=5)")
    group_figure_options = parser.add_argument_group('figure output options')
    group_figure_options.add_argument(
        '-l',
        dest='length',
        metavar='INT',
        type=int,
        default=7,
        help="length (inches) for figure output (default=7)")
    group_figure_options.add_argument(
        '-w',
        dest='width',
        metavar='INT',
        type=int,
        default=9,
        help="width (inches) for figure output (default=9)")
    group_figure_options.add_argument(
        '--dpi',
        dest='dpi',
        metavar='INT',
        type=int,
        default=300,
        help="dpi for figure output (default=300)")
    group_figure_options.add_argument(
        '--show',
        dest='show_fig',
        action='store_true',
        help="if --show is active, also shows the figure (requires graphical interface) (default=False)")
    
    args = parser.parse_args()

    return args
# --------------------------------------------------
def _add_labels(axis_arg, x: list, y: list) -> None:
    """
    Given a matplotlib axis, add value labels

    Parameters:
        axis_arg: matplotlib.axes._subplots.AxesSubplot
            the matplotlib axis to plot the value labels to
        x: list
            list of x values
        y: list
            list of y values
    
    Returns:
        (None)

    """
    low = int(min(y)*.1)
    for i, _ in enumerate(x):
        val = round(y[i], 2)
        axis_arg.text(i, low, val, ha = 'center', color='w')

    return None
def _generate_bins(args: Namespace, bin_num: int) -> list:
    """
    Given a number of bins, generate a list of bins for histogram plotting

    Parameters:
        args: Namespace
            the arguments produced by argparse, lazy
        bin_num: int
            the number of bins

    Returns:
        (list)
            list of bins for plotting
    """
    step = (args.size_max - args.size_min)/bin_num
    bins = [args.size_min + ((i+1)*step) for i in range(bin_num)]
    bins.insert(0, args.size_min)
    return bins
def _filter_fragments(args: Namespace, fragment_dict_arg: dict):
    """
    Given a dictionary of fragments per chromosome, filter if specified.

    Parameters:
        args: Namespace
            the arguments produced by argparse, lazy
        fragment_dict_arg: dict
            the dictionary of fragment lengths per chromosome
    
    Returns:
        (dict)
            if specified: dictionary of filtered number of fragments per chromosome, else: fragments per chromosome
    """
    total_per_chromosome: dict = {}

    for chromosome in fragment_dict_arg:
        fragment_list = fragment_dict_arg[chromosome]['fragments_list']
        filtered_fragments =\
            [fragment for fragment in fragment_list if args.size_min < fragment < args.size_max] if args.filter \
            else fragment_list
        total_per_chromosome[chromosome] = filtered_fragments

    return total_per_chromosome
def _generate_total_per_chromosome(args: Namespace, fragment_dict_arg: dict, relative: bool = False) -> dict:
    """
    Given a dictionary of fragment lengths, return a dictionary of fragments per chromosome

    Parameters:
        args: Namespace
            the arguments produced by argparse, lazy
        fragment_dict_arg: dict
            the dictionary of fragment lengths per chromosome
        relative: bool
            determine whether to get fragments relative to chromosome length or number of fragments

    Returns:
        (dict)
            dictionary of either representation % or number of fragments per chromosome
    """
    total_per_chromosome: dict = {}
    for chromosome in fragment_dict_arg:
        fragment_list = fragment_dict_arg[chromosome]['fragments_list']
        filtered_fragments =\
            [fragment for fragment in fragment_list if args.size_min < fragment < args.size_max] if args.filter \
            else fragment_list
        total_per_chromosome[chromosome] =\
            sum(filtered_fragments)/fragment_dict_arg[chromosome]['total_len']*100 if relative \
            else int(len(filtered_fragments)/1)
    return total_per_chromosome
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """
    args = get_args()

    fragment_dict: dict = pickle.load(open(args.input_path, 'rb'))
    total_per_chromosome: dict = _generate_total_per_chromosome(args, fragment_dict)
    rel_total_per_chromosome: dict = _generate_total_per_chromosome(args, fragment_dict, relative=True)
    
    # generate plot
    subplots_num = 2 + len(fragment_dict.keys())
    fig, axs = plt.subplots(subplots_num, constrained_layout=True)
    fig.suptitle(args.input_path.stem)

    # plot the axes for genomic representation summary
    for i, fragments_dict in enumerate([total_per_chromosome, rel_total_per_chromosome]):
        axs[i].bar(list(fragments_dict.keys()), list(fragments_dict.values()))
        _add_labels(axs[i], list(fragments_dict.keys()), list(fragments_dict.values()))
    
    # set the axes titles for genomic representation summary
    if not args.filter: axs[0].set_title("Number of fragments")
    else: axs[0].set_title(f"Number of fragments (filtered to {args.size_min} bp - {args.size_max} bp)")
    axs[1].set_title("Representation (%) per chromosome")
    
    # set the plots per chromosome
    for i, chrom in enumerate(fragment_dict.keys()):
        bin_value: list = _generate_bins(args, args.bin_num)
        range_value = (args.size_min, args.size_max) if args.filter else None

        adjusted_ax = axs[i+2]
        
        # plot histogram
        adjusted_ax.hist(_filter_fragments(args, fragment_dict)[chrom], bins=bin_value, range=range_value)

        # share the axes for everything below the genomic representation summary
        if i > 0:
            adjusted_ax.sharex(axs[i+1])
            adjusted_ax.sharey(axs[i+1])
        
        # hide the axes for better spacing
        if i < len(fragment_dict.keys())-1:
            adjusted_ax.set_xticks(bin_value)
            adjusted_ax.axes.tick_params(axis='x', width=0, labelbottom=False)

        # put the chromosome number on the right side so that plots can be compressed
        right_axis = adjusted_ax.twinx()
        right_axis.axes.tick_params(axis='y', width=0, labelright=False)
        right_axis.set_ylabel(chrom, rotation=0, labelpad=10)

    # set size of figure and save
    plt.gcf().set_size_inches(args.length, args.width)
    plt.savefig(args.output_path, dpi=args.dpi)
    if args.show_fig: plt.show()

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()