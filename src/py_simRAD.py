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
import _export
import _catalysis
import _print_summary
from Bio.Restriction import Restriction
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    
    exploratory_subparser = parser.add_subparsers(
        title='analysis options (choose one)',
        dest='options',
        required=True,
        metavar='',
        description="options for performing or visualizing exploratory restriction digests")

    # --------------------------------------------------
    parser_catalyze = exploratory_subparser.add_parser('catalyze',
        help='generate restriction fragment positions using restriction enzyme combinations (RUN FIRST)')
    parser_catalyze.add_argument('input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser_catalyze.add_argument('--concat', dest='concat',
        action='store_true',
        help="concatenate sequences together in fasta (DEFUNCT)")
    
    group_rst_enz_parser = parser_catalyze.add_argument_group('restriction enzyme options')
    group_rst_enz_parser.add_argument('--enzymes', dest='enzymes',
        metavar='STR',
        type=str,
        default="SbfI;EcoRI;SphI;PstI;MspI;MseI",
        help='enzymes to generate combinations for double-restriction, ex: "EcoRI;BamHI" (default="SbfI;EcoRI;SphI;PstI;MspI;MseI")')
    group_rst_enz_parser.add_argument('--combinations', dest='combinations',
        metavar='INT',
        type=str,
        default="2",
        help='produce combinations of size n (default=2)')
    group_rst_enz_parser.add_argument('--as-is', dest='as_is',
        action='store_true',
        help='use enzymes list as is (double-restrictions allowed); ex: "EcoRI;BamHI;EcoRI-BamHI" (default=False)')
    group_rst_enz_parser.add_argument('--fast', dest='use_fast',
        action='store_true',
        help="use the 'fast' algorithm (but slightly inaccurate; not that much faster) (default=False)")
    group_rst_enz_parser.add_argument('--force', dest='force_new',
        action='store_true',
        help="force program to generate new positions even if previously generated positions exist (default=False)")
    # --------------------------------------------------
    parser_export = exploratory_subparser.add_parser('export',
        help='given a set of restriction fragment positions, export them into other data types')
    parser_export.add_argument('input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser_export.add_argument('output_path',
        type=Path,
        help="path of directory to output the files (see --type argument)")
    parser_export.add_argument('positions_paths',
        type=Path,
        nargs='+',
        help="path of input positions files for output (.pos)")
    
    parser_export.add_argument('-m', '--min', dest='size_min',
        metavar='INT',
        type=int,
        default=0,
        help="min fragment size for filtering (default=None)")
    parser_export.add_argument('-M', '--max', dest='size_max',
        metavar='INT',
        type=int,
        default=10_000,
        help="max fragment size for filtering (default=None)")
    parser_export.add_argument('--select-adapt', dest='select_adapt',
        metavar='str',
        type=str,
        default="AB;BA",
        help="fragment ends for adapter selection (ex: EcoRI + BamHI ends = 'AB;BA') (default='AB;BA')")
    
    group_export_options = parser_export.add_argument_group(title='export options')
    group_export_options.add_argument('--type', dest='export_type',
        metavar='STR',
        type=str,
        choices=['fasta', 'gff'],
        default='fasta',
        required=True,
        help="export types: ['fasta', 'gff'] (default='fasta')")
    # --------------------------------------------------
    parser_genome_rep = exploratory_subparser.add_parser('summary',
        help='given a set of restriction fragment positions, print out percent genomic representation')
    parser_genome_rep.add_argument('positions_paths',
        type=Path,
        nargs='+',
        help="path of input positions files for output")
    
    group_filters = parser_genome_rep.add_argument_group(title='filtering')
    group_filters.add_argument('-m', '--min', dest='size_min',
        metavar='n',
        type=int,
        default=0,
        help="min fragment size (bp) for filtering (default=0)")
    group_filters.add_argument('-M', '--max', dest='size_max',
        metavar='n',
        type=int,
        default=10_000,
        help="max fragment size (bp) for filtering (default=None)")
    group_filters.add_argument('-r', '--min_rep', dest='rep_min',
        metavar='n',
        type=int,
        default=0,
        help="min representation (%%) for filtering (default=0)")
    group_filters.add_argument('-R', '--max_rep', dest='rep_max',
        metavar='n',
        type=int,
        default=10_000,
        help="max representation (%%) for filtering (default=None)")
    group_filters.add_argument('--select-adapt', dest='select_adapt',
        metavar='str',
        type=str,
        default="",
        help="fragment ends for adapter selection (ex: EcoRI + BamHI ends = 'AB;BA') (default='AA' or 'AB;BA' ...)")
    group_filters.add_argument('--exclude', dest='exclude',
        metavar='str',
        type=str,
        default="",
        help="fragments cut by this list are excluded (ex: in a restriction with EcoRI+MseI, excluding MseI, the following are exluded 'MseI-EcoRI;EcoRI-MseI;MseI-MseI') (default='')")

    group_adjusts = parser_genome_rep.add_argument_group(title='adjustments')
    group_adjusts.add_argument( '-p', '--ploidy', dest='ploidy',
        metavar='N',
        type=int,
        default=1,
        help="ploidy to make output slightly more realistic (default=2)")

    group_export_delimiter = parser_genome_rep.add_argument_group(title='format options')
    group_export_delimiter.add_argument('--delimiter', dest='delimiter',
        metavar='STR',
        type=str,
        choices=['tab', ','],
        default='tab',
        help="delimiter for printing (use comma to redirect to .csv): ['tab', ','] (default='tab')")
    group_export_delimiter.add_argument('--type', dest='summary_type',
        metavar='STR',
        type=str,
        choices=['genomic_rep', 'fragment_num'],
        default='fragment_num',
        help="output type: ['genomic_rep', 'fragment_num'] (default='fragment_num')")
    group_export_delimiter.add_argument('--just-totals', dest='chromosome_output',
        action='store_false',
        help="only print totals, not chromosomal breakdown")
    # --------------------------------------------------
    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if args.options=='catalyze':
        if not args.as_is:
            invalid_enzymes: list = [enzyme for enzyme in args.enzymes.split(';') if enzyme not in Restriction.AllEnzymes.elements()]
        else:
            flattened_enzymes: list = sum([enzyme.split('-') for enzyme in args.enzymes.split(';')], [])
            invalid_enzymes: list = [enzyme for enzyme in flattened_enzymes if enzyme not in Restriction.AllEnzymes.elements()]

        if invalid_enzymes: parser.error(f"Couldn't process the following enzymes: {' '.join(invalid_enzymes)}")
        if not args.input_path.is_file():
            parser.error('Invalid input path to a genomic fasta (.fna)!')
        invalid_combinations: list = [num for num in args.combinations.split(';') if num not in "0123456789"]
        if invalid_combinations: parser.error(f"Couldn't process the following enzymes: {' '.join(invalid_enzymes)}")
    return args
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    if args.options=='catalyze': _catalysis._mp_catalysis(args)
    if args.options=='export':
        if args.export_type=='fasta': _export._mp_export_fasta(args)
        if args.export_type=='gff': _export._mp_export_gff(args)
    if args.options=='summary': _print_summary._mp_print_genomic_representation(args)

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()
