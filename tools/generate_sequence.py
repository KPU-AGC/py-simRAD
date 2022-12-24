#!/usr/bin/env python3
__description__ =\
"""
Purpose: Generate a random sequence with a given size (--size) and gc content (--gc).
Use '>' to pipe the output of this script into a fasta file. 
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "stable"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
# --------------------------------------------------
import random
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """
    parser = ArgumentParser(
        usage="%(prog)s [options] (> fasta_file.fasta)",
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        '--gc',
        dest='gc_content',
        metavar='%',
        type=float,
        default=0.5,
        help='GC content (%%) (default=50)')
    parser.add_argument(
        '-n',
        '--name',
        dest='fasta_header',
        metavar='str',
        type=str,
        help='name for fasta header (default="# bp sequence with # GC content")')
    parser.add_argument(
        '-w',
        '--width',
        dest='width',
        metavar='N',
        type=int,
        default=None,
        help='character widths for sequence output (default=None)')
    group_sequence_length = parser.add_argument_group('REQUIRED')
    length_synonyms = group_sequence_length.add_mutually_exclusive_group(required=True)
    length_synonyms.add_argument(
        '-l',
        '--length',
        dest='length',
        metavar='bp',
        type=str,
        help="size of sequence (bp), SI prefixes ['k', 'M', 'G'] allowed, \nexamples: 100 = 100, 1G = 1000000000, 1.2Mbp = 1200000")
    length_synonyms.add_argument(
        '-s',
        '--size',
        dest='size',
        metavar='bp',
        type=str,
        help="synonym of --length")

    args = parser.parse_args()
    # parser errors and processing
    # --------------------------------------------------
    if not args.size: args.size = args.length
    if not args.fasta_header: args.fasta_header = f"{args.size} bp sequence with {args.gc_content} GC content"
    if args.gc_content > 1: args.gc_content = args.gc_content/100
    if args.gc_content > 100: parser.error("GC content should be [0-100]")
    if args.size.endswith('bp'): args.size = args.size.replace('bp', '')
    if all([False for char in args.size if char not in "1234567890."]) and args.size.count('.') <= 1: args.size = float(args.size)
    else:
        multiplier = {'k': 1000, 'M': 1000000, 'G': 1e+9}
        prefix = args.size[-1]
        if prefix not in multiplier: parser.error(f"couldn't process {args.size}, check if {prefix} is an SI prefix")
        base = float(args.size[:-1])
        args.size = int(base*multiplier[prefix])

    return args
# --------------------------------------------------
def _generate_sequence(_size: int, _gc_content: float) -> str:
    """
    Return a DNA sequence with a given size (bp) and GC content (%)
    
    Parameters:
        _size: int
            the size/length of the sequence
        _gc_content: float
            GC content to consider when generating the sequence
    
    Return:
        sequence: str
    """
    # choose from list of nucleotides with certain weights and do this k number of times,
    # then join the list together as a single string
    return ''.join(random.choices(['A', 'T', 'C', 'G'], weights=[(1-(_gc_content/2)), (1-(_gc_content/2)), (_gc_content/2), (_gc_content/2)], k=int(_size)))
def main() -> None:
    """ Insert docstring here """
    args = get_args()
    sequence: str = _generate_sequence(args.size, args.gc_content)
    # if width is specified, add newlines whenever the sequence hits the width
    if args.width: sequence: str = ''.join([char if not (i_char+1) % args.width == 0 else char+"\n" for i_char, char in enumerate(sequence)])
    print(f">{args.fasta_header}\n{sequence}")
    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()