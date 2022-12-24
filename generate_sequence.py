#!/usr/bin/env python3
__description__ =\
"""
Purpose: Generate a random sequence with a given size (--size) and gc content (--gc).
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
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        '--size',
        dest='size',
        metavar='bp',
        type=str,
        required=True,
        help="size of sequence (bp), SI prefixes ['k', 'M', 'G'] allowed")
    parser.add_argument(
        '--gc',
        dest='gc_content',
        metavar='%',
        type=float,
        default=0.5,
        help='GC content (%%)')
    parser.add_argument(
        '--name',
        dest='fasta_header',
        metavar='str',
        type=str,
        help='name for fasta header (default: "# bp sequence with # GC content)')

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if not args.fasta_header: args.fasta_header = f"{args.size} bp sequence with {args.gc_content} GC content"
    if args.gc_content > 1: args.gc_content = args.gc_content/100
    if args.gc_content > 100: parser.error("GC content should be [0-100]")
    if all([False for char in args.size if char not in "1234567890."]) and args.size.count('.') <= 1: args.size = float(args.size)
    else:
        multiplier = {'k': 1000, 'M': 1000000, 'G': 1e+9}
        prefix = args.size[-1]
        if prefix not in multiplier: parser.error(f"couldn't process {args.size}, check if {prefix} is an SI prefix")
        base = float(args.size[:-1])
        args.size = base*multiplier[prefix]

    return args
# --------------------------------------------------
def _generate_sequence(_gc_content: float, _size) -> str:
    """"""
    return ''.join(random.choices(['A', 'T', 'C', 'G'], weights=[(1-(_gc_content/2)), (1-(_gc_content/2)), (_gc_content/2), (_gc_content/2)], k=int(_size)))
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()
    sequence = _generate_sequence(args.gc_content, args.size)
    print(f">{args.fasta_header}")
    print(sequence)

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()