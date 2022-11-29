#!/usr/bin/env python3
__description__ =\
"""
Purpose: Perform double restriction digest on a given genome.
"""
__author__ = "Erick Samera; Michael Ke"
__version__ = "3.0.0"
__comments__ = "stable; better implementation"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
from itertools import combinations
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import Restriction
import pandas as pd
import pickle
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)
    
    exploratory_subparser = parser.add_subparsers(
        title='analysis options',
        dest='options',
        required=True,
        description="options for performing or visualizing exploratory restriction digests")

    parser_catalyze = exploratory_subparser.add_parser('catalyze', help='Given a genome, generate restriction fragments using combinations of restriction enzymes (RUN FIRST)')
    parser_catalyze.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    group_rst_enz_parser = parser_catalyze.add_argument_group('restriction enzyme options')
    group_rst_enz_parser.add_argument(
        '--enzymes',
        dest='enzymes',
        metavar='STR',
        type=str,
        default="SbfI;EcoRI;SphI;PstI;MspI;MseI",
        help='enzymes to generate combinations for double-restriction, ex: "EcoRI;BamHI" (default: "SbfI;EcoRI;SphI;PstI;MspI;MseI")')
    group_rst_enz_parser.add_argument(
        '--no_combination',
        dest='no_combination',
        action='store_true',
        help='do not generate combinations for double-restriction, use enzyme list as is')
    group_rst_enz_parser.add_argument(
        '--fast',
        dest='use_fast',
        action='store_true',
        help="use the fast algorithm (but slightly inaccurate)")

    parser_export = exploratory_subparser.add_parser('export', help='Given a set of restriction fragment positions, export them into other data types')
    parser_export.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser_export.add_argument(
        'output_path',
        type=Path,
        help="path of output file (see --type argument)")
    parser_export.add_argument(
        'positions_paths',
        type=Path,
        nargs='+',
        help="path of input positions files for output (.pos)")
    parser_export.add_argument(
        '-m',
        '--min',
        dest='size_min',
        metavar='INT',
        type=int,
        default=0,
        help="min fragment size for filtering (default=None)")
    parser_export.add_argument(
        '-M',
        '--max',
        dest='size_max',
        metavar='INT',
        type=int,
        default=None,
        help="max fragment size for filtering (default=None)")
    group_export_options = parser_export.add_argument_group(title='asdfasdf')
    group_export_options.add_argument(
        '--type',
        dest='export_type',
        metavar='STR',
        type=str,
        choices=['fasta', 'gff'],
        default='fasta',
        help="export types: ['fasta', 'gff'] (default='fasta')")

    parser_genome_rep = exploratory_subparser.add_parser('genome_rep', help='Given a set of restriction fragment positions, print out percent genomic representation')
    parser_genome_rep.add_argument(
        'positions_paths',
        type=Path,
        nargs='+',
        help="path of input positions files for output")
    parser_genome_rep.add_argument(
        '-m',
        '--min',
        dest='size_min',
        metavar='INT',
        type=int,
        default=0,
        help="min fragment size for filtering (default=None)")
    parser_genome_rep.add_argument(
        '-M',
        '--max',
        dest='size_max',
        metavar='INT',
        type=int,
        default=None,
        help="max fragment size for filtering (default=None)")
    parser_genome_rep.add_argument(
        '-r',
        '--min_rep',
        dest='rep_min',
        metavar='INT',
        type=int,
        default=0,
        help="min fragment size for filtering (default=None)")
    parser_genome_rep.add_argument(
        '-R',
        '--max_rep',
        dest='rep_max',
        metavar='INT',
        type=int,
        default=None,
        help="max fragment size for filtering (default=None)")
    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if args.options=='catalyze':
        if not args.no_combination:
            invalid_enzymes: list = [enzyme for enzyme in args.enzymes.split(';') if enzyme not in Restriction.AllEnzymes.elements()]
        else:
            flattened_enzymes: list = sum([enzyme.split('-') for enzyme in args.enzymes.split(';')], [])
            invalid_enzymes: list = [enzyme for enzyme in flattened_enzymes if enzyme not in Restriction.AllEnzymes.elements()]

        if invalid_enzymes: parser.error(f"Couldn't process the following enzymes: {' '.join(invalid_enzymes)}")
    return args
# --------------------------------------------------
def _perform_catalysis(args: Namespace) -> None:
    """

    """

    def _generate_restriction_fragments(seq_arg: SeqRecord, restriction_enzymes_arg: Restriction.RestrictionBatch) -> dict:
        """
        From a list of cutting positions, do the cutting and generate a list of fragment sizes, \
        but this way is more biologically accurate, technically. It cuts with the first enzyme, \
        then cuts again with the second enzyme (if applicable).

        Parameters:
            seq_arg: Seq
                chromosomal sequence to cut
            restriction_enzymes_arg: RestrictionBatch
                a set of enzymes to cut it with
            buffer_arg: int
                biological buffer, assuming that enzyme can't latch on to edge sites
            fasta_output_path_arg: path
                path to output fastas if given
        
        Returns:
            (dict)
                list of restriction restriction fragments of given lengths
        """

        enzymes = [i for i in restriction_enzymes_arg]
        first_restriction = enzymes[0]
        second_restriction = enzymes[1] if len(enzymes)>1 else None

        # generate restriction fragments using the first enzyme,
        first_pass = first_restriction.catalyse(seq_arg)

        total_positions = []
        first_pass_offsets: list = [0]
        for fragment in first_pass:
            first_pass_offsets.append(first_pass_offsets[-1] + len(fragment))
        
        if second_restriction:
            for i, fragment in enumerate(first_pass):
                first_pass_offset = first_pass_offsets[i]
                second_restriction_results = second_restriction.catalyse(fragment)

                second_pass_offsets: list = [first_pass_offset]
                for second_fragment in second_restriction_results:
                    second_pass_offsets.append(second_pass_offsets[-1] + len(second_fragment))
                second_pass_offsets.pop()
                total_positions += second_pass_offsets
            total_positions.append(len(seq_arg))

        else:
            first_pass_offsets.append(len(seq_arg))
            total_positions = first_pass_offsets

        return _generate_restriction_fragments_fast(sorted(set(total_positions)))
    def _generate_restriction_fragments_fast(slice_positions_arg: list) -> dict:
        """
        From a list of cutting positions, do the cutting and generate a list of fragment sizes.

        Parameters:
            slice_positions_arg: list
                list of slice positions from the restriction enzymes

        Returns:
            (list)
                list of restriction restriction fragments of given lengths
        """
        fragment_positions: list = []
        len_slice_positions: int = len(slice_positions_arg)-1
        for i_pos, _ in enumerate(slice_positions_arg):
            if i_pos < len_slice_positions:
                start_pos = slice_positions_arg[i_pos]
                end_pos = slice_positions_arg[i_pos + 1]
                fragment_positions.append((start_pos, end_pos))
        return {'fragment_positions': fragment_positions}

    _method = 'fast' if args.use_fast else 'standard'
    
    output_dir = args.input_path.parent.joinpath('.'+str(args.input_path.stem).replace(' ', '_'))
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # set list of restriction enzymes to query
    if not args.no_combination:
        restriction_enzymes_list: list = args.enzymes.split(';') if isinstance(args.enzymes, str) else args.enzymes
        rst_enz_combinations: list = \
        sorted([enzyme_pair for enzyme_pair in combinations(restriction_enzymes_list, 2)] +\
        [(enzyme, ) for enzyme in restriction_enzymes_list])
    else:
        restriction_enzymes_list: list = [tuple(enzyme_pair.split('-')) for enzyme_pair in args.enzymes.split(';')]
        rst_enz_combinations: list = sorted(restriction_enzymes_list)

    for combination in rst_enz_combinations:
        combination_str = '-'.join(list(combination))

        fragments_per_chrom: dict = {
            'metadata': {
                'version': __version__,
                'method': _method
            }
        }
        total = 0

        for chr in SeqIO.parse(args.input_path, 'fasta'):
            # ignore the mitochrondrial genome, only do nuclear genome
            if 'mitochondrion' in chr.description: continue

            # generate a list of fragments per chromosome
            #chr_str: str = f'chr{chr.description[-1]}'
            fragments_per_chrom[chr.id] = {
                'fragments_list': []
                }
            
            restriction_fragments_positions: dict = {}
            if args.use_fast:
                # get all unique slice positions
                restriction_batch = Restriction.RestrictionBatch(list(combination))
                restriction_result = restriction_batch.search(chr.seq.upper())
                # using (end - beginning), add 0 position so that first fragment is the length from the start to that position
                # also add end position -- the entire length of the chromosome
                # and remove duplicate positions from isoschizomers or similar cut sites
                slice_positions: list = sorted(set([0] + [slice_pos for _, enzyme_slice_list in restriction_result.items() for slice_pos in enzyme_slice_list] + [len(chr.seq)]))

                # generate "fragments" by cutting between slice positions
                restriction_fragments_positions = _generate_restriction_fragments_fast(
                        slice_positions_arg=slice_positions)
            elif not args.use_fast:
                restriction_enzymes = Restriction.RestrictionBatch(list(combination))
                restriction_fragments_positions = _generate_restriction_fragments(
                    seq_arg=chr.seq,
                    restriction_enzymes_arg=restriction_enzymes)
            
            fragments_per_chrom[chr.id] = restriction_fragments_positions
            total += len(*restriction_fragments_positions.values())
        pickle.dump(fragments_per_chrom, open(output_dir.joinpath(f'{combination_str}.pos'), 'wb'))

    return None
def _export_fasta(args: Namespace) -> None:
    """
    """

    for file in args.positions_paths:
        positions_dict: dict = pickle.load(open(file, 'rb'))

        list_of_seqs: list = []
        for chr in SeqIO.parse(args.input_path, 'fasta'):
            if 'mitochondrion' in chr.description: continue
            
            for position in positions_dict[chr.id]['fragment_positions']:

                fragment_length = position[1] - position[0]
                if fragment_length < args.size_min: continue
                if args.size_max: 
                    if fragment_length > args.size_max: continue

                fragment_SeqRecord = SeqRecord(
                    seq=chr.seq[position[0]:position[1]],
                    id=chr.id,
                    name='',
                    description=chr.description + f' {file.stem} {position[0]}-{position[1]}')
                list_of_seqs.append(fragment_SeqRecord)
        SeqIO.write(list_of_seqs, args.output_path.joinpath(f'{file.stem}.fasta'), 'fasta')
def _export_gff(args: Namespace) -> None:
    """
    """
    for file in args.positions_paths:
        with open(args.output_path.joinpath(f'{file.stem}.gff'), mode='w', encoding='utf-8') as output_gff:
            positions_dict: dict = pickle.load(open(file, 'rb'))

            output_gff.write('##gff-version 3\n')
            for chr in SeqIO.parse(args.input_path, 'fasta'):
                if 'mitochondrion' in chr.description: continue
                
                for position in positions_dict[chr.id]['fragment_positions']:
                    gff_str = f"{chr.id}\tpy-simRAD\trestriction_fragment\t{position[0]}\t{position[1]}\t.\t+\t.\n"
                    output_gff.write(gff_str)
    return None
def _print_genomic_representation(args: Namespace) -> None:
    """
    """

    chr_list = [key for key in pickle.load(open(args.positions_paths[0], 'rb')).keys() if key != 'metadata']
    chr_header: str = '\t'.join(chr_list)

    print(f'enzyme\ttotal repr (%)\t{chr_header}')

    for file in args.positions_paths:
        positions_dict: dict = pickle.load(open(file, 'rb'))

        total_genome_length: int = 0
        genome_represented: int = 0
        per_chromosome_rep_list: list = []
        
        chromosomes_list: list = [key for key in positions_dict.keys() if key != 'metadata']

        for chr in chromosomes_list:
            
            chr_length = positions_dict[chr]['fragment_positions'][-1][1]
            total_genome_length += chr_length
        
            chromosome_represented: int = 0
            for position in positions_dict[chr]['fragment_positions']:

                fragment_length = position[1] - position[0]
                if fragment_length < args.size_min: continue
                if args.size_max: 
                    if fragment_length > args.size_max: continue
                chromosome_represented += fragment_length

            genome_represented += chromosome_represented
            
            perc_chromosome_represented: float = round(chromosome_represented/chr_length*100, 3)
            per_chromosome_rep_list.append(perc_chromosome_represented)
        perc_genome_represented: float = round(genome_represented/total_genome_length*100, 3)

        if perc_genome_represented < args.rep_min: continue
        if args.rep_max: 
            if perc_genome_represented > args.rep_max: continue
        per_chromosome_str: str = '\t'.join([f'{i}' for i in per_chromosome_rep_list])
        print(f'{file.stem}\t{perc_genome_represented}\t{per_chromosome_str}')
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    if args.options=='catalyze': _perform_catalysis(args)
    if args.options=='export':
        if args.export_type=='fasta': _export_fasta(args)
        if args.export_type=='gff': _export_gff(args)
    if args.options=='genome_rep': _print_genomic_representation(args)

# --------------------------------------------------
if __name__ == '__main__':
    main()