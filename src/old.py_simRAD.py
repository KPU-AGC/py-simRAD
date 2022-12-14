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
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import Restriction
from itertools import combinations, product
from multiprocessing import Pool
import pickle
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

    parser_catalyze = exploratory_subparser.add_parser('catalyze', help='generate restriction fragment positions using restriction enzyme combinations (RUN FIRST)')
    parser_catalyze.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser_catalyze.add_argument(
        '--concat',
        dest='concat',
        action='store_true',
        help="concatenate sequences together in fasta (DEFUNCT)")
    group_rst_enz_parser = parser_catalyze.add_argument_group('restriction enzyme options')
    group_rst_enz_parser.add_argument(
        '--enzymes',
        dest='enzymes',
        metavar='STR',
        type=str,
        default="SbfI;EcoRI;SphI;PstI;MspI;MseI",
        help='enzymes to generate combinations for double-restriction, ex: "EcoRI;BamHI" (default="SbfI;EcoRI;SphI;PstI;MspI;MseI")')
    group_rst_enz_parser.add_argument(
        '--combinations',
        dest='combinations',
        metavar='INT',
        type=str,
        default="2",
        help='produce combinations of size n (default=2)')
    group_rst_enz_parser.add_argument(
        '--as-is',
        dest='as_is',
        action='store_true',
        help='use enzymes list as is (double-restrictions allowed); ex: "EcoRI;BamHI;EcoRI-BamHI" (default=False)')
    group_rst_enz_parser.add_argument(
        '--fast',
        dest='use_fast',
        action='store_true',
        help="use the 'fast' algorithm (but slightly inaccurate; not that much faster) (default=False)")
    group_rst_enz_parser.add_argument(
        '--force',
        dest='force_new',
        action='store_true',
        help="force program to generate new positions even if previously generated positions exist (default=False)")

    parser_export = exploratory_subparser.add_parser('export', help='given a set of restriction fragment positions, export them into other data types')
    parser_export.add_argument(
        'input_path',
        type=Path,
        help="path of input genomic file (.fna)")
    parser_export.add_argument(
        'output_path',
        type=Path,
        help="path of directory to output the files (see --type argument)")
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
    parser_export.add_argument(
        '--select-adapt',
        dest='select_adapt',
        metavar='str',
        type=str,
        default="AB;BA",
        help="fragment ends for adapter selection (ex: EcoRI + BamHI ends = 'AB;BA') (default='AB;BA')")
    group_export_options = parser_export.add_argument_group(title='export options')
    group_export_options.add_argument(
        '--type',
        dest='export_type',
        metavar='STR',
        type=str,
        choices=['fasta', 'gff'],
        default='fasta',
        required=True,
        help="export types: ['fasta', 'gff'] (default='fasta')")

    parser_genome_rep = exploratory_subparser.add_parser('summary', help='given a set of restriction fragment positions, print out percent genomic representation')
    parser_genome_rep.add_argument(
        'positions_paths',
        type=Path,
        nargs='+',
        help="path of input positions files for output")
    group_filters = parser_genome_rep.add_argument_group(title='filtering')
    group_filters.add_argument(
        '-m',
        '--min',
        dest='size_min',
        metavar='n',
        type=int,
        default=0,
        help="min fragment size (bp) for filtering (default=0)")
    group_filters.add_argument(
        '-M',
        '--max',
        dest='size_max',
        metavar='n',
        type=int,
        default=None,
        help="max fragment size (bp) for filtering (default=None)")
    group_filters.add_argument(
        '-r',
        '--min_rep',
        dest='rep_min',
        metavar='n',
        type=int,
        default=0,
        help="min representation (%%) for filtering (default=0)")
    group_filters.add_argument(
        '-R',
        '--max_rep',
        dest='rep_max',
        metavar='n',
        type=int,
        default=None,
        help="max representation (%%) for filtering (default=None)")
    group_filters.add_argument(
        '--select-adapt',
        dest='select_adapt',
        metavar='str',
        type=str,
        default="AB;BA",
        help="fragment ends for adapter selection (ex: EcoRI + BamHI ends = 'AB;BA') (default='AB;BA')")
    group_filters.add_argument(
        '--exclude',
        dest='exclude',
        metavar='str',
        type=str,
        default="",
        help="fragments cut by this list are excluded (ex: in a restriction with EcoRI+MseI, excluding MseI, the following are exluded 'MseI-EcoRI;EcoRI-MseI;MseI-MseI') (default='')")

    group_adjusts = parser_genome_rep.add_argument_group(title='adjustments')
    group_adjusts.add_argument(
        '-p',
        '--ploidy',
        dest='ploidy',
        metavar='N',
        type=int,
        default=1,
        help="ploidy to make output slightly more realistic (default=2)")

    group_export_delimiter = parser_genome_rep.add_argument_group(title='format options')
    group_export_delimiter.add_argument(
        '--delimiter',
        dest='delimiter',
        metavar='STR',
        type=str,
        choices=['tab', ','],
        default='tab',
        help="delimiter for printing (use comma to redirect to .csv): ['tab', ','] (default='tab')")
    group_export_delimiter.add_argument(
        '--type',
        dest='summary_type',
        metavar='STR',
        type=str,
        choices=['genomic_rep', 'fragment_num'],
        default='fragment_num',
        help="output type: ['genomic_rep', 'fragment_num'] (default='fragment_num')")

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
def _generate_restriction_fragments3(_input_seq: SeqRecord, _restriction_batch: Restriction.RestrictionBatch) -> dict:
    """
    From a list of cutting positions, do the cutting and generate a list of fragment sizes,
    but this way is more biologically accurate, technically. This is meant to handle multiple
    sizes of restrictions batches -- i.e., triple restriction is technically doable.

    Parameters:
        _input_seq: Seq
            chromosomal sequence to cut
        _restriction_batch: RestrictionBatch
            a set of enzymes to cut it with

    Returns:
        (dict)
            list of restriction restriction fragments of given lengths
    """

    enzymes = [i for i in _restriction_batch]

    intermediate_fragments: list = [_input_seq.upper()]
    intermediate_enzymes: list = [("None", "None")]
    for i_enzyme, enzyme in enumerate(enzymes):
        intermediate_fragment_results = []
        intermediate_enzyme_results = []
        for i_int_frag, intermediate_fragment in enumerate(intermediate_fragments):
            catalysis_results = enzyme.catalyse(intermediate_fragment)
            

            intermediate_enzyme_tags = []
            for i_enzyme_frag, _ in enumerate(catalysis_results):
                left_enzyme = str(enzyme)
                right_enzyme = str(enzyme)
                
                if i_enzyme_frag == 0: left_enzyme = intermediate_enzymes[i_int_frag][0]
                if i_enzyme_frag == len(catalysis_results)-1: right_enzyme = intermediate_enzymes[i_int_frag][1]
                intermediate_enzyme_tags.append((left_enzyme, right_enzyme))
            
            intermediate_enzyme_results += intermediate_enzyme_tags
            intermediate_fragment_results += catalysis_results
        intermediate_fragments = intermediate_fragment_results
        intermediate_enzymes = intermediate_enzyme_results
    
    adjusted_fragment_positions = []
    fragment_lengths = [0] + [len(fragment) for fragment in intermediate_fragments]
    i = 0
    while i < len(fragment_lengths):
        if i == 0: previous_value = 0
        else: previous_value = adjusted_fragment_positions[i-1][1]
        adjusted_fragment_positions.append((previous_value, previous_value + fragment_lengths[i]))
        i += 1
    adjusted_fragment_positions.pop(0)

    return {'fragment_positions': adjusted_fragment_positions, 'enzyme_end_positions': intermediate_enzymes}
def _generate_restriction_fragments_fast(_slice_positions: list) -> dict:
        """
        From a list of cutting positions, do the cutting and generate a list of fragment sizes.

        Parameters:
            _slice_positions: list
                list of slice positions from the restriction enzymes

        Returns:
            (list)
                list of restriction restriction fragments of given lengths
        """
        fragment_positions: list = []
        len_slice_positions: int = len(_slice_positions)-1
        for i_pos, _ in enumerate(_slice_positions):
            if i_pos < len_slice_positions:
                start_pos = _slice_positions[i_pos]
                end_pos = _slice_positions[i_pos + 1]
                fragment_positions.append((start_pos, end_pos))
        return {'fragment_positions': fragment_positions}
def _mp_catalysis(args: Namespace) -> None:
    """ Multiprocess wrapper for the _catalysis function. Also generates a list of enzymes for processing. """

    output_dir = args.input_path.parent.joinpath('.'+str(args.input_path.stem).replace(' ', '_'))
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # set list of restriction enzymes to query
    # generate restriction fragment combinations, including singletons
    if not args.as_is:
        restriction_enzymes_list: list = [enzyme.strip() for enzyme in args.enzymes.split(';') if enzyme]
        rst_enz_combinations: list = []
        for i in ([int(num) for num in args.combinations.split(';')]): rst_enz_combinations += [enzyme_pair for enzyme_pair in combinations(restriction_enzymes_list, i)]
        rst_enz_combinations: list = sorted(rst_enz_combinations)
    # generate a list of restriction fragment (singletons/pairs) for analysis
    else:
        restriction_enzymes_list: list = [tuple(enzyme_pair.split('-')) for enzyme_pair in args.enzymes.split(';')]
        rst_enz_combinations: list = sorted(restriction_enzymes_list)

    map_offset = len(rst_enz_combinations)
    map_args = tuple(zip(
        [args.input_path]*map_offset,
        rst_enz_combinations,
        [output_dir]*map_offset,
        [args.use_fast]*map_offset,
        [args.force_new]*map_offset,
        ))
    with Pool() as pool:
        pool.starmap(_catalysis, map_args)

    return None
def _catalysis(_input_path: Path, _enzyme_combination: tuple, _output_dir: Path, _use_fast: bool, _force_new: bool) -> None:
    """
    Function performs catalysis on a given genome with a given cocmbination of restriction enzymes.

    Parameters:
        _input_path: Path
            the path of the genomic fasta (.fna)
        _enzyme_combination: tuple
            tuple containing the combination of restriction enzymes
        _output_dir: Path
            the path for output
        _use_fast: bool
            determines whether the fast algorithm is used
    
    Returns
        (None)
    """

    _method = 'fast' if _use_fast else 'comprehensive'
    _metadata = {
        'metadata': {
            'version': __version__,
            'method': _method,
        }
    }

    combination_str: str = '-'.join(list(_enzyme_combination))
    output_file: Path = _output_dir.joinpath(f'{combination_str}.pos')

    # check if there is a previous output that exists and matches the same version and method,
    # and don't bother redoing it if so.

    if output_file.exists() and not _force_new:
        previous_session = pickle.load(open(output_file, 'rb'))
        if all([
                previous_session['metadata']['version'].split('.')[:2] == __version__.split('.')[:2],   # check for MAJOR.MINOR, not PATCH
                previous_session['metadata']['method'] == _method,]):                                   # check that method is the same
            return None

    fragments_per_chrom: dict = {}
    fragments_per_chrom.update(_metadata)

    for chr in SeqIO.parse(_input_path, 'fasta'):
        # ignore the mitochrondrial genome, only do nuclear genome
        if 'mitochondrion' in chr.description: continue
        
        restriction_fragments_positions: dict = {}
        if _use_fast:
            # get all unique slice positions
            restriction_batch = Restriction.RestrictionBatch(list(_enzyme_combination))
            restriction_result = restriction_batch.search(chr.seq.upper())
            # using (end - beginning), add 0 position so that first fragment is the length from the start to that position
            # also add end position -- the entire length of the chromosome
            # and remove duplicate positions from isoschizomers or similar cut sites
            slice_positions: list = sorted(set([0] + [slice_pos for _, enzyme_slice_list in restriction_result.items() for slice_pos in enzyme_slice_list] + [len(chr.seq)]))

            # generate "fragments" by cutting between slice positions
            restriction_fragments_positions = _generate_restriction_fragments_fast(
                    _slice_positions=slice_positions)
        elif not _use_fast:
            restriction_enzymes = Restriction.RestrictionBatch(list(_enzyme_combination))
            restriction_fragments_positions = _generate_restriction_fragments3(
                _input_seq=chr.seq,
                _restriction_batch=restriction_enzymes)
        
        fragments_per_chrom[chr.id] = restriction_fragments_positions
    # output fragments per chromosome
    pickle.dump(fragments_per_chrom, open(output_file, 'wb'))
    return None
# --------------------------------------------------
def _mp_export_fasta(args: Namespace) -> None:
    """ Multiprocess wrapper for the _export_fasta function. """
    map_offset = len(args.positions_paths)
    map_args = tuple(zip(
        [(args.input_path)]*map_offset,
        args.positions_paths,
        [(args.output_path)]*map_offset,
        [(args.size_min, args.size_max)]*map_offset
        ))
    with Pool() as pool: pool.starmap(_export_fasta, map_args)
    return None
def _export_fasta(_input_path: Path, _positions_path: Path, _output_path: Path, _size_filter: tuple) -> None:
    """
    Function outputs a fasta (.fasta) file of (all/filtered) fragments generated from a given genome.

    Parameters:
        _input_path: Path
            path of a given genomic fasta (.fna)
        _positions_path: Path
            path of restriction fragment position data (.pos)
        _output_path: Path
            path of output directory for fasta (.fasta) sequences
        _size_filter: tuple(min_size, max_size)
            tuple of filtering options for each fragment
    
    Returns
        (None)
    """

    positions_dict: dict = pickle.load(open(_positions_path, 'rb'))

    list_of_seqs: list = []
    for chr in SeqIO.parse(_input_path, 'fasta'):
        if 'mitochondrion' in chr.description: continue
        
        for position in positions_dict[chr.id]['fragment_positions']:

            fragment_length = position[1] - position[0]
            if fragment_length < _size_filter[0]: continue
            if _size_filter[1]: 
                if fragment_length > _size_filter[1]: continue

            fragment_SeqRecord = SeqRecord(
                seq=chr.seq[position[0]:position[1]],
                id=chr.id,
                name='',
                description=chr.description + f' {_positions_path.stem} {position[0]}-{position[1]}')
            list_of_seqs.append(fragment_SeqRecord)
    SeqIO.write(list_of_seqs, _output_path.joinpath(f'{_positions_path.stem}.fasta'), 'fasta')
    return None
# --------------------------------------------------
def _mp_export_gff(args: Namespace) -> None:
    """ Multiprocess wrapper for the _export_gff function. """
    map_offset = len(args.positions_paths)
    map_args = tuple(zip(
        [(args.input_path)]*map_offset,
        args.positions_paths,
        [(args.output_path)]*map_offset,
        [(args.size_min, args.size_max)]*map_offset,
        [(args.select_adapt)]*map_offset,
        ))
    with Pool() as pool: pool.starmap(_export_gff, map_args)
    return None
def _export_gff(_input_path: Path, _positions_path: Path, _output_path: Path, _size_filter: tuple, _select_adapt_str: str) -> None:
    """
    Function outputs a general feature format (.gff) file to describe the fragments generated from a given genome.

    Parameters:
        _input_path: Path
            path of a given genomic fasta (.fna)
        _positions_path: Path
            path of restriction fragment position data (.pos)
    
    Returns:
        (None)
    """

    _select_adapt_allowed = _select_adapt_str.split(';')
    enzyme_to_letter_conversion = {enzyme: letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')}

    with open(_output_path.joinpath(f'{_positions_path.stem}.gff'), mode='w', encoding='utf-8') as output_gff:
        positions_dict: dict = pickle.load(open(_positions_path, 'rb'))

        output_gff.write('##gff-version 3\n')
        for chr in SeqIO.parse(_input_path, 'fasta'):
            if 'mitochondrion' in chr.description: continue
            
            for fragment_position, enzyme_ends in zip(positions_dict[chr.id]['fragment_positions'], positions_dict[chr.id]['enzyme_end_positions']):
                fragment_length = fragment_position[1] - fragment_position[0]
                if fragment_length < _size_filter[0]: continue
                if _size_filter[1]: 
                    if fragment_length > _size_filter[1]: continue
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
                gff_str = f"{chr.id}\tpy-simRADv{__version__}\trestriction_fragment\t{fragment_position[0]}\t{fragment_position[1]}\t.\t+\t0\t.\n"
                output_gff.write(gff_str)
    return None
# --------------------------------------------------
def _mp_print_genomic_representation(args: Namespace) -> None:
    """ Multiprocess wrapper for the _export_gff function. """

    chr_list = [key for key in pickle.load(open(args.positions_paths[0], 'rb')).keys() if key != 'metadata']

    delimiter = '\t' if args.delimiter == 'tab' else ','

    map_offset = len(args.positions_paths)
    
    map_args = tuple(zip(
        args.positions_paths,
        [(args.ploidy)]*map_offset,
        [(args.size_min, args.size_max)]*map_offset,
        [(args.rep_min, args.rep_max)]*map_offset,
        [args.summary_type]*map_offset,
        [(args.select_adapt)]*map_offset,
        ))

    with Pool() as pool: 
        genomic_representation_list: list = pool.starmap(_print_genomic_representation, map_args)

    total_val: str = ""
    if args.summary_type == "genomic_rep": total_val = 'total repr (%)'
    elif args.summary_type == "fragment_num": total_val = "total number of fragments"

    headers: list = ['enzmye', total_val] + chr_list
    print(delimiter.join(headers))
    for row in sorted([row for row in genomic_representation_list if row[0]], key=lambda x: x[1]):
        row = [str(item) for item in row]
        print(delimiter.join(row))
    return None
def _print_genomic_representation(_positions_path: Path, _ploidy:int, _size_filters: tuple, _representation_filter: tuple, _print_type: str, _select_adapt_str: str) -> list:
    """
    Print genomic representation (%) to the console, given whatever filtering options.

    Parameters:
        _positions_path: Path
            path of restriction fragment position data (.pos)
        _size_filter: tuple(min_size, max_size)
            tuple of filtering options for each fragment
        _representation_filter: tuple(min_rep, max_rep)
            tuple of filtering options for genomic representation percentage
    
    Returns:
        (list)
            list of [enzyme combination, percent genome represented] + [representation per chromosome]
    """
    positions_dict: dict = pickle.load(open(_positions_path, 'rb'))

    _select_adapt_allowed = _select_adapt_str.split(';')
    enzyme_to_letter_conversion = {enzyme: letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')}
    
    def _genomic_representation() -> list:
        total_genome_length: int = 0
        genome_represented: int = 0
        per_chromosome_rep_list: list = []
        
        chromosomes_list: list = [key for key in positions_dict.keys() if key != 'metadata']

        for chr in chromosomes_list:
            chr_length = positions_dict[chr]['fragment_positions'][-1][1]
            total_genome_length += chr_length
        
            chromosome_represented: int = 0
            for fragment_position, enzyme_ends in zip(positions_dict[chr]['fragment_positions'], positions_dict[chr]['enzyme_end_positions']):
                fragment_length = fragment_position[1] - fragment_position[0]
                if fragment_length < _size_filters[0]: continue
                if _size_filters[1]: 
                    if fragment_length > _size_filters[1]: continue
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
                chromosome_represented += fragment_length

            genome_represented += chromosome_represented
            perc_chromosome_represented: float = round(chromosome_represented/chr_length*100, 3)
            per_chromosome_rep_list.append(perc_chromosome_represented)

        perc_genome_represented: float = round(genome_represented/total_genome_length*100, 3)
        if perc_genome_represented < _representation_filter[0]: return [None]
        if _representation_filter[1]: 
            if perc_genome_represented > _representation_filter[1]: return [None]

        return [_positions_path.stem, perc_genome_represented] + per_chromosome_rep_list
    def _fragments_per_chromosome() -> list:
        """"""

        fragments_per_chromosome = []
        for chr in positions_dict:
            if chr == 'metadata': continue
            fragments_passing_filter = 0
            for fragment_position, enzyme_ends in zip(positions_dict[chr]['fragment_positions'], positions_dict[chr]['enzyme_end_positions']):
                fragment_length = fragment_position[1] - fragment_position[0]
                if fragment_length < _size_filters[0]: continue
                if _size_filters[1]: 
                    if fragment_length > _size_filters[1]: continue
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
                fragments_passing_filter += 1
            fragments_per_chromosome.append(int(fragments_passing_filter*_ploidy))
        total_fragment_count: int = sum(fragments_per_chromosome)
        return [_positions_path.stem, total_fragment_count, *fragments_per_chromosome]

    if _print_type == 'genomic_rep': return _genomic_representation()
    elif _print_type == 'fragment_num': return _fragments_per_chromosome()
    else: return []
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    import timeit

    args = get_args()

    if args.options=='catalyze': _mp_catalysis(args)
    if args.options=='export':
        if args.export_type=='fasta': _mp_export_fasta(args)
        if args.export_type=='gff': _mp_export_gff(args)
    if args.options=='summary': _mp_print_genomic_representation(args)

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()
