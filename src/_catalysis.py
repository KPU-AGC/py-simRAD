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
from itertools import combinations
from multiprocessing import Pool
import pickle
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