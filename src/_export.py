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
from multiprocessing import Pool
import pickle
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
                if not _size_filter[0] < fragment_length < _size_filter[1]: continue
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
                gff_str = f"{chr.id}\tpy-simRADv{__version__}\trestriction_fragment\t{fragment_position[0]}\t{fragment_position[1]}\t.\t+\t0\t.\n"
                output_gff.write(gff_str)
    return None
# --------------------------------------------------