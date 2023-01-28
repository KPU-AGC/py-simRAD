#!/usr/bin/env python3
__description__ =\
"""
Purpose: Perform double restriction digest on a given genome.
"""
__author__ = "Erick Samera; Michael Ke"
__version__ = "5.1.0"
__comments__ = "stable; multi-processing"
# --------------------------------------------------
from argparse import Namespace
from pathlib import Path
# --------------------------------------------------
from itertools import product, permutations
from multiprocessing import Pool
import pickle
# --------------------------------------------------

def _mp_print_summary(args: Namespace) -> None:
    """ Multiprocess wrapper for the print summary function. """

    chr_list = [key for key in pickle.load(open(args.positions_paths[0], 'rb')).keys() if key != 'metadata']

    delimiter = '\t' if args.delimiter == 'tab' else ','

    map_offset = len(args.positions_paths)
    map_args = tuple(zip(
        args.positions_paths,
        [(args.ploidy)]*map_offset,
        [(args.size_min, args.size_max)]*map_offset,
        [(args.rep_min, args.rep_max)]*map_offset,
        [(args.select_adapt)]*map_offset,
        ))

    summary_functions = {
        'genomic_rep': {
            'function': _print_genomic_representation,
            'headers': []},
        'fragment_num': {
            'function': _print_fragment_num
            },
    }

    with Pool() as pool: 
        summary_to_print: list = pool.starmap(summary_functions[args.summary_type]['function'], map_args)

    total_val: str = ""
    if args.summary_type == "genomic_rep": total_val = 'total repr (%)'
    elif args.summary_type == "fragment_num": total_val = "total number of fragments"

    headers: list = ['enzmye', total_val] + chr_list if args.chromosome_output else ['enzmye', total_val]
    print(delimiter.join(headers))
    for row in sorted([row for row in summary_to_print if row[0]], key=lambda x: x[1]):
        row = [str(item) for item in row] if args.chromosome_output else [str(item) for item in row][:2]
        print(delimiter.join(row))
    return None

def _print_genomic_representation(_positions_path: Path, _ploidy:int, _size_filter: tuple, _representation_filter: tuple, _select_adapt_str: str) -> list:
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

    if _select_adapt_str:
        letters = [letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')]
        _select_adapt_str = ';'.join([''.join(i) for i in (permutations(letters, len(_positions_path.stem.split('-'))))])
        print(_select_adapt_str)
        _select_adapt_allowed = _select_adapt_str.split(';')
        enzyme_to_letter_conversion = {enzyme: letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')}

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
            if not _size_filter[0] < fragment_length < _size_filter[1]: continue
            if _select_adapt_str:
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
            chromosome_represented += fragment_length

        genome_represented += chromosome_represented
        perc_chromosome_represented: float = round(chromosome_represented/chr_length*100, 3)
        per_chromosome_rep_list.append(perc_chromosome_represented)

    perc_genome_represented: float = round(genome_represented/total_genome_length*100, 3)
    if not _representation_filter[0] < perc_genome_represented < _representation_filter[1]: return [None]

    return [_positions_path.stem, perc_genome_represented] + per_chromosome_rep_list


def _print_fragment_num(_positions_path: Path, _ploidy:int, _size_filter: tuple, _representation_filter: tuple, _select_adapt_str: str) -> list:
    positions_dict: dict = pickle.load(open(_positions_path, 'rb'))

    if _select_adapt_str:
        letters = [letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')]
        if len(_positions_path.stem.split('-')) > 1:
            _select_adapt_str = ';'.join([''.join(i) for i in (permutations(letters, len(_positions_path.stem.split('-'))))])
        elif len(_positions_path.stem.split('-')) == 1:
            _select_adapt_str = ';'.join([''.join(i*2) for i in (permutations(letters, len(_positions_path.stem.split('-'))))])
        _select_adapt_allowed = _select_adapt_str.split(';')
        enzyme_to_letter_conversion = {enzyme: letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')}

    fragments_per_chromosome = []
    for chr in positions_dict:
        if chr == 'metadata': continue
        fragments_passing_filter = 0
        for fragment_position, enzyme_ends in zip(positions_dict[chr]['fragment_positions'], positions_dict[chr]['enzyme_end_positions']):
            fragment_length = fragment_position[1] - fragment_position[0]
            if not _size_filter[0] < fragment_length < _size_filter[1]: continue
            if _select_adapt_str:
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
            fragments_passing_filter += 1
        fragments_per_chromosome.append(int(fragments_passing_filter*_ploidy))
    total_fragment_count: int = sum(fragments_per_chromosome)
    return [_positions_path.stem, total_fragment_count, *fragments_per_chromosome]


def _mp_print_genomic_representation(args: Namespace) -> None:


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

    headers: list = ['enzmye', total_val] + chr_list if args.chromosome_output else ['enzmye', total_val]
    print(delimiter.join(headers))
    for row in sorted([row for row in genomic_representation_list if row[0]], key=lambda x: x[1]):
        row = [str(item) for item in row] if args.chromosome_output else [str(item) for item in row][:2]
        print(delimiter.join(row))
    return None

def _print_genomic_representation2(_positions_path: Path, _ploidy:int, _size_filter: tuple, _representation_filter: tuple, _print_type: str, _select_adapt_str: str) -> list:
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

    if not _select_adapt_str:
        letters = [letter for enzyme, letter in zip(_positions_path.stem.split('-'), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')]
        _select_adapt_str = ';'.join([''.join(i) for i in (permutations(letters, 2))])
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
                if not _size_filter[0] < fragment_length < _size_filter[1]: continue
                if 'None' in enzyme_ends: continue
                processed_ends = ''.join([enzyme_to_letter_conversion[i] for i in enzyme_ends])
                if processed_ends not in _select_adapt_allowed: continue
                chromosome_represented += fragment_length

            genome_represented += chromosome_represented
            perc_chromosome_represented: float = round(chromosome_represented/chr_length*100, 3)
            per_chromosome_rep_list.append(perc_chromosome_represented)

        perc_genome_represented: float = round(genome_represented/total_genome_length*100, 3)
        if not _representation_filter[0] < perc_genome_represented < _representation_filter[1]: return [None]

        return [_positions_path.stem, perc_genome_represented] + per_chromosome_rep_list
    def _fragments_per_chromosome() -> list:
        """"""

        fragments_per_chromosome = []
        for chr in positions_dict:
            if chr == 'metadata': continue
            fragments_passing_filter = 0
            for fragment_position, enzyme_ends in zip(positions_dict[chr]['fragment_positions'], positions_dict[chr]['enzyme_end_positions']):
                fragment_length = fragment_position[1] - fragment_position[0]
                if not _size_filter[0] < fragment_length < _size_filter[1]: continue
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
