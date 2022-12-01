# py-simRAD
A parallel implementation of simRAD (Lepais & Weir, 2014) in Python using BioPython (Cock et al., 2009) for simulated prediction of loci expected in RADseq.

## Table of contents
* [Installation](#installation)
* [Usage](#usage)
  * [1. Catalysis: Generating fragment positions](#1-catalysis-generating-fragment-positions)
  * [2A. Export FASTA](#2a-export-fasta)
  * [2B. Export GFF](#2b-export-gff)
  * [3. Summary genomic representation](#3-summary-genomic-representation)
  * [4. Batched functionality](#4-batched-functionality)
* [Example](#example)
* [References](#references)

## Installation
The installation is quite easy. The only requirement is BioPython which can be installed using the [BioPython official documentation](https://biopython.org/wiki/Packages).
They generally recommend using the `pip` package manager:
```
pip install biopython
```

For use with `conda`, BioPython can be installed from conda packages:
```
conda install -c conda-forge biopython
```

## Usage
### 1. Catalysis: Generating fragment positions
Before visualization, exporting, or any other summarization of the restriction fragment data, restriction fragment positions must first be generated. 
```
python src/py_simRAD.py catalyze \
    "raw-data/Trichoderma_atroviride_genomic.fna" \
    --enzymes "HhaI;HindIII;NotI"
```

This will generate a directory containing files (`.pos`) of restriction fragment position data for each enzyme and enzyme pair. These files are meant to be used as input for each of the export or summary options.

### 2A. Export FASTA
After generation of fragment positions, these positions can be further filtered and exported into FASTA format.
```
python src/py_simRAD.py export \
    "raw-data/Trichoderma_atroviride_genomic.fna" \
    "output-dir" \
    "raw-data/.Trichoderma_atroviride_genomic/HhaI-HindIII.pos" \
    --min 300 \
    --max 600 \
    --type 'fasta'
```

The result is a FASTA file containing sequences of restriction fragments from the given restriction combination used on the given genome.
```
HEADER FORMAT     >{ID} {DESCRIPTION} {ENZYME_COMBINATION} {POSITIONS}
EXAMPLE           >CP084935.1 Trichoderma atroviride strain P1 chromosome 1 HhaI-HindIII 0-403
```

### 2B. Export GFF
Again, this command is meant to be performed after generation of fragment positions. These features can be filtered according to fragment lengths before export into GFF format and can be used with IGV.
```
python src/py_simRAD.py export \
    "raw-data/Trichoderma_atroviride_genomic.fna" \
    "output-dir" \
    "raw-data/.Trichoderma_atroviride_genomic/HhaI-HindIII.pos" \
    --min 300 \
    --max 600 \
    --type 'gff'
```

```
FORMAT            {chromosome}	{py-simRAD version}	restriction_fragment	{start} {end}	.	+	.
EXAMPLE           CP084935.1	py-simRADv4.1.2	restriction_fragment	0	403	.	+	.
```
### 3. Summary genomic representation
Use this command to print out a tab-delimited report of genomic representation.
```
python src/py_simRAD.py summary \
    "raw-data/.Trichoderma_atroviride_genomic/HhaI-HindIII.pos" \
    --min 300 \
    --max 600 \
```

```
enzmye  total repr (%)  CP084935.1      CP084936.1      CP084937.1      CP084938.1      CP084939.1      CP084940.1      CP084941.1
HhaI-HindIII    37.26   37.367  36.282  36.528  36.439  39.51   37.94   37.962
```

The default behaviour is to print to console with human-readable tab delimiting, but this behaviour may be changed. As delimited output, the console printout of this function can be easily redirected to a file:
```
python src/py_simRAD.py summary \
    "raw-data/.Trichoderma_atroviride_genomic/HhaI-HindIII.pos" \
    --min 300 \
    --max 600 \
    --delimiter ',' \
    > target_file.csv
```

### 4. Batched functionality
Each of the export and summary functions have been written so that wildcard (`*`) file input is allowed.

#### Summary genomic representation
The following program takes an input of every restriction enzyme combination generated (`*`) and filters printout according to fragment size and percent genomic representation.
```
python src/py_simRAD.py summary \
    "raw-data/.Trichoderma_atroviride_genomic/*" \
    --min 300 \
    --max 600 \
    --delimiter ',' \
    --min_rep 5 \
    --max_rep 15 \
    > target_file.csv
```

## Example


## References

> Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163


> Lepais, O., & Weir, J. T. (2014). SimRAD: an R package for simulation-based prediction of the number of loci expected in RADseq and similar genotyping by sequencing approaches. Molecular ecology resources, 14(6), 1314–1321. https://doi.org/10.1111/1755-0998.12273