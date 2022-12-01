# py-simRAD
A parallel implementation of simRAD (Lepais & Weir, 2014) in Python using BioPython (Cock et al., 2009) for simulated prediction of loci expected in RADseq.

## Table of contents
* [Installation](#installation)
* [Usage](#usage)
  * [1. Catalysis: Generating fragment positions](#1-catalysis-generating-fragment-positions)
* [Example](#example)
* [References](#references)

## Installation
The installation is quite easy.

## Usage
### 1. Catalysis: Generating fragment positions
Before visualization, exporting, or any other summarization of restriction fragment data, restriction fragment positions must first be generated.
```
python src/py_simRAD.py catalyze \
    "raw-data/Trichoderma_atroviride_genomic.fna" \
    --enzymes "HhaI;HindIII;NotI"
```

## Example

## References

> Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163


> Lepais, O., & Weir, J. T. (2014). SimRAD: an R package for simulation-based prediction of the number of loci expected in RADseq and similar genotyping by sequencing approaches. Molecular ecology resources, 14(6), 1314–1321. https://doi.org/10.1111/1755-0998.12273