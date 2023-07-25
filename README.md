
# HiC-BERG

A brief description of what this project does and who it's for


## Badges

Add badges from somewhere like: [shields.io](https://shields.io/)

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![AGPL License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0)

## Table of contents



## Environment and dependencies

### Environement 

Create environment by using following command : 

```
conda create -f hicberg_env.yaml
```


### Dependencies

To ensure that HiC-BERG is correctly working, Bowtie2 and Samtools have to be installed. These can be install through : 

```bash

mamba install bowtie2 -c bioconda;
mamba install samtools -c bioconda;
```


## Installation

Install my-project with pip

```bash
  pip install -e .
```

## pip

Install HiC-BERG locally by using

```bash

  pip install -e .

```

### Conda

TO BE COMPLETED

### Docker


TO BE COMPLETED



    
## Usage/Examples

```bash
import Component from 'my-project'

function App() {
  return <Component />
}
```

## Full pipeline

<img src="/docs/images/Workflow.png" alt="HiC-BERG" style="height: 100px; width:100px;"/>

All components of the pipeline can be run at once using the hicberg pipeline command. This allows to generate a contact matrix and its reconstruction from reads in a single command.\
By default, the output is in COOL format. More detailed documentation can be found on the readthedocs website:

WEBSITE TO BE ADDED

```bash

hicberg pipeline --genome=FILE --fq-for=FILE --fq-rev=FILE [--enzyme=["DpnII", "HinfI"]]
[--rate=1.0] [--cpus=1] [--mode="full"] [--max-alignments=None] [--sensitivity="very-sensitive"] 
[--bins=2000] [--circular=""] [--start-stage="fastq"] [--exit-stage=None] [--output=DIR] [--force]

```

For example, to run the pipeline using 8 threads using ARIMA Hi-C kit enzymes (DpnII & HinfI) and generate a matrix and its reconstruction in the directory out: 

```bash
hicberg pipeline -g genome.fa --fq-for reads_for.fq --fq_rev rev_reads.fq 
-e DpnII -e HinfI --cpus 8 -o out/

```

## Individual components

### Component 1

#### Sub component 1

#### Sub component 2

### Component 2

## Library

All components of the hicberg program can be used as python modules. See the documentation on reathedocs. The expected contact map format for the library is a simple COOL file, and the objects handled by the library are simple numpy arrays through Cooler. The various submodules of hicberg contain various utilities.

```python

import hicberg.io #Functions for I/O and folder management.
import hicberg.align #Functions for sequence alignment steps
import hicberg.utils #Functions for handling reads and alignment
import hicberg.statistics #Functions for extract and create statistical models
import hicberg.pipeline #Functions to run end to end Hi-C map reconstruction.

```

## Connecting the modules

## File formats

* pair files: This format is used for all intermediate files in the pipeline and is also used by hicberg build_matrix. It is a tab-separated format holding informations about Hi-C pairs. It has an official specification defined by the 4D Nucleome data coordination and integration center.


* cool files: This format is used to store genomic interaction data such as Hi-C contact matrices. These file can be handled using `cooler` Python library.

* npy files: This format is used to store dictionaries containing information about genomic coordinates, binning or statistical laws. Dictionaries are stores with * chromosome * as key and * arrays* as values. Such file can be handled using `numpy` Python library. 


* bt2l files: Thi format is used to store index of genomes performer using Bowtie2. 

* bam files: This format is used to built analyses on, by several functions of hicberg. It is a compressed standard alignement format file providing multiple informations about read alignments performer by [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). Such files can be handled through [Samtools](http://www.htslib.org/doc/) and it's Python wrapper [PySam](https://pysam.readthedocs.io/en/latest/api.html). More details about SAM and BAM format can be found [here](https://en.wikipedia.org/wiki/SAM_(file_format)).

* fragments_fixed_sizes.txt: 

  * chrom: Chromosome identifier. Order should be the same as in pairs files.

  * start: 0-based start of fragment, in base pairs.

  * end: 0-based end of fragment, in base pairs.

```
chrom start end
chr1   0       2000
chr1   2000    4000
chr1   4000    6000
...
chr1   14000   16000
...
chr2   0   2000
chr2  2000  4000
...
```


## Contributing

All contributions are welcome, in the form of bug reports, suggestions, documentation or pull requests. We use the numpy standard for docstrings when documenting functions.

The code formatting standard we use is black, with --line-length=79 to follow PEP8 recommendations. We use pytest with the pytest-doctest and pytest-pylint plugins as our testing framework. Ideally, new functions should have associated unit tests, placed in the tests folder. To test the code, you can run:

```bash

PUT CODE HERE
```


## Authors

- [@sebgra](https://www.github.com/sebgra)

## Citation
TO BE COMPLETED