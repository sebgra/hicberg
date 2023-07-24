
# Project Title

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

All components of the pipeline can be run at once using the hicberg pipeline command. This allows to generate a contact matrix and its reconstruction from reads in a single command.\
By default, the output is in COOL format. More detailed documentation can be found on the readthedocs website:

WEBSITE TO BE ADDED

```bash

hicberg pipeline --genome=FILE --fq-for=FILE --fq-rev=FILE [--enzyme=["DpnII", "HinfI"]]
[--rate=1.0] [--cpus=1] [--mode="full"] [--max-alignments=None] [--sensitivity="very-sensitive"] 
[--bins=2000] [--circular=""] [--start-stage="fastq"] [--exit-stage=None] [--output=DIR] [--force]

```

For example, to run the pipeline using 8 threads using ARIMA Hi-C kit enzymes (DpnII & HinfI) and generate a matrix and its reconstruction in the directory out : 

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

* pair files


* cool files





## Contributing

Contributions are always welcome!

See `contributing.md` for ways to get started.

Please adhere to this project's `code of conduct`.


## Authors

- [@sebgra](https://www.github.com/sebgra)

## Citation