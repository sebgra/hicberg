
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

<img src="/docs/images/Workflow.png" alt="HiC-BERG"/>

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

### I/O

#### Create folder

```bash
hicberg create_folder --output=DIR [--name="folder_name"] [--force]
```

For example to create a folder named "test" on the desktop:

```bash
hicberg create_folder -o ~/Desktop/ -n test
```

### Preprocessing

After having created a folder with the previous command mentionned in **create folder**, the gennome can be processed to generate fragment file __*fragment_fixed_sizes.txt*__ and the dictionary of chromosomes' sizes __*chromosome_sizes.npy*__  using the following command:

```bash
hicberg get-tables --output=DIR --genome=FILE [--bins=2000] 
``` 
For example to these files in a folder named "test" previously created on the desktop with a binning size of 2000 bp :

```bash
hicberg get-tables -o ~/Desktop/test/ -g genome.fa --bins 2000
```

### Alignment

After having created a folder with the previous command mentionned in **create folder** and performed the creation of fragment file __*fragment_fixed_sizes.txt*__ and the dictionary of chromosomes' sizes __*chromosome_sizes.npy*__ , the reads can be aligned using the following command:

```bash
hicberg alignment --genome=FILE --fq-for=FILE --fq-rev=FILE --output=DIR [--cpus=1] [--max-alignments=None] [--sensitivity="very-sensitive"] 
  [--verbosity]
```

For example to align reads in a folder named "test" previously created on the desktop with 8 threads:

```bash
hicberg alignment -g genome.fa -f reads_for.fq -r rev_reads.fq -o ~/Desktop/test/ --cpus 8
```

### Pairs and matrix building

#### Build pairs

After having aligned the reads, the pairs file __*group1.pairs*__ can be built using the following command:

```bash
hicberg build-pairs --output=DIR [--recover]
```

If the flag argument *recover* is used, the pairs file will be built from the last step of the analyis e.g. after having computed the statistics and reattributed reads from **group2** bam files.

Considering the previous example, to build the matrix in a folder named "test" previously created on the desktop:

```bash
hicberg build-pairs -o ~/Desktop/test/ 
```

If the pairs file has to be built after reads of **group2** reassignament, the following command can be used:

```bash
hicberg build-pairs -o ~/Desktop/test/ --recover
```
Thus, the built pairs file will be  __*all_group.pairs*__.


#### Build matrix

After having aligned the reads and built the pairs file __*group1.pairs*__, the cooler matrix  __*unrescued_map.cool*__ can be built using the following command:

```bash
hicberg build-matrix  --output=DIR [--recover]
```

If the flag argument *recover* is used, the matrix file will be built from the last step of the analyis e.g. after having computed the statistics and reattributed reads from **group2** bam files.

Considering the previous example, to build the matrix in a folder named "test" previously created on the desktop:

```bash
hicberg build-pamatrixirs -o ~/Desktop/test/ 
```

If the cooler file has to be built after reads of **group2** reassignament, the following command can be used:

```bash
hicberg build-matrix -o ~/Desktop/test/ --recover
```

Thus, the built matrix file will be  __*rescued_map.cool*__.


### Classification

```bash
hicberg classify --output=DIR
```

Considering the previous example, to classify the reads in a folder named "test" previously created on the desktop:

```bash
hicberg classify -o ~/Desktop/test/
```

### Statistics

After having aligned the reads and built the pairs file __*group1.pairs*__, the cooler matrix  __*unrescued_map.cool*__,  the statistical laws for the reassignment of the reads from **group2** can be learnt by using the following command:

```bash
hicberg statistics --genome=FILE --output=DIR [--bins=bins_number] [--circular=""] [--rate=1.0] 
```
Considering the previous example, to get the statistical laws (with respect of ARIMA kit enzymes), without subsampling th restriction map and considering "chrM" as circular in a folder named "test" previously created on the desktop:

```bash
hicberg statistics -g genome.fa -e DpnII -e HinfI -c "chrM" -o ~/Desktop/test/ 
```

The statistical laws are going to be saved as:

- __*xs.npy*__ : dictionary containing the log binned genome as dictionary such as ```{chromosome: [log bins]}```
- __*uncuts.npy*__, __*loops*__, __*weirds*__ : dictionary containing the distribution of uncuts, loops and weirds as dictionary such as ```{chromosome: [distribution]}```
- __*coverage*__: dictionary containing the coverage of the genome as dictionary such as ```{chromosome: [coverage]}```
- __*d1d2.npy*__: np.array containing the d1d2 law as dictionary such as ```[distribution]```


### Reconstruction

After having learnt the statistical laws (based on reads of **group1**), the reads from **group2** can be reassigned using the following command:

```bash
hicberg rescue --genome=FILE  --output=DIR  [--enzyme=["DpnII", "HinfI"]] [--mode="full"] [--cpus=1]
```

Considering the previous example, to reassign the reads from **group2** in a folder named "test" previously created on the desktop:

```bash

hicberg rescue -g genome.fa -e DpnII -e HinfI -o ~/Desktop/test/ 
```

### Plot

```bash

```

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