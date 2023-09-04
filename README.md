
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

### Conda / Mamba

We highly recommand installing HiC-BERG through [Mamba]()

```bash

wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
mamba create -n hicberg python=3.11.4
mamba activate hicberg
mamba install -c bioconda samtools

```

TO BE COMPLETED

### Docker

TO BE COMPLETED
    
## Usage/Examples


## Full pipeline

<img src="/docs/images/Workflow.png" alt="HiC-BERG"/>

All components of the pipeline can be run at once using the hicberg pipeline command. This allows to generate a contact matrix and its reconstruction from reads in a single command.\
By default, the output is in COOL format. More detailed documentation can be found on the readthedocs website:

WEBSITE TO BE ADDED

```bash

hicberg pipeline --genome=FILE --fq-for=FILE --fq-rev=FILE [--enzyme=["DpnII", "HinfI"]]
[--rate=1.0] [--cpus=1] [--mode="full"] [--max-alignments=None] [--sensitivity="very-sensitive"] 
[--bins=2000] [--circular=""] [--mapq=35] [--rounds=1] [--magnitude=1.0] [--start-stage="fastq"] [--exit-stage=None] [--output=DIR] [--force]

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
hicberg create-folder --output=DIR [--name="folder_name"] [--force]
```

For example to create a folder named "test" on the desktop:

```bash
hicberg create-folder -o ~/Desktop/ -n test
```

The folders architecture will be the following:

```bash
├── output
│   ├── alignments
│   ├── contacts
│   ├── index
│   ├── plots
│   ├── statistics
```

### Preprocessing

After having created a folder with the previous command mentioned in **create folder**, the genome can be processed to generate fragment file __*fragment_fixed_sizes.txt*__ and the dictionary of chromosomes' sizes __*chromosome_sizes.npy*__  using the following command:

```bash
hicberg get-tables --output=DIR --genome=FILE [--bins=2000] 
``` 
For example to these files in a folder named "test" previously created on the desktop with a binning size of 2000 bp :

```bash
hicberg get-tables -o ~/Desktop/test/ -g genome.fa --bins 2000
```

The files __*fragment_fixed_sizes.txt*__ and __*chromosome_sizes.npy*__ will be generated in the folder **output/**.

### Alignment

After having created a folder with the previous command mentioned in **create folder** and performed the creation of fragment file __*fragment_fixed_sizes.txt*__ and the dictionary of chromosomes' sizes __*chromosome_sizes.npy*__ , the reads can be aligned using the following command:

```bash
hicberg alignment --genome=FILE --fq-for=FILE --fq-rev=FILE --output=DIR [--cpus=1] [--max-alignments=None] [--sensitivity="very-sensitive"] 
  [--verbosity]
```

For example to align reads in a folder named "test" previously created on the desktop with 8 threads:

```bash
hicberg alignment -g genome.fa --fq-for reads_for.fq --fq-rev rev_reads.fq -o ~/Desktop/test/ --cpus 8
```

The files __*XXX.btl2*__, __*1.sorted.bam*__ and __*2.sorted.bam*__ will be created.

### Classification

```bash
hicberg classify --output=DIR [--mapq=35] 
```

Considering the previous example, to classify the reads in a folder named "test" previously created on the desktop:

```bash
hicberg classify -o ~/Desktop/test/
```

The files created are: 

- __*group0.1.bam*__ and __*group0.2.bam*__ : bam files containing the reads of **group0** i.e. where at least one read of the pair is unaligned.
- __*group1.1.bam*__ and __*group1.2.bam*__ : bam files containing the reads of **group1** i.e. where both reads of the pair are aligned only one time.
- __*group2.1.bam*__ and __*group2.2.bam*__ : bam files containing the reads of **group2** i.e. where at least one reads of the pair are aligned more than one time.

### Pairs and matrix building

#### Build pairs

After having aligned the reads, the pairs file __*group1.pairs*__ can be built using the following command:

```bash
hicberg build-pairs --output=DIR [--recover]
```

If the flag argument *recover* is used, the pairs file will be built from the last step of the analysis e.g. after having computed the statistics and reattributed reads from **group2** bam files.

Considering the previous example, to build the matrix in a folder named "test" previously created on the desktop:

```bash
hicberg build-pairs -o ~/Desktop/test/ 
```

The file __*group1.pairs*__ will be created.

If the pairs file has to be built after reads of **group2** reassignment, the following command can be used:

```bash
hicberg build-pairs -o ~/Desktop/test/ --recover
```
Thus, the built pairs file will be  __*all_group.pairs*__.


#### Build matrix

After having aligned the reads and built the pairs file __*group1.pairs*__, the cooler matrix  __*unrescued_map.cool*__ can be built using the following command:

```bash
hicberg build-matrix  --output=DIR [--recover]
```

If the flag argument *recover* is used, the matrix file will be built from the last step of the analyis e.g. after having computed the statistics and re-attributed reads from **group2** bam files.

Considering the previous example, to build the matrix in a folder named "test" previously created on the desktop:

```bash
hicberg build-matrix -o ~/Desktop/test/ 
```

The file __*unrescued_map.cool*__ will be created.

If the cooler file has to be built after reads of **group2** re-assignament, the following command can be used:

```bash
hicberg build-matrix -o ~/Desktop/test/ --recover
```

Thus, the built matrix file will be  __*rescued_map.cool*__.



### Statistics

After having aligned the reads and built the pairs file __*group1.pairs*__, the cooler matrix  __*unrescued_map.cool*__,  the statistical laws for the reassignment of the reads from **group2** can be learnt by using the following command:

```bash
hicberg statistics --genome=FILE --output=DIR [--bins=bins_number] [--circular=""] [--rate=1.0] 
```
Considering the previous example, to get the statistical laws (with respect of ARIMA kit enzymes), without sub-sampling th restriction map and considering "chrM" as circular in a folder named "test" previously created on the desktop:

```bash
hicberg statistics -g genome.fa -e DpnII -e HinfI -c "chrM" -o ~/Desktop/test/ 
```

The statistical laws are going to be saved as:

- __*xs.npy*__ : dictionary containing the log binned genome as dictionary such as ```{chromosome: [log bins]}```
- __*uncuts.npy*__, __*loops*__, __*weirds*__ : dictionary containing the distribution of uncuts, loops and weirds as dictionary such as ```{chromosome: [distribution]}```
- __*pseudo_ps.pdf*__ : plot of the distribution of the pseudo ps, i.e. ps equivalent for trans-chromosomal cases, extracted from the reads of **group1**.
- __*coverage*__: dictionary containing the coverage of the genome as dictionary such as ```{chromosome: [coverage]}```
- __*d1d2.npy*__: np.array containing the d1d2 law as dictionary such as ```[distribution]```
- __*density_map.npy*__ : dictionary containing the density map as dictionary such as ```{chromosome_pair: [density map]}```


### Reconstruction

After having learnt the statistical laws (based on reads of **group1**), the reads from **group2** can be reassigned using the following command:

```bash
hicberg rescue --genome=FILE  --output=DIR  [--enzyme=["DpnII", "HinfI"]] [--mode="full"] [--cpus=1]
```

Considering the previous example, to reassign the reads from **group2** in a folder named "test" previously created on the desktop:

```bash

hicberg rescue -g genome.fa -e DpnII -e HinfI -o ~/Desktop/test/ 
```

The files __*group2.1.rescued.bam*__ and __*group2.2.rescued.bam*__ will be created.

### Plot

To plot all the information about the analysis, the following command can be used:

```bash
hicberg plot --genome=FILE --output=DIR  [--bins=2000]
```

Considering all the previous analysis, with 2000bp as bin size to plot all the information in a folder named "test" previously created on the desktop:

```bash
hicberg plot -g genome.fa -o ~/Desktop/test/ -b 2000
```

The plots created are:

- __*patterns_distribution_X.pdf*__ : plot of the distribution of the different patterns extracted from the reads of **group1**.
- __*coverage_X.png*__ : plot of the genome coverage extracted from the reads of **group1**.
- __*d1d2.pdf*__ : plot of the d1d2 law extracted from the reads of **group1**.
- __*density_X-Y.pdf*__ : plot of the density map extracted from the reads of **group1**.
- __*Couple_sizes_distribution.pdf*__ : plot of the distribution of the plausible couple sizes extracted from the reads of **group2**.
- __*chr_X.pdf*__ : plot of the original map and reconstructed one for each chromosome.


### Tidy folder

To tidy the folder, the following command can be used:

```bash
hicberg tidy --output=DIR
```

Considering all the previous analysis, to tidy the folder in a folder named "test" previously created on the desktop:

```bash
hicberg tidy -o ~/Desktop/test/ 
```

After tidying the folders architecture will be the following:

```bash
/home/sardine/Bureau/sample_name
├── alignments
│   ├── group0.1.bam
│   ├── group0.2.bam
│   ├── group1.1.bam
│   ├── group1.2.bam
│   ├── group2.1.bam
│   ├── group2.1.rescued.bam
│   ├── group2.2.bam
│   └── group2.2.rescued.bam
├── chunks
│   ├── chunk_for_X.bam
│   └── chunk_rev_X.bam
├── contacts
│   ├── matricies
│   │   ├── rescued_map.cool
│   │   └── unrescued_map.cool
│   └── pairs
│       ├── all_group.pairs
│       └── group1.pairs
├── fragments_fixed_sizes.txt
├── index
│   ├── index.1.bt2l
│   ├── index.2.bt2l
│   ├── index.3.bt2l
│   ├── index.4.bt2l
│   ├── index.rev.1.bt2l
│   └── index.rev.2.bt2l
├── plots
│   ├── chr_X.pdf
│   ├── Couple_sizes_distribution.pdf
│   ├── coverage_X.pdf
│   ├── patterns_distribution_X.pdf
│   ├── pseudo_ps.pdf
│   └── density_X-Y.pdf
└── statistics
    ├── chromosome_sizes.npy
    ├── coverage.npy
    ├── d1d2.npy
    ├── dist.frag.npy
    ├── loops.npy
    ├── restriction_map.npy
    ├── trans_ps.npy
    ├── uncuts.npy
    ├── weirds.npy
    └── xs.npy
```

## Evaluating the model

HiC-BERG provide a method to evaluate the infered reconstructed maps. The evaluation is based on first a split of the original uniquely mapping reads into two sets :
  
  - __*group1.X.out.bam*__ : alignment files where selected reads are complementary with the __*group1.X.in.bam*__ from the alignment files. Thus the reads are uniquely mapped (as the orginal alignment files) and used to learn the statistics for read couple inference.
  - __*group1.X_in.bam*__: alignment files where selected reads are duplicated  between all possible genomic intervals defined by the user. Thus ambiguity is introduced in the alignment of the reads.

Hence, the most plausible couple from fakely duplicated reads in __*group1.X.in.bam*__ is inferred and the corresponding Hi-C contact matrix is built and compared to the one built from the original reads in __*group1.X.bam*__ (__*unrescued_map.cool*__). The two matrices are then compared (modified bins through duplication) compared using the Pearson correlation coefficient that relates the quality of the reconstruction. The closest the coefficient is to 1, the better the reconstruction is.

The evaluation pipeline can be illustrated as follow : 

<img src="/docs/images/Evaluation.png" alt="HiC-BERG Evaluation"/>

The genomic intervals used to duplicate the reads are defined by the user through the definition of source and target intervals.The source interval is set through the parameters __*--chromosome*__ , __*--position*__ and __*--bins*__.  The target intervals are set through the parameters __*--strides*__ and eventually  __*--trans_chromosome*__ with __*--trans_position*__. 

So in a basic example considering only one chromosome and two artificial duplicated sequence, it is necessary to define a source interval corresponding to the chromosome of interest and a target interval corresponding to the duplicated sequence. The source interval is defined by the chromosome name (__*chromosome*__), the position (__*--position*__) and the width of the interval in number of bins (__*bins*__). 

Thus the source interval is defined as $[chromosome:position-bins*bin size ; chromosome:position+bins * binsize]$ and the target interval as $[chromosome:(position-bins*binsize) + stride ; chromosome:(position+bins*binsize) + stride]$. 


For example, if the source interval is __chromosome 1__, position _68000_ and strides set as __[0, 50000]__ with a bin size of __2000bp__, the source interval is defined as _chr1:68000-70000_ and the target interval is defined as _chr1:118000-120000_. 

The files __*group1.1.in.bam*__, __*group1.2.in.bam*__, __*group1.1.out.bam*__ and __*group1.2.out.bam*__ will be created.

The duplicated aligned reads should look like this :

__*group1.1.in.bam*__ :

```
NS500150:487:HNLLNBGXC:1:11101:1071:2862        0       chr1    69227   255     35M     *       0       0       ATCTGTTGTGNNGAAGGATACTCCCAGAACTCGTT     AAAAAEEEAE##EEEEEEEEEEEEEEEEEEEEEEE     AS:i:-2 XN:i:0  XM:i:2  XO:i:0  NM:i:2  MD:Z:10G0A23    YT:Z:UU   XG:i:230218
NS500150:487:HNLLNBGXC:1:11101:1071:2862        0       chr1    119227  255     35M     *       0       0       ATCTGTTGTGNNGAAGGATACTCCCAGAACTCGTT     AAAAAEEEAE##EEEEEEEEEEEEEEEEEEEEEEE     AS:i:-2 XN:i:0  XM:i:2  XO:i:0  NM:i:2  MD:Z:10G0A23    YT:Z:UU   XG:i:230218     XF:Z:Fake
NS500150:487:HNLLNBGXC:1:11101:3001:19423       16      chr1    118866  255     35M     *       0       0       GAAAAAGGATTGGTCCAATAAGTGGGAAAAAAGAT     EEAEEEAEE/EAE/EEEEEEEE/EEEEEE6AAAAA     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218
NS500150:487:HNLLNBGXC:1:11101:3001:19423       16      chr1    68866   255     35M     *       0       0       GAAAAAGGATTGGTCCAATAAGTGGGAAAAAAGAT     EEAEEEAEE/EAE/EEEEEEEE/EEEEEE6AAAAA     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218       XF:Z:Fake
NS500150:487:HNLLNBGXC:1:11101:4986:15168       16      chr1    69239   255     35M     *       0       0       GAAGGATACTCCCAGAACTCGTTACTGTCTGGACT     EEEEEEEEEEEEEEEEEEEEEEEAEEEAEEAAAAA     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218
NS500150:487:HNLLNBGXC:1:11101:4986:15168       16      chr1    119239  255     35M     *       0       0       GAAGGATACTCCCAGAACTCGTTACTGTCTGGACT     EEEEEEEEEEEEEEEEEEEEEEEAEEEAEEAAAAA     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218       XF:Z:Fake
...
```


__*group1.2.in.bam*__ :

``` 
NS500150:487:HNLLNBGXC:1:11101:1071:2862        16      chr1    103994  255     35M     *       0       0       TGCTTTTTTGGGATTGGGAATGATTTTTCCTCCTT     EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218
NS500150:487:HNLLNBGXC:1:11101:3001:19423       16      chr1    121776  255     35M     *       0       0       GGTCAAGAAATGGTTTTCACAGGCGAAATCATTGG     EEEEEEEEEEE<EEEE/EEEEEEEEEEAEEAAAAA     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218
NS500150:487:HNLLNBGXC:1:11101:4986:15168       0       chr1    86626   255     35M     *       0       0       GATCTAGGGGTACCTCCTCGGGAAACATCCAGCCC     AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE     AS:i:0  XN:i:0  XM:i:0  XO:i:0  NM:i:0  MD:Z:35 YT:Z:UU XG:i:230218
...
```

_The XF:Z:Fake signed read duplication._

In the case of trans chromosomal duplications, the user has to specify the names of trans chromosomes and the relative positions for each trans chromosome selected. The user has to provide as many position as the number of chromosome names selected. 

For example, if the source interval is __chromosome 1__, position __100000__ and strides set as __[0, 50000]__ with a bin size of __2000bp__ and the specified trans chromosomes and trans positions are respectively __[chr2, chr8]__ and __[70000, 130000]__, the source interval is defined as _chr1:100000-102000_ and the target intervals are defined as _chr1:150000-152000_, _chr2:70000-72000_ and _chr8:130000-132000_.

The stride is the number of bins between the first bin of the source interval and the first bin of the target interval. The stride can be negative or positive. If the stride is negative, the target interval is located before the source interval. If the stride is positive, the target interval is located after the source interval. The stride can be set to 0, in this case the target interval is the same as the source interval. The target interval can be located on the same chromosome as the source interval or on another chromosome. In this case, the chromosome name and the position of the first bin of the target interval must be specified. All the parameters __*--position*__, __*--strides*__, __*--trans-chromosome*__ and __*--trans-position*__ should be provided as coma separated lists.

The benchmark can be performed considering several modes. The modes are defined by the parameter __*--modes*__. The modes are defined as a list of strings separated by comas. The modes are the same as the one used for the reconstruction : 

- full
- ps_only
- cover_only
- d1d2_only
- density_only
- no_ps
- no_cover
- no_d1d2

_N.B : depending on the modes selected for the benchmark, if one of the mode is not included in the list of modes selected for the reconstruction, the reconstruction will not be performed for this mode, and the corresponding statistics will not be computed._ 



The evaluation can be run using the following command :


```bash
hicberg benchmark  --output=DIR [--chromosome] [--position] [--trans-chromosome] [--trans-position] [--stride] [--bins] [--auto] [--modes] 
```

Considering a benchmark with 4 artificially duplicated sequences set at __chr1:100000-102000 (source)__, __chr1:200000-202000 (target 1)__ __chr4:50000-52000 (target 2)__ and __chr7:300000-302000 (target 3)__, with 2000bp as bin size and considering __full and ps_only modes__ to get the performance of the reconstructions considering a folder named "test" previously created on the desktop containing the original alignment files and the unreconstructed maps, the command line is the following : 


```bash
hicberg benchmark  -o ~/Desktop/test/ -c chr1 -p 100000 -s 0,100000 -C chr4,chr7 -P 50000,30000 -m full,ps_only
```

It is also possible to let the source and target intervals being picked at random. However in such cases, the empty bins are not considered in the evaluation. The random mode is activated by setting the parameter __*--auto*__ to the number of desired artificially duplicated sequences. 

Thus, considering a benchmark with __100 artificially duplicated sequences__ , with 2000bp as bin size and considering full and ps_only modes to get the performance of the reconstructions considering a folder named "test" previously created on the desktop containing the original alignment files and the unreconstructed maps, the command line is the following : 

```bash
hicberg benchmark  -o ~/Desktop/test/ -a 100 -m full,ps_only
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

All the steps described here are handled automatically when running the ```hicberg``` pipeline. But if you want to connect the different modules manually, the intermediate input and output files can be processed using some python scripting.

## File formats

* pair files: This format is used for all intermediate files in the pipeline and is also used by hicberg build_matrix. It is a tab-separated format holding information about Hi-C pairs. It has an official specification defined by the 4D Nucleome data coordination and integration center.


* cool files: This format is used to store genomic interaction data such as Hi-C contact matrices. These file can be handled using `cooler` Python library.

* npy files: This format is used to store dictionaries containing information about genomic coordinates, binning or statistical laws. Dictionaries are stores with * chromosome * as key and * arrays* as values. Such file can be handled using `numpy` Python library. 

  * *chromosome_sizes.npy* : This file is used to store the size of each chromosome. Structure is the following : ```{chromosome: size}```
  * *xs.npy* : This file is used to store the log binned genome. Structure is the following : ```{chromosome: [log bins]}``` with log bins a list of integers.
  * *uncuts.npy* : This file is used to store the distribution of uncuts. Structure is the following : ```{chromosome: [distribution]}``` with distribution a list of integers.
  * *loops.npy* : This file is used to store the distribution of loops. Structure is the following : ```{chromosome: [distribution]}``` with distribution a list of integers.
  * *weirds.npy* : This file is used to store the distribution of weirds. Structure is the following : ```{chromosome: [distribution]}``` with distribution a list of integers.
  * *pseudo_ps.npy* : This file is used to store the distribution of pseudo ps. Structure is the following : ```{(chrom1, chrom2): [map]}``` with (chrom1, chrom2) a tuple of chromosomes where chom1 is different than chrom2 and map a float value.
  * *coverage.npy* : This file is used to store the coverage of the genome. Structure is the following : ```{chromosome: [coverage]}``` with coverage a list of integers.
  * *d1d2.npy* : This file is used to store the d1d2 law. Structure is the following : ```[distribution]``` with distribution a list of integers.
  * *density_map.npy* : This file is used to store the density map. Structure is the following : ```{(chrom1, chrom2): [density map]}``` with (chrom1, chrom2) a tuple of chromosomes density map a 2D numpy array.

* bt2l files: Thi format is used to store index of genomes performer using Bowtie2. 

* bam files: This format is used to built analyses on, by several functions of hicberg. It is a compressed standard alignement format file providing multiple informations about read alignments performer by [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml). Such files can be handled through [Samtools](http://www.htslib.org/doc/) and it's Python wrapper [PySam](https://pysam.readthedocs.io/en/latest/api.html). More details about SAM and BAM format can be found [here](https://en.wikipedia.org/wiki/SAM_(file_format)).

* fragments_fixed_sizes.txt: 

  * *chrom*: Chromosome identifier. Order should be the same as in pairs files.

  * *start*: 0-based start of fragment, in base pairs.

  * *end*: 0-based end of fragment, in base pairs.

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