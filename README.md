# xipher
Tools for calling X-inactivation in scRNA-seq data

## Introduction

TO COME

## Installation

```remotes::install_github("broadinstitute/xipher", subdir="src/R/packages/xipher")```

The above installs the non-necessarily stable main branch.  
We will soon create tagged branches for stable releases.

## Usage

### Metadata preparation

#### GTF preparation
1. Obtain a GTF file for the reference genome with which your reads are aligned.
2. Convert GTF to reduced GTF using [Drop-seq](https://github.com/broadinstitute/Drop-seq/) ReduceGTF.
You'll need a sequence dictionary for the reference genome.  You can create one using
[Picard](https://github.com/broadinstitute/picard)'s CreateSequenceDictionary or another tool.    
You can edit your sequence dictionary so that it only contains the X contig to make the reduced 
GTF smaller.
3. Run `xipher::prepareGtfClp(outAnnotationsPath, inReducedGtfPath)` to prepare the GTF for xipher.

#### GnomAD preparation (optional)
TO COME

### Running xipher

#### Genotyping
TO COME

