# xipher
Tools for calling X-inactivation in scRNA-seq data
TODO: sburger

## Introduction

TODO: sburger

## Installation

```
remotes::install_github("broadinstitute/xipher", subdir="src/R/packages/xipher")
```

The above installs the not-necessarily-stable main branch.  
We will soon create tagged branches for stable releases.
TODO: alecw

## Other tools you will need
- [Drop-seq](https://github.com/broadinstitute/Drop-seq/) for ReduceGTF and GatherDigitalAlleleCounts.
- [GATK](https://gatk.broadinstitute.org/hc/en-us) for VariantsToTable.
- [Picard](https://github.com/broadinstitute/picard) for CreateSequenceDictionary (optional). 

## Usage

### Metadata preparation

#### GTF preparation
1. Obtain a GTF file for the reference genome with which your reads are aligned.
2. Convert GTF to reduced GTF using ReduceGTF.
You'll need a sequence dictionary for the reference genome.  You can create one using
Picard's CreateSequenceDictionary or another tool.    
You can edit your sequence dictionary so that it only contains the X contig to make the reduced 
GTF smaller.
3. Prepare the GTF for xipher.
```
xipher::prepareGtfClp(outAnnotationsPath, inReducedGtfPath)
```
#### GnomAD preparation (optional)
TODO: sburger explain why to do this
1. Download the genome chrX gnomAD VCF file.  See [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) for more information.
2. Run GATK VariantsToTable to convert the VCF to a table.
```
java -jar GenomeAnalysisTK.jar /
-R your.fasta /
-T VariantsToTable /
-V gnomad.genomes.v?.sites.chrX.vcf.bgz /
-F CHROM -F POS -F REF -F ALT -F TYPE /
-F AC -F AN -F AF /
-o gnomad_chrX_variants_table.txt
```
- `your.fasta` is the reference genome fasta file you used to align your reads.
- `gnomad.genomes.v?.sites.chrX.vcf.bgz` is the gnomAD VCF file you downloaded.
- `gnomad_chrX_variants_table.txt` is the output file.
3. Prepare the gnomAD table for xipher.
```
xipher::prepareGnomAdClp(outGnomADPath, inGnomADVariantsTablePath)
```
- `outGnomADPath` is the path to the output file.
- `inGnomADVariantsTablePath` is the path to the table created by GATK VariantsToTable.

### Running xipher

#### Genotyping
TODO: sburger, alecw

#### Prepare VCF
Load the VCF created in Genotyping step, optionally filter based on gnomAD, and prepare it for xipher.
```
xipher::prepareVcfClp(outVcfPath, inVcfPath, gnomAdPath)
```

#### Prepare digital allele counts (DAC)
1. Run `GatherDigitalAlleleCounts` on the BAM file for each reaction.
```
GatherDigitalAlleleCounts INPUT=your.bam \
    OUTPUT=your.dac.txt.gz \
    VCF=your.vcf \
    CELL_BC_FILE=your.cell_barcodes.txt
	HET_SNPS_ONLY=true \
	SINGLE_VARIANT_READS=true \
	MULTI_GENES_PER_READ=false \
	LOCUS_FUNCTION_LIST=INTRONIC \
	POLYMORPHIC_SNPS_ONLY=false \
	GQ_THRESHOLD=20
```
- `CELL_BC_FILE` is a file with the cell barcodes for the reaction.
- `VCF` is TODO: sburger
- `INPUT`: the BAM file for the reaction
- `OUTPUT`: the output file.  It must have extension .dac.txt.gz,
  and the filename before the extension will be used to prefix the cell barcodes to disambiguate across reactions.
2. Aggregate the DAC files for all reactions into a single file and prepare for xipher.
`xipher::prepareDacClp(outDacPath, dacPaths, annotationsPath, vcfPath)`
- `dacPaths` is a vector of paths to the DAC files created above
- `annotationsPath` is the path to the annotations file created in GTF preparation step
- `vcfPath` is the path to the VCF file created in Prepare VCF step

#### Phase X SNPs
`xipher::phaseClp(outPhasePath,preparedDacPath)`
- preparedDacPath is the file created by `xipher::prepareDacClp()`

#### Call X-inactivation
`xipher::callActiveXClp(outErrorRatePath, outXCallsPath, preparedDacPath, phasePath)`
- `preparedDacPath` is the file created by `xipher::prepareDacClp()
- `phasePath` is the file created by `xipher::phaseClp()`

## Understanding the output
TODO: sburger
