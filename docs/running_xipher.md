---
title: Running xipher
---
## Genotyping
TODO: sburger, alecw

#### Prepare VCF
Load the VCF created in Genotyping step, optionally filter based on gnomAD, and prepare it for xipher.
```
xipher::prepareVcfClp(outVcfPath, inVcfPath, gnomAdPath)
```

## Prepare digital allele counts (DAC)
1. Run [Drop-seq](https://github.com/broadinstitute/Drop-seq/) `GatherDigitalAlleleCounts` on the BAM file for each reaction.
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
- `VCF` is the file produced as part of genotyping step. TODO: sburger, is this right?
- `INPUT`: the BAM file for the reaction.
- `OUTPUT`: the output file.  It must have extension .dac.txt.gz,
  and the filename before the extension will be used to prefix the cell barcodes to disambiguate across reactions.
2. Aggregate the DAC files for all reactions into a single file and prepare for xipher.
```
xipher::prepareDacClp(outDacPath, dacPaths, annotationsPath, vcfPath)
```
- `dacPaths` is a vector of paths to the DAC files created above.
- `annotationsPath` is the path to the annotations file created in GTF preparation step.
- `vcfPath` is the path to the VCF file created in Prepare VCF step.

## Phase X SNPs
```
xipher::phaseClp(outPhasePath,preparedDacPath)
```
- preparedDacPath is the file created by `xipher::prepareDacClp()`

## Call X-inactivation
```
xipher::callActiveXClp(outErrorRatePath, outXCallsPath, preparedDacPath, phasePath)
```
- `preparedDacPath` is the file created by `xipher::prepareDacClp()
- `phasePath` is the file created by `xipher::phaseClp()`

