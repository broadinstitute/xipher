---
title: Metadata preparation
---
## GTF preparation
1. Obtain a GTF file for the reference genome with which your reads are aligned.
2. Convert GTF to reduced GTF using [Drop-seq](https://github.com/broadinstitute/Drop-seq/) `ReduceGTF`.
   You'll need a sequence dictionary for the reference genome.  You can create one using
   [Picard](https://github.com/broadinstitute/picard) CreateSequenceDictionary or another tool.    
   You can edit your sequence dictionary so that it only contains the X contig to make the reduced
   GTF smaller (optional).
3. Prepare the GTF for xipher.
```
xipher::prepareGtfClp(outAnnotationsPath, inReducedGtfPath)
```
## GnomAD preparation (optional)
TODO: sburger explain why to do this
1. Download the genome chrX gnomAD VCF file.  See [gnomAD downloads](https://gnomad.broadinstitute.org/downloads) for more information.
2. Run [GATK](https://gatk.broadinstitute.org/hc/en-us) VariantsToTable to convert the VCF to a table.
```
java -jar GenomeAnalysisTK.jar \
VariantsToTable \
-R your.fasta \
-V gnomad.genomes.v?.sites.chrX.vcf.bgz \
-F CHROM -F POS -F REF -F ALT -F TYPE \
-F AC -F AN -F AF \
-O gnomad_chrX_variants_table.txt
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

