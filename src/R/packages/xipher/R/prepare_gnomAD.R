#' Prepare gnomAD variants table for xipher::prepareVcfClp.
#' @param outGnomADPath path to write gnomAD variants table.  You can then pass this to prepareVcfClp() with gnomAdPath=
#' @param inGnomADVariantsTablePath path to gnomAD variants table produced by GATK VariantsToTable
#' @param maf.min minimum allele frequency
#' @export
#' @import data.table
prepareGnomAdClp<-function(outGnomADPath,
                          inGnomADVariantsTablePath,
                          maf.min = 0.001) {
  # Silence R CMD check warnings for data.table column references
  chrom = type = af = maf = NULL
  gnomad <- data.table::fread(inGnomADVariantsTablePath)

  colnames( gnomad ) = tolower( colnames( gnomad ) )
  
  gnomad[ , maf := min( af, 1 - af ), by = 1:nrow( gnomad ) ]
  gnomad = gnomad[ maf >= maf.min ]
  
  gnomad = gnomad[ type == "SNP" ]
  
  gnomad = gnomad[ chrom == "X" | chrom == "chrX" ]
  
  setcolorder( gnomad, c( "chrom", "pos", "ref", "alt", "type", "ac", "an", "af", "maf" ) )
  write_table_helper(outGnomADPath, gnomad)
}