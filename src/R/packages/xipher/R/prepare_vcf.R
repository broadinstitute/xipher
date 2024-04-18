# MIT License
#
# Copyright 2024 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


#' Prepare sample VCF for x inference, including gnomAD filtering if provided.
#' @param outVcfPath path to write output vcf
#' @param inVcfPath path to input vcf
#' @param isMouse logical, is the input vcf from a mouse?
#' @param gnomAdPath path to gnomAD vcf
#' @param X_contig_name name of the X contig in the reference genome. Default: `r toString(default_X_contig_name)`
#' @export
#' @import data.table
prepareVcfClp <- function(outVcfPath, inVcfPath, isMouse=FALSE, gnomAdPath=NULL,
                          X_contig_name=default_X_contig_name) {
  vcf <- data.table::fread(inVcfPath )
  if ( !is.null(gnomAdPath) ) {
  gnomad.vcf = data.table::fread( gnomAdPath )
  } else {
    gnomad.vcf = NULL
  }
  outVcf <- prepareVcf(vcf, gnomad.vcf, isMouse, X_contig_name=X_contig_name)
  write_table_helper(outVcfPath, outVcf)
}

prepareVcf<-function(vcf, gnomad.vcf = NULL, isMouse=FALSE,
                     X_contig_name=default_X_contig_name) {
  # Silence R CMD check warnings for data.table column references
  CHROM = TYPE = gt = ref = alt = het = chrom = type = NULL

  vcf[ , CHROM := gsub( "chr", "", CHROM ) ]
  vcf = vcf[ CHROM == X_contig_name & TYPE == "SNP" ]
  vcf[ , c( "CHROM", "TYPE" ) := NULL ] 
  colnames( vcf ) = do.call( "rbind", strsplit( colnames( vcf ), ".", fixed = T ) )[ , 2 ]
  colnames( vcf ) = tolower( colnames( vcf ) )
  
  # add correct gt column for mouse
  if ( isMouse ) { vcf[ , gt := paste0( ref, "|", alt ) ] }
  
  # add het column if it doesn't already exist
  if ( length( vcf$het ) == 0 ) {
    
    genotypes = vcf[ , tstrsplit( gt, "[[:punct:]]" ) ]
    vcf$het = as.numeric( genotypes[ , 1 ] != genotypes[ , 2 ] )
    
  }
  
  # ensure only heterozygous variants are present
  n.rows = nrow( vcf )
  vcf = vcf[ het == 1 ]
  vcf$het = NULL
  message( nrow( vcf ), " heterozygous variants ( ", round( nrow( vcf ) / n.rows, 2 ), " )" )
  
  # format gnomad vcf
  if ( !is.null(gnomad.vcf) ) {
    
    gnomad.vcf[ , chrom := gsub( "chr", "", chrom ) ]
    gnomad.vcf = gnomad.vcf[ chrom == X_contig_name & type == "SNP" ]
    gnomad.vcf[ , c( "chrom", "type" ) := NULL ]
    colnames( gnomad.vcf )[ -1 ] = paste0( "gnomad.", colnames( gnomad.vcf )[ -1 ] )
    vcf <- merge( vcf, gnomad.vcf, by.x = c( "pos", "ref", "alt" ), by.y = c( "pos", "gnomad.ref", "gnomad.alt" ) )
  }
  return(vcf)
}