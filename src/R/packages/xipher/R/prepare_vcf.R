#' load donor vcf, prepare for x inference, including gnomAD filtering if provided.
#' @param outVcfPath path to write output vcf
#' @param inVcfPath path to input vcf
#' @param isMouse logical, is the input vcf from a mouse?
#' @param gnomAdPath path to gnomAD vcf
#' @export
#' @import data.table, R.utils, utils
prepareVcfClp <- function(outVcfPath, inVcfPath, isMouse=FALSE, gnomAdPath=NULL) {
  vcf <- data.table::fread(inVcfPath )
  if ( !is.null(gnomAdPath) ) {
  gnomad.vcf = data.table::fread( gnomAdPath )
  } else {
    gnomad.vcf = NULL
  }
  outVcf <- prepareVcf(vcf, gnomad.vcf, isMouse)
  outconn <- open_conn(outVcfPath, "w")
  utils::write.table(outVcf, outconn, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  close(outconn)
}

prepareVcf<-function(vcf, gnomad.vcf = NULL, isMouse=FALSE) {
  
  vcf[ , CHROM := gsub( "chr", "", CHROM ) ]
  vcf = vcf[ CHROM == "X" & TYPE == "SNP" ]
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
    gnomad.vcf = gnomad.vcf[ chrom == "X" & type == "SNP" ]
    gnomad.vcf[ , c( "chrom", "type" ) := NULL ]
    colnames( gnomad.vcf )[ -1 ] = paste0( "gnomad.", colnames( gnomad.vcf )[ -1 ] )
    vcf <- merge( vcf, gnomad.vcf, by.x = c( "pos", "ref", "alt" ), by.y = c( "pos", "gnomad.ref", "gnomad.alt" ) )
  }
  return(vcf)
}