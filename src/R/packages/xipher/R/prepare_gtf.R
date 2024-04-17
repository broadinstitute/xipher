#' Prepare GTF gene annotations for xipher
#' @param outAnnotationsPath path to write gene annotations in format used by prepareDacClp()
#' @param inReducedGtfPath path to reduced GTF file produced by ReduceGtf (see https://github.com/broadinstitute/Drop-seq/).
#'                         Note that you can speed ReduceGtf and make it's output smaller by providing a sequence dictionary
#'                         that contains only the X contig.
#' @param X_contig_name name of the X contig. Default: `r toString(default_X_contig_name)`
#' @export
#' @import data.table
prepareGtfClp<-function(outAnnotationsPath,
                        inReducedGtfPath,
                        X_contig_name=default_X_contig_name) {
  gtf = data.table::fread( inReducedGtfPath )
  
  gtf = unique( gtf[ , .( chr, start, end, gene = gene_name, annotation = annotationType ) ] )
  gtf = gtf[ chr == X_contig_name ]
  gtf = gtf[ annotation %in% c( "exon", "intron" ) ]
  write_table_helper(outAnnotationsPath, gtf)
}