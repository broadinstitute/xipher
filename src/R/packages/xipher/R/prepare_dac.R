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

#' Join DACs and VCF, annotate with GTF, and prepare for phasing
#' @param outDacPath path to write output DAC
#' @param dacPaths paths to DAC files, one for each reaction.  Filenames are used to disambiguate CBCs
#' @param annotationsPath reduced GTF with columns chr,start,end,gene,annotation
#' @param vcfPath X-only VCF created with prepareVcfClp(), with columns pos	ref,alt,ps,gt,ad,dp,gq,pl (and if gnomAD:	gnomad.ac,gnomad.an,gnomad.af,gnomad.maf)
#' @param intronicOnly logical, keep only intronic reads?
#' @param X_contig_name name of the X contig
#' @export 
#' @import data.table
prepareDacClp <- function(outDacPath,
                     dacPaths, 
                     annotationsPath, 
                     vcfPath, 
                     intronicOnly=TRUE,
                     X_contig_name=default_X_contig_name) {
  dac <- load_and_merge_dacs(dacPaths)
  vcf <- data.table::fread( vcfPath )
  gtf <- data.table::fread( annotationsPath )
  dac[ , chr := gsub( "chr", "", chr ) ]
  # TODO: Steve: What does this line do?
  dac = dac[ chr == X_contig_name ]
  dac = merge( dac, vcf, by = "pos" )
  unique.pos = unique( dac[ , .( pos, gene ) ] )
  pos.gene.annotations = data.table::data.table()
  for ( i in 1:nrow( unique.pos ) ) {
    
    pos.use = unique.pos[ i, pos ]
    gene.use = unique.pos[ i, gene ]
    
    annotation = gtf[ gene == gene.use & pos.use > start & pos.use < end, annotation ]
    annotation = unique( annotation )
    annotation = sort( annotation )
    annotation = paste0( annotation, collapse = "," )
    
    pos.gene.annotations = rbind( pos.gene.annotations, data.table::data.table( pos = pos.use, gene = gene.use, annotation ) )
    
  }
  shared.cols = intersect( colnames( dac ), colnames( pos.gene.annotations ) )
  dac = merge( dac, pos.gene.annotations, by = shared.cols )
  
  # HERE: skip this section to include exonic positions
  if ( intronicOnly ) {
    
    # keep only intronic positions
    n.rows = nrow( dac )
    message( "keeping only intronic reads" )
    dac = dac[ annotation == "intron" ]
    n.lost = n.rows - nrow( dac )
    message( "removed ", n.lost, " ( ", round( n.lost / n.rows, 2 ), " ) rows from dac" )
    
  }
  
  # add haplotype ( allele ) counts
  dac = cbind( dac, data.table::data.table( t( apply( dac, 1, assign.bases.to.haplotypes ) ) ) )
  
  write_table_helper(outDacPath, dac)
}