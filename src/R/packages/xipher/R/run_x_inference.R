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

# setup -----------------------------------------------------------------------------------------------------------


# currently all human analyses go through gnomad phasing
if ( !yaml$is.mouse ) {
  
  message( "using gnomad filtering" )
  filter.using.gnomad = TRUE
  
} else { 
  
  message( "not using gnomad filtering ( non-human )" )
  filter.using.gnomad = FALSE 
  
}



# run initial x phasing -------------------------------------------------------------------------------------------

# load donor dac
dac = load.raw.dac( yaml$dac.dir )
dac[ , chr := gsub( "chr", "", chr ) ]
dac = dac[ chr == "X" ]

# load donor vcf
vcf = fread( yaml$variant.table.path )
vcf[ , CHROM := gsub( "chr", "", CHROM ) ]
vcf = vcf[ CHROM == "X" & TYPE == "SNP" ]
vcf[ , c( "CHROM", "TYPE" ) := NULL ] 
colnames( vcf ) = do.call( "rbind", strsplit( colnames( vcf ), ".", fixed = T ) )[ , 2 ]
colnames( vcf ) = tolower( colnames( vcf ) )

# add correct gt column for mouse
if ( yaml$is.mouse ) { vcf[ , gt := paste0( ref, "|", alt ) ] }

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

# merge dac and vcf
dac = merge( dac, vcf, by = "pos" )

# format gnomad vcf
if ( filter.using.gnomad ) {
  
  gnomad.vcf = readRDS( yaml$gnomad.path )
  gnomad.vcf[ , chrom := gsub( "chr", "", chrom ) ]
  gnomad.vcf = gnomad.vcf[ chrom == "X" & type == "SNP" ]
  gnomad.vcf[ , c( "chrom", "type" ) := NULL ]
  colnames( gnomad.vcf )[ -1 ] = paste0( "gnomad.", colnames( gnomad.vcf )[ -1 ] )
  
  # merge dac and gnomad.vcf ( triallelic sites may be lost, but this shouldn't be an issue )
  dac = merge( dac, gnomad.vcf, by.x = c( "pos", "ref", "alt" ), by.y = c( "pos", "gnomad.ref", "gnomad.alt" ) )
  
}


# add position-gene annotations ( exonic, intronic ) to dac
gtf = readRDS( yaml$gtf.path )
unique.pos = unique( dac[ , .( pos, gene ) ] )
pos.gene.annotations = data.table()
for ( i in 1:nrow( unique.pos ) ) {
  
  pos.use = unique.pos[ i, pos ]
  gene.use = unique.pos[ i, gene ]
  
  annotation = gtf[ gene == gene.use & pos.use > start & pos.use < end, annotation ]
  annotation = unique( annotation )
  annotation = sort( annotation )
  annotation = paste0( annotation, collapse = "," )
  
  pos.gene.annotations = rbind( pos.gene.annotations, data.table( pos = pos.use, gene = gene.use, annotation ) )
  
}
shared.cols = intersect( colnames( dac ), colnames( pos.gene.annotations ) )
dac = merge( dac, pos.gene.annotations, by = shared.cols )

# HERE: skip this section to include exonic positions
if ( yaml$intronic.only ) {
  
  # keep only intronic positions
  n.rows = nrow( dac )
  message( "keeping only intronic reads" )
  dac = dac[ annotation == "intron" ]
  n.lost = n.rows - nrow( dac )
  message( "removed ", n.lost, " ( ", round( n.lost / n.rows, 2 ), " ) rows from dac" )
  
}

# add haplotype ( allele ) counts
dac = cbind( dac, data.table( t( apply( dac, 1, assign.bases.to.haplotypes ) ) ) )

saveRDS( dac, "unphased_x_dac.RDS" )

start = Sys.time()
flip.record = phase.x.variants( dac, min.cells = yaml$min.cells, max.xi.proportion = yaml$max.xi.proportion, 
                                time.limit = 2, max.gen.units = yaml$max.gen.units )
end = Sys.time()
end - start

saveRDS( flip.record, "flip_record.RDS" )



# get error rates and x calls -------------------------------------------------------------------------------------

# print libraries
print( table( sapply( unique( dac$cell ), get.cell.barcode.prefix ) ) )

phased.positions = unlist( strsplit( flip.record[ which.max( pct.total ), contains ], ",", fixed = TRUE ) )
message( length( phased.positions ), " phased postions" )
flip.positions = flip.record[ pos %in% phased.positions & flip, pos ]
message( length( flip.positions ), " flip postions ( ", round( length( flip.positions ) / length( phased.positions ), 2 ), " )" )

if ( length( phased.positions ) < 1 ) { 
  
  message( "insufficient initial phasing ( quitting )" )
  quit()
  
}

x.calls = call.active.x( dac, phased.positions, flip.positions )
saveRDS( x.calls, "x_calls_initial.RDS" )


position.stats = generate.error.rates.and.update.phasing( dac, phased.positions, flip.positions, min.cell = yaml$min.cells.final )

phased.positions = position.stats[ ( ( error.rate <= yaml$max.error.rate ) | ( error.rate >= ( 1 - yaml$max.error.rate ) ) ) & total >= yaml$min.cells.final, pos ]
message( length( phased.positions ), " phased postions" )
flip.positions = position.stats[ pos %in% phased.positions & error.rate > 0.5, pos ]
message( length( flip.positions ), " flip postions ( ", round( length( flip.positions ) / length( phased.positions ), 2 ), " )" )
position.uncertainties = position.stats[ , .( pos, gene, uncertainty ) ]

x.calls = call.active.x( dac, phased.positions, flip.positions, position.uncertainties = position.uncertainties,
                         calling.purity.min = yaml$calling.purity.min, likelihood.min = yaml$likelihood.min )

saveRDS( position.stats, "error_rates.RDS" )
saveRDS( x.calls, "x_calls.RDS" )

message( "done!" )

