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

# general ---------------------------------------------------------------------------------------------------------

#'
#'
get.cell.barcode.prefix = function( cell.name ) {
  
  name.fragments = unlist( strsplit( cell.name, "_", fixed = TRUE ) )
  barcode.length = nchar( name.fragments[ length( name.fragments ) ] )
  barcode.prefix = substr( cell.name, 1, nchar( cell.name ) - ( barcode.length + 1 ) )
  
  return( barcode.prefix )
  
}


#'
#'
get.cell.barcode.suffix = function( cell.name ) {
  
  name.fragments = unlist( strsplit( cell.name, "_", fixed = TRUE ) )
  barcode.length = nchar( name.fragments[ length( name.fragments ) ] )
  barcode.suffix = substr( cell.name, nchar( cell.name ) - ( barcode.length - 1 ), nchar( cell.name ) )
  
  return( barcode.suffix )
  
}


#'
#'
pivot.wider = function( data.table, rows.from, cols.from, values.from ) {
  
  out = pivot_wider( data.table, id_cols = rows.from, names_from = all_of( cols.from ), values_from = values.from )
  
  row.names = out[[ 1 ]]
  out[[ 1 ]] = NULL
  out = as.matrix( out )
  rownames( out ) = row.names
  
  return( out )
  
}


#'
#'
get.corresponding.labels = function( labels, from, to ) {
  
  from = as.character( from )
  labels = as.character( labels )
  names( to ) = from
  new.labels = to[ labels ]
  names( new.labels ) = NULL
  
  return( new.labels )
  
}


#'
#'
get.corresponding.cell.names = function( cell.names.a, cell.names.b, overlap.threshold = 10 ^ 5, plot = TRUE ) {
  
  # overlap threshold should be 10 ^ 5 for 10X V3 libraries, 10 ^ 4 for 10X V2, 10 ^ 3 for drop-seq
  
  set.a = data.table::data.table( library.a = get.cell.barcode.prefix( cell.names.a ), cbc = get.cell.barcode.suffix( cell.names.a ) )
  set.a.table = set.a[ , table( library.a ) ]
  range( set.a.table )
  if ( any( set.a.table < 1000 ) ) { message( "WARNING: libraries with fewer than 1000 cells present in set.a" ) }
  
  set.b = data.table::data.table( library.b = get.cell.barcode.prefix( cell.names.b ), cbc = get.cell.barcode.suffix( cell.names.b ) )
  set.b.table = set.b[ , table( library.b ) ]
  range( set.b.table )
  if ( any( set.b.table < 1000 ) ) { message( "WARNING: libraries with fewer than 1000 cells present in set.b" ) }
  
  cbc.length = median( c( nchar( set.a$cbc ), nchar( set.b$cbc ) ) )
  
  merged.sets = merge( set.a, set.b, by = "cbc", all = TRUE )
  
  overlap.table.counts = merged.sets[ , table( library.a, library.b, useNA = "a" ) ]
  overlap.table.proportions = overlap.table.counts
  for ( i in 1:nrow( overlap.table.proportions ) ) {
    
    for ( j in 1:ncol( overlap.table.proportions ) ) {
      
      set.a.counts = set.a[ library.a == rownames( overlap.table.proportions )[ i ], .N ]
      set.b.counts = set.b[ library.b == colnames( overlap.table.proportions )[ j ], .N ]
      
      expected.overlap = set.a.counts * set.b.counts * ( 0.25 ^ cbc.length )
      
      overlap.table.proportions[ i , j ] = overlap.table.counts[ i, j ] / expected.overlap
      
    }
    
  }
  
  overlap.table.proportions = overlap.table.proportions[ !is.na( rownames( overlap.table.proportions ) ), 
                                                         !is.na( colnames( overlap.table.proportions ) ), drop = FALSE ]
  # overlap.table.proportions = round( overlap.table.proportions, 2 )
  
  if ( plot ) { 
    
    x.min = -1
    
    hist( log10( unlist( overlap.table.proportions ) ), 
          ylim = c( 0, min( dim( overlap.table.proportions ) ) ), xlim = c( x.min, 7 ),
          main = "", col = "grey90", border = "white",
          breaks = seq( x.min, 7, 0.1 ), xlab = "log10( observed / expected )" )
    
    hist( log10( unlist( overlap.table.proportions[ overlap.table.proportions > overlap.threshold ] ) ), 
          col = "mediumseagreen", border = "white",
          breaks = seq( x.min, 7, 0.1 ), add = TRUE )
    
    abline( v = log10( overlap.threshold ), lty = 2 )
    
  }
  
  binary.overlap.table = overlap.table.proportions > overlap.threshold
  
  
  set.a.missing = rownames( binary.overlap.table )[ !apply( binary.overlap.table, 1, any, na.rm = TRUE ) ]
  set.a.missing = set.a.missing[ !is.na( set.a.missing ) ]
  set.b.missing = colnames( binary.overlap.table )[ !apply( binary.overlap.table, 2, any, na.rm = TRUE ) ]
  set.b.missing = set.b.missing[ !is.na( set.b.missing ) ]
  
  if ( length( set.a.missing > 0 ) ) {
    
    message( "cell.names.a libraries with no match in cell.names.b" )
    print( set.a.missing )
    
  }
  
  if ( length( set.a.missing > 0 ) ) {
    
    message( "cell.names.b libraries with no match in cell.names.a" )
    print( set.b.missing )
    
  }
  
  
  array.indices = which( binary.overlap.table, arr.ind = TRUE )
  colnames( array.indices ) = c( "row", "col" )
  
  libraries.a = character()
  libraries.b = character()
  for ( row.use in 1:nrow( array.indices ) ) {
    
    library.a = rownames( binary.overlap.table )[ array.indices[ row.use, 1 ] ]
    library.b = colnames( binary.overlap.table )[ array.indices[ row.use, 2 ] ]
    
    # message( library.a, " = ", library.b )
    
    libraries.a = c( libraries.a, library.a )
    libraries.b = c( libraries.b, library.b )
    
  }
  
  if( any( duplicated( libraries.a ) ) | any( duplicated( libraries.b ) ) ) { 
    
    stop( "ERROR: multiple matches found by get.corresponding.cell.names" ) 
    
  }
  
  set.a[ library.a %in% set.a.missing, library := NA  ]
  set.a$library.renamed = get.corresponding.labels( set.a$library.a, libraries.a, libraries.b  )
  
  # print( set.a[ , table( library.renamed, useNA = "a" ) ] )
  
  set.a[ !is.na( library.renamed ), library.renamed := paste0( library.renamed, "_", cbc ) ]
  new.cell.names = set.a$library.renamed
  
  outs = list( new.cell.names = new.cell.names, set.a.only = set.a.missing, set.b.only = set.b.missing )
  
  return( outs )
  
}


#'
#'
convert.server.and.windows.paths = function( paths, path.conversions ) {
  
  all.windows = all( grepl( ":", paths ) )
  any.windows = any( grepl( ":", paths ) )
  
  if ( all.windows ) {
    
    for ( i in 1:length( path.conversions ) ) { paths = gsub( names( path.conversions )[ i ], path.conversions[ i ], paths, fixed = TRUE ) }
    
  } else if ( !any.windows ) {
    
    for ( i in 1:length( path.conversions ) ) { paths = gsub( path.conversions[ i ], names( path.conversions )[ i ], paths, fixed = TRUE ) }
    
  } else { message( "WARNING: mixture of server and windows paths ( quitting )" ); break }
  
  return( paths )
  
}


#'
#'
create.loading.colors = function( data, center.col = "white", neg.col = "blue", pos.col = "red", 
                                  center.at = 0, min = NULL, max = NULL ) {
  
  data = data - center.at
  
  if ( is.null( min ) ) { min = -max( abs( data ) ) }
  if ( is.null( max ) ) { max = max( abs( data ) ) }
  
  if ( abs( min ) != abs( max ) ) { message( "WARNING: min and max are different distances from center" ) }
  
  if ( any( data < min ) | any( data > max ) ) { message( "WARNING: data falls outside of min / max range" ) }
  data[ data < min ] = min
  data[ data > max ] = max
  
  center = which( data == 0 )
  negative = which( data < 0 )
  positive = which( data > 0 )
  
  
  center.colors = rep( center.col, length( center ) )
  
  negative.gradient = colorRampPalette( c( neg.col, center.col ) )
  negative.colors = negative.gradient( 100 )[ as.numeric( cut( c( 0, min, data[ negative ] ), breaks = 100 ) ) ]
  negative.colors = negative.colors[ -c( 1:2 ) ]
  
  positive.gradient = colorRampPalette( c( center.col, pos.col ) )
  positive.colors = positive.gradient( 100 )[ as.numeric( cut( c( 0, max, data[ positive ] ), breaks = 100 ) ) ]
  positive.colors = positive.colors[ -c( 1:2 ) ]
  
  colors = vector( length = length( data ) )
  
  colors[ center ] = center.colors
  colors[ negative ] = negative.colors
  colors[ positive ] = positive.colors
  
  data.order = order( abs( data ) )
  out = list( colors = colors, order = data.order )

  return( out )
  
}



# x analyses ------------------------------------------------------------------------------------------------------

#'
#'
load.raw.dac = function( dac.directory ) {
  
  library( data.table )
  
  files = dir( dac.directory, pattern = ".dac.txt.gz" )
  replicate.names = do.call( "rbind", strsplit( files, ".dac.txt.gz" ) )
  dac = lapply( paste0( dac.directory, "/", files ), 
                function( path ) { message( path ); fread( path ) } )
  dac = mapply( function( table, tag ) { table[ , cell := paste( tag, cell, sep = "_" ) ] }, 
                table = dac, tag = replicate.names, SIMPLIFY = FALSE )
  dac = rbindlist( dac )
  # EDIT - include read columns ( i.e. rA ... rN ) too?
  col.keep = c( "chr", "pos", "gene", "cell", "A", "C", "G", "T", "N" )
  dac = dac[ , col.keep, with = FALSE ]
  data.table::setkey( dac, pos )
  
}

#' Load one or more DACs, disambiguate the CBCs based on the file names, and merge them into a single data.table.
#' @param files one or more files with extension .dac.txt.gz, produced by GatherDigitalAlleleCounts.  The colnames should be
#' chr	pos	gene	cell	ref_base	alt_base	rA	rC	rG	rT	rN	A	C	G	T	N	num_umi	umi_mean_purity	read_pval	read_ratio	read_ci_low	read_ci_high	umi_pval	umi_ratio	umi_ci_low	umi_ci_high
#' @import data.table
load_and_merge_dacs<-function(files) {
  replicate.names = do.call( "rbind", strsplit( basename(files), ".dac.txt.gz" ) )
  dac = lapply( files, function( path ) { data.table::fread( path ) } )
  dac = mapply( function( table, tag ) { table[ , cell := paste( tag, cell, sep = "_" ) ] }, 
                table = dac, tag = replicate.names, SIMPLIFY = FALSE )
  dac = data.table::rbindlist( dac )
  col.keep = c( "chr", "pos", "gene", "cell", "A", "C", "G", "T", "N" )
  dac = dac[ , col.keep, with = FALSE ]
  data.table::setkey( dac, pos )
  return(dac)
}

#'
#'  
assign.bases.to.haplotypes = function( row, gt.col = "gt" ) {
  
  row = unlist( row )
  total = sum( as.numeric( row[ c( "A", "C", "G", "T", "N" ) ] ) )
  gt.vector = unlist( strsplit( row[ gt.col ], "[[:punct:]]" ) )
  h1 = as.numeric( row[ gt.vector[ 1 ] ] )
  h2 = as.numeric( row[ gt.vector[ 2 ] ] )
  other = total - sum( h1 + h2 )
  out = c( "H1" = h1, "H2" = h2, "other" = other, "total" = total )
  
}


#'
#'
phase.x.variants = function( unphased.dac, gen.unit.use = "pos", xi.genes = NULL, genes.exclude = NULL, 
                             cells.exclude = NULL, min.cells = 10, max.xi.proportion = 0.2, time.limit = 48,
                             max.gen.units = 1000 ) {
  
  first.start = Sys.time()
  current.time = Sys.time()
  
  copy.unphased.dac = data.table::copy( unphased.dac )
  # copy.unphased.dac = copy( dac )
  copy.unphased.dac[ chr == "X" ]
  
  colnames( copy.unphased.dac ) = gsub( gen.unit.use, "gen.unit", colnames( copy.unphased.dac ) )
  gen.units = copy.unphased.dac[ , sort( unique( gen.unit ) ) ]
  flip.record = data.table::data.table( gen.unit = gen.units, flip = FALSE, contains = as.character( gen.units ), touched = FALSE )
  
  # flip xi.genes and remove unwanted genes
  copy.unphased.dac[ gene %in% xi.genes, `:=`( H1 = H2, H2 = H1 ) ]
  copy.unphased.dac = copy.unphased.dac[ !( gene %in% genes.exclude ) ]
  
  # aggregate to cell.by.gen.unit
  cell.by.gen.unit = copy.unphased.dac[ , lapply( .SD, sum ), by = .( cell, gen.unit ), .SDcols = c( "H1", "H2" ) ] 
  
  # remove unwanted cells 
  cell.by.gen.unit = cell.by.gen.unit[ !( cell %in% cells.exclude ) ] 
  
  # remove cells / gen.unit combos with tied H1 / H2 observations ( which also removes cells with no H1 | H2 observations )
  cell.by.gen.unit = cell.by.gen.unit[ H1 != H2 ]  
  
  # remove cells with single gen.unit ( they're no help )
  cell.by.gen.unit[ , n.gen.unit := .N, by = cell ]
  cell.by.gen.unit = cell.by.gen.unit[ n.gen.unit > 1 ]
  
  # remove gen.units in fewer than min.cells cells ( striking a balance between "completeness" and computational feasibility )  
  cell.by.gen.unit[ , length( unique( gen.unit ) ) ]
  cell.by.gen.unit[ , n.cell := .N, by = gen.unit ]
  cell.by.gen.unit = cell.by.gen.unit[ n.cell >= min.cells ]
  cell.by.gen.unit[ , length( unique( gen.unit ) ) ]
  
  # limit the total number of gen.units
  tmp = cell.by.gen.unit[ , .N, by = gen.unit ]
  colnames( tmp )[ 2 ] = "cells"
  cells.per.gen.unit = tmp$cells
  names( cells.per.gen.unit ) = tmp$gen.unit
  cells.per.gen.unit = sort( cells.per.gen.unit,  decreasing = TRUE )
  if ( length( cells.per.gen.unit ) > max.gen.units ) {
    
    gen.units.use = names( cells.per.gen.unit )[ 1:max.gen.units ]
    
  } else {
    
    gen.units.use = names( cells.per.gen.unit )
    
  }
  cell.by.gen.unit = cell.by.gen.unit[ gen.unit %in% gen.units.use ]
  
  # add option to select the best n represented gen.units? ( instead of a min.cells cutoff )
  # you can do it *after* the list generation steps below ( they're pretty quick )
  
  cell.by.gen.unit = cell.by.gen.unit[ , cell:H2 ]
  
  # split cell.by.gen.unit into list of data.tables by gen.unit and active.x    
  cell.by.gen.unit[ H1 > H2, active.x := "H1" ]
  cell.by.gen.unit[ H1 < H2, active.x := "H2" ]
  data.table::setkey( cell.by.gen.unit, cell, gen.unit, active.x )
  list.by.gen.unit = split( cell.by.gen.unit, by = "gen.unit" )
  list.by.gen.unit.and.active.x = lapply( list.by.gen.unit, function( x ) { split( x, by = "active.x" ) } )
  
  # try using # of observations instead of # of cells? ( or is it not worth it? )
  
  gen.units = cell.by.gen.unit[ , sort( unique( gen.unit ) ) ]
  n.gen.units = length( gen.units )
  
  message( "attempting to phase ", n.gen.units, " gen.units" )
  gen.unit.pairs = combn( gen.units, 2 )
  
  message( "generating initial concordance matrix ( slow )" )
  start = Sys.time()
  concordance.vector = apply( gen.unit.pairs, 2, get.pair.concordance, 
                              list.a = list.by.gen.unit, 
                              list.b = list.by.gen.unit.and.active.x, 
                              min.cells = min.cells,
                              max.xi.proportion = max.xi.proportion )
  end = Sys.time()
  print( end - start )
  
  # format concordance vector into concordance matrix
  concordance.matrix = matrix( nrow = n.gen.units, ncol = n.gen.units, dimnames = list( gen.units, gen.units ) )
  concordance.matrix[ lower.tri( concordance.matrix, diag = FALSE ) ] = concordance.vector
  concordance.matrix = t( concordance.matrix )
  
  message( "starting loop" )
  termination.message = "merged all gen.units"
  start = Sys.time()
  while( length( unique ( cell.by.gen.unit$gen.unit ) ) > 1 & ( as.numeric( current.time - first.start, units = "hours" ) < time.limit ) ) {
    
    # check if concordance matrix is all null
    if ( all( is.na( concordance.matrix ) ) ) {
      
      termination.message = "no evidence remaining"; break 
      
    }
    
    # get index of smallest value in concordance matrix and names of corresponding gen.units
    min.index = which( concordance.matrix == min( concordance.matrix, na.rm = T ), arr.ind = TRUE )
    min.index = matrix( min.index[ 1, ], nrow = 1 )
    gen.unit.pair.use = c( colnames( concordance.matrix )[ min.index[ 1 ] ], rownames( concordance.matrix )[ min.index[ 2 ] ] )
    gen.unit.pair.use = as.integer( gen.unit.pair.use )
    gen.units.contained = unlist( strsplit( flip.record[ gen.unit == gen.unit.pair.use[ 2 ], contains ], ",", fixed = TRUE ) )
    
    # get concordance.differential for selected pair
    pair = as.character( gen.unit.pair.use )
    shared = sum( list.by.gen.unit[[ pair[ 1 ] ]]$cell %in% list.by.gen.unit[[ pair[ 2 ] ]]$cell )
    H1 = sum( list.by.gen.unit.and.active.x[[ pair[ 1 ] ]]$H1$cell %in% list.by.gen.unit.and.active.x[[ pair[ 2 ] ]]$H1$cell )
    H2 = sum( list.by.gen.unit.and.active.x[[ pair[ 1 ] ]]$H2$cell %in% list.by.gen.unit.and.active.x[[ pair[ 2 ] ]]$H2$cell )
    concordant = H1 + H2
    discordant = shared - concordant
    concordance.differential = concordant - discordant
    
    # update pair phase in cell.by.gen.unit
    if( concordance.differential < 0 ) { 
      
      cell.by.gen.unit[ gen.unit == gen.unit.pair.use[ 2 ], ':='( H1 = H2, H2 = H1 ) ]
      flip.record[ gen.unit %in% gen.units.contained, flip := !flip ]
      message.tail = " - flipped!"
      
    } else {  message.tail = "" }
    cell.by.gen.unit[ gen.unit %in% gen.unit.pair.use[ 2 ], gen.unit := gen.unit.pair.use[ 1 ] ]
    
    # update flip.record
    flip.record[ gen.unit == gen.unit.pair.use[ 1 ], contains := paste( c( contains, gen.units.contained ), collapse = "," ) ]
    flip.record[ gen.unit == gen.unit.pair.use[ 2 ], contains := NA ]
    flip.record[ gen.unit %in% gen.unit.pair.use, touched := TRUE ]
    
    message( gen.unit.use, " ", gen.unit.pair.use[ 1 ], " & ", gen.unit.pair.use[ 2 ], " merged w/ ", concordant, " versus ", 
             discordant, " ( ", round( min( concordant, discordant ) / sum( concordant, discordant ), 2 ), " )", message.tail )
    
    # aggregate the now merged gen.units observations in cell.by.gen.unit
    cell.by.gen.unit = cell.by.gen.unit[ , lapply( .SD, sum ), by = .( cell, gen.unit ), .SDcols = c( "H1", "H2" ) ]
    
    # split cell.by.gen.unit into list of data.tables by gen.unit and active.x ( updating the previous version of these lists )
    P1 = as.character( gen.unit.pair.use[ 1 ] )
    P2 = as.character( gen.unit.pair.use[ 2 ] )
    cell.by.gen.unit.subset = cell.by.gen.unit[ gen.unit == gen.unit.pair.use[ 1 ] ]
    cell.by.gen.unit.subset[ H1 > H2, active.x := "H1" ]
    cell.by.gen.unit.subset[ H1 < H2, active.x := "H2" ]
    data.table::setkey( cell.by.gen.unit.subset, cell, gen.unit, active.x )
    list.by.gen.unit[[ P1 ]] = cell.by.gen.unit.subset
    list.by.gen.unit[[ P2 ]] = NULL
    list.by.gen.unit.and.active.x[[ P1 ]] = split( cell.by.gen.unit.subset, by = "active.x" )
    list.by.gen.unit.and.active.x[[ P2 ]] = NULL
    
    if ( nrow( concordance.matrix ) == 3 ) { 
      
      termination.message = "reached end of concordance matrix"; break
      
    }
    
    # update the concordance matrix for start of next iteration
    concordance.matrix = concordance.matrix[ , !( colnames( concordance.matrix ) %in% as.character( gen.unit.pair.use ) ) ]
    concordance.matrix = concordance.matrix[ !( rownames( concordance.matrix ) %in% as.character( gen.unit.pair.use ) ), ]
    new.gen.unit.pairs = rbind( rep( gen.unit.pair.use[ 1 ], nrow( concordance.matrix ) ), rownames( concordance.matrix ) )
    new.concordance.vector = apply( new.gen.unit.pairs, 2, get.pair.concordance, 
                                    list.a = list.by.gen.unit, 
                                    list.b = list.by.gen.unit.and.active.x, 
                                    min.cells = min.cells,
                                    max.xi.proportion = max.xi.proportion )
    concordance.matrix = cbind( concordance.matrix, new.concordance.vector )
    concordance.matrix = rbind( concordance.matrix, NA )
    colnames( concordance.matrix )[ ncol( concordance.matrix ) ] = gen.unit.pair.use[ 1 ]
    rownames( concordance.matrix )[ nrow( concordance.matrix ) ] = gen.unit.pair.use[ 1 ]
    
    current.time = Sys.time()
    
  }
  
  if( as.numeric( current.time - first.start, units = "hours" ) > time.limit ) { termination.message = "time.limit reached" }
  
  message( termination.message )
  
  end = Sys.time()
  print( end - start )
  
  # additions to flip.record
  flip.record[ which( !touched ), flip := NA ]
  # count gen.units contained by each
  flip.record[ , num.gen.units := stringr::str_count( contains, pattern = "," ) + 1 ]
  
  # EDIT: make this faster by first making a n.obs summary table for each position
  
  total.counts.by.gen.unit = copy.unphased.dac[ , sum( total ), by = "gen.unit" ]
  colnames( total.counts.by.gen.unit ) = c( "gen.unit", "total" )
  
  tmp.function = function( contains ) {
    
    gen.units.contained = strsplit( contains, ",", fixed = TRUE )[[ 1 ]]
    total.observations = total.counts.by.gen.unit[ gen.unit %in% gen.units.contained, sum( total ) ]
    
  }
  flip.record[ !is.na( contains ), total.observations := tmp.function( contains ), by = "gen.unit" ] # WARNING: fix this! ( not robust to other gen.units )
  flip.record[ , pct.total := round( total.observations / sum( total.observations, na.rm = TRUE ), 3 ) ]
  colnames( flip.record )[ 1 ] = gen.unit.use
  
  return( flip.record )
  
}


#'
#'
get.pair.concordance = function( pair, list.a, list.b, min.cells, max.xi.proportion ) {
  
  pair = as.character( pair )
  
  shared = sum( list.a[[ pair[ 1 ] ]]$cell %in% list.a[[ pair[ 2 ] ]]$cell )
  if ( shared < min.cells ) { return( NA ) }
  
  H1 = sum( list.b[[ pair[ 1 ] ]]$H1$cell %in% list.b[[ pair[ 2 ] ]]$H1$cell )
  H2 = sum( list.b[[ pair[ 1 ] ]]$H2$cell %in% list.b[[ pair[ 2 ] ]]$H2$cell )
  
  concordant = H1 + H2
  discordant = shared - concordant
  successes = min( concordant, discordant )
  
  if ( ( successes / shared ) > max.xi.proportion ) { return( NA ) }
  
  binom::binom.confint( successes, shared, methods = "exact" )[ , "upper" ]
  
}


#'
#'
generate.error.rates.and.update.phasing = function( dac, phased.positions, flip.positions, 
                                                    min.cells = 1, min.umi = 1, 
                                                    calling.count.min = 1, calling.purity.min = 0.8, 
                                                    likelihood.min = 0.8, position.uncertainties = NULL,
                                                    flip.xist = FALSE ) {
  
  dac.copy = data.table::copy( dac )
  dac.copy$total = NULL
  dac.copy[ , total := H1 + H2 ]
  dac.copy = dac.copy[ total > 0 ]
  
  if ( flip.xist ) { dac.copy[ gene %in% c( "XIST", "TSIX" ), `:=`( H1 = H2, H2 = H1 ) ] }
  
  # remove any positions with insufficient observation
  cells.per.position = dac.copy[ , .( total = sum( total ) ), by = c( "cell", "pos" ) ]
  cells.per.position = cells.per.position[ , .( cells = sum( total > 0 ) ), by = "pos" ]
  positions.use = sort( cells.per.position[ cells >= min.cells, pos ] )
  umis.per.position = dac.copy[ , .( total = sum( total ) ), by = "pos" ]
  positions.use = positions.use[ positions.use %in% umis.per.position[ total >= min.umi, pos ] ]
  
  n.pos = length( positions.use )
  message( "using ", n.pos, " positions ( ", round( n.pos / length( unique( dac.copy$pos ) ), 2 ), " )" )
  
  
  position.stats = data.table::data.table()
  
  # progress.bar = txtProgressBar( min = 0, max = length( positions.use ), initial = 0, style = 3 )
  # i = 1
  
  for ( position.use in positions.use ) {
    
    position.dac = data.table::copy( dac.copy )
    position.dac = position.dac[ pos == position.use ]
    position.dac = position.dac[ , .( pos, ref, alt, gene, cell, A1 = H1, A2 = H2 ) ]
    
    ref = position.dac[ 1, ref ]
    alt = position.dac[ 1, alt ]
    gene = position.dac[ 1, gene ]
    
    dac.subset = data.table::copy( dac.copy )
    dac.subset = dac.subset[ pos != position.use ]
    
    x.calls = call.active.x( dac.subset, phased.positions, flip.positions, position.uncertainties = position.uncertainties, 
                             calling.count.min = calling.count.min, calling.purity.min = calling.purity.min, 
                             likelihood.min = likelihood.min, flip.xist = flip.xist, verbose = FALSE )
    
    x.calls = x.calls[ , .( cell, active.x ) ]
    
    position.dac = merge( position.dac, x.calls[ !is.na( active.x ) ], by = "cell" )
    
    n.cells = length( unique( position.dac$cell ) )
    
    
    # calculate phasing misses and matches for position
    
    A1.mis = position.dac[ active.x == "X2", sum( A1 ) ]
    A1.match = position.dac[ active.x == "X1", sum( A1 ) ]
    A1.total = A1.mis + A1.match
    
    A2.mis = position.dac[ active.x == "X1", sum( A2 ) ]
    A2.match = position.dac[ active.x == "X2", sum( A2 ) ]
    A2.total = A2.mis + A2.match
    
    mis = A1.mis + A2.mis
    match = A1.match + A2.match
    total = A1.total + A2.total
    error.rate = mis / total
    
    if ( is.nan( error.rate ) ) { error.rate = NA }
    
    
    # EDIT - use allele specific error.rates? ( some example code included below )
    # x.ratio = X1.cells / X2.cells
    # A1.error.rate = ( A1.mis * x.ratio ) / ( A1.match + A1.mis * x.ratio )
    # A2.error.rate = ( A2.mis / x.ratio ) / ( A2.match + A2.mis / x.ratio )

    
    previously.phased = position.use %in% phased.positions
    
    if ( previously.phased ) {
      
      previously.flip = position.use %in% flip.positions
      
    } else { previously.flip = NA }
    
    
    new.row = data.table::data.table( pos = position.use, ref, alt, gene, 
                          previously.phased, previously.flip,
                          n.cells, 
                          mis, match, total, error.rate )
    
    position.stats = rbind( position.stats, new.row )
    
    
    # setTxtProgressBar( progress.bar, i )
    # i = i + 1 
    
  }
  
  # EDIT - try min.uncertainty alternatives? ( probably not )
  #   currently using 1 / ( n.obs + 1 )

  position.stats[ total == 0, min.uncertainty := 0.5, by = "pos" ]
  position.stats[ total != 0, min.uncertainty := 1 / ( total + 1 ), by = "pos" ]
  
  position.stats[ error.rate < 0.5, uncertainty := max( error.rate, min.uncertainty, na.rm = TRUE ), by = "pos" ]
  position.stats[ error.rate > 0.5, uncertainty := max( 1 - error.rate, min.uncertainty, na.rm = TRUE ), by = "pos" ]
  position.stats[ error.rate == 0.5, uncertainty := 0.5 ]
  position.stats[ is.na( error.rate ), uncertainty := 0.5 ]
  
  position.stats[ uncertainty > 0.5, ] # should be empty

  return( position.stats )
  
}


#'
#'
call.active.x = function( dac, phased.positions, flip.positions, position.uncertainties = NULL, 
                          calling.count.min = 1, calling.purity.min = 0.8, likelihood.min = 0.8,
                          flip.xist = FALSE, verbose = TRUE ) {
  
  # using likelihood and weighted.purity for x calls when available ( i.e. when position.uncertainties known )
  
  dac.copy = data.table::copy( dac )
  colnames( dac.copy ) = gsub( "H1", "X1", colnames( dac.copy ), fixed = TRUE )
  colnames( dac.copy ) = gsub( "H2", "X2", colnames( dac.copy ), fixed = TRUE )
  
  if ( flip.xist ) { dac.copy[ gene %in% c( "XIST", "TSIX" ), `:=`( X1 = X2, X2 = X1 ) ] }
  
  dac.copy = dac.copy[ X1 + X2 > 0 ]
  dac.copy = dac.copy[ pos %in% phased.positions ]
  
  if ( is.null( position.uncertainties ) ) {
    
    if ( verbose ) { message( "using counts ( no position uncertainties provided )" ) }
    
    dac.copy[ pos %in% flip.positions, `:=`( X1 = X2, X2 = X1 ) ]
    
    x.calls = dac.copy[ , lapply( .SD, sum ), by = .( cell ), .SDcols = c( "X1", "X2" ) ]
    
    x.calls[ , total := X1 + X2 ]
    x.calls[ , purity := round( max( X1, X2 ) / total, 2 ), by = 1:nrow( x.calls ) ]
    
    x.calls[ X1 > X2 & X1 >= calling.count.min & purity >= calling.purity.min, active.x := "X1" ]
    x.calls[ X2 > X1 & X2 >= calling.count.min & purity >= calling.purity.min, active.x := "X2" ]
    
    return( x.calls )

  } else if ( !is.null( position.uncertainties ) ) {
    
    if ( verbose ) { message( "using position uncertainties" ) }

    dac.copy = merge( dac.copy, position.uncertainties[ !is.na( uncertainty ) ], by = c( "pos", "gene" ) )
    
    dac.copy[ , X1.weighted := X1 * ( 1 - uncertainty ), by = 1:nrow( dac.copy ) ]
    dac.copy[ , X2.weighted := X2 * ( 1 - uncertainty ), by = 1:nrow( dac.copy ) ]
    
    dac.copy[ , X1.probability := prod( c( rep( 1 - uncertainty, X1 ), rep( uncertainty, X2 ) ) ), by = 1:nrow( dac.copy ) ]
    dac.copy[ , X2.probability := prod( c( rep( uncertainty, X1 ), rep( 1 - uncertainty, X2 ) ) ), by = 1:nrow( dac.copy ) ]
    
    dac.copy[ pos %in% flip.positions, `:=`( X1 = X2, X2 = X1 ) ]
    dac.copy[ pos %in% flip.positions, `:=`( X1.weighted = X2.weighted, X2.weighted = X1.weighted ) ]
    dac.copy[ pos %in% flip.positions, `:=`( X1.probability = X2.probability, X2.probability = X1.probability ) ]
    
    x.calls.sum = dac.copy[ , lapply( .SD, sum ), by = .( cell ), .SDcols = c( "X1", "X2", "X1.weighted", "X2.weighted" ) ]
    x.calls.prod = dac.copy[ , lapply( .SD, prod ), by = .( cell ), .SDcols = c( "X1.probability", "X2.probability" ) ]
    x.calls = merge( x.calls.sum, x.calls.prod, by = "cell" )
    
    x.calls[ , X1.likelihood := X1.probability / ( X1.probability + X2.probability ), by = 1:nrow( x.calls ) ]

    x.calls[ , purity := max( X1, X2 ) / ( X1 + X2 ), by = "cell" ]
    x.calls[ , purity.weighted := max( X1.weighted, X2.weighted ) / ( X1.weighted + X2.weighted ), by = "cell" ]
    
    x.calls[ X1.likelihood > likelihood.min & purity.weighted > calling.purity.min, active.x := "X1" ]
    x.calls[ ( 1 - X1.likelihood ) > likelihood.min & purity.weighted > calling.purity.min, active.x := "X2" ]
    
    return( x.calls )

  }

}


#'
#'
get.corrected.skew.variances = function( skew.matrix, total.matrix ) {
  
  error.matrix = skew.matrix * ( 1 - skew.matrix ) / total.matrix
  
  correction = apply( error.matrix, 2, sum, na.rm = TRUE )
  correction = correction / ( nrow( error.matrix ) - 1 )
  
  variance = apply( skew.matrix, 2, var, na.rm = TRUE )
  corrected.variance = variance - correction
  
  return( list( variances = variance, corrected.variances = corrected.variance ) )
  
}


#'
#'
get.corrected.skew.correlations = function( skew.matrix, total.matrix, cor.method.use = "pearson" ) {
  
  var.out = get.corrected.skew.variances( skew.matrix, total.matrix )
  
  variance = var.out$variances
  corrected.variance = var.out$corrected.variances
  
  
  cor.matrix = cor( skew.matrix, method = cor.method.use, use = "pairwise.complete.obs" )
  
  corrected.cor.matrix = cor.matrix
  
  for ( p in colnames( cor.matrix ) ) {
    
    for ( q in colnames( cor.matrix ) ) {
      
      cor.correction = sqrt( variance[ p ] * variance[ q ] ) / sqrt( corrected.variance[ p ] * corrected.variance[ q ] )
      corrected.cor.matrix[ p, q ] = cor.matrix[ p, q ] * cor.correction
      
    }
    
  }
  
  return( list( correlations = cor.matrix, corrected.correlations = corrected.cor.matrix ) )
  
}


#'
#'
process.seurat.diff.exp.output = function( diff.exp, chr.table ) {
  
  diff.exp = data.table::as.data.table( diff.exp, keep.rownames = "gene" )
  diff.exp = merge( diff.exp, chr.table[ , .( chr, gene ) ], by = "gene", all.x = TRUE )
  colnames( diff.exp ) = c( "gene", "p.value", "log.2.fc", "pct.1", "pct.2", "p.value.adj", "chr" )
  setorderv( diff.exp, c( "p.value" ) )
  diff.exp[ , fold.change := 2 ^ log.2.fc ]
  diff.exp$q.value = qvalue( diff.exp$p.value )$qvalues
  diff.exp = diff.exp[ , .( chr, gene, pct.1, pct.2, fold.change, log.2.fc, p.value, q.value, bonferroni = p.value.adj ) ]
  diff.exp[ , chr := gsub( "chr", "", chr, fixed = TRUE ) ]
  diff.exp
  
  return( diff.exp )
  
}


#'
#'
equalize.x.proportions.by.cluster = function( input.table, plot = TRUE ) {
  
  by.cluster = data.table::copy( input.table )
  by.cluster[ active.x == "X1", X1.cells := 1 ]
  by.cluster[ active.x == "X2", X2.cells := 1 ]
  by.cluster = by.cluster[ , lapply( .SD, sum, na.rm = TRUE ), by = .( cluster ), .SDcols = ( X1.cells:X2.cells ) ]
  setorder( by.cluster, cluster )
  
  if ( nrow( by.cluster ) == 1 ) {
    
    message( "only one cluster present ( skiping )" )
    return( input.table$cell )
    
  }
  
  by.cluster[ , x.total := X1.cells + X2.cells ]
  
  target.proportion = by.cluster[ , sum( X1.cells ) / sum( x.total ) ]
  message( "target proportion is ", round( target.proportion, 3 ) )
  
  equalized.table = data.table::copy( input.table )
  equalized.table$remove = FALSE
  for ( i in 1:nrow( by.cluster ) ) {
    
    cluster.to.equalize = by.cluster[ i, cluster ]
    current.X1.count = equalized.table[ cluster == cluster.to.equalize & active.x == "X1", .N ]
    current.X2.count = equalized.table[ cluster == cluster.to.equalize & active.x == "X2", .N ]
    
    if ( by.cluster[ i, X1.cells / x.total ] < target.proportion ) {
      
      n.to.remove = round( ( current.X1.count + current.X2.count ) - ( current.X1.count / target.proportion ) )
      to.remove = equalized.table[ cluster == cluster.to.equalize & active.x == "X2" ][ sample( 1:.N, n.to.remove ), cell ]
      equalized.table[ cell %in% to.remove, remove := TRUE ]
      
    } else if ( by.cluster[ i, X1.cells / x.total ] > target.proportion ) {
      
      n.to.remove = round( ( current.X1.count - target.proportion * current.X1.count - target.proportion * current.X2.count ) / ( 1 - target.proportion ) )
      to.remove = equalized.table[ cluster == cluster.to.equalize & active.x == "X1" ][ sample( 1:.N, n.to.remove ), cell ]
      equalized.table[ cell %in% to.remove, remove := TRUE ]
      
    } else { message( "already equal" ) }
    
  }
  
  removal.table = equalized.table[ , table( cluster, remove ) ]
  removal.table.prop = round( removal.table / apply( removal.table, 1, sum ), 2 )
  removal.counts = apply( removal.table, 2, sum )
  message( "removing ", round( removal.counts[ 2 ] / sum( removal.counts ), 2 ), " of data" )
  
  equalized.table = equalized.table[ !( remove ), ]
  
  if ( plot ) {
    
    by.cluster.equalized = data.table::copy( equalized.table )
    by.cluster.equalized[ active.x == "X1", X1.cells := 1 ]
    by.cluster.equalized[ active.x == "X2", X2.cells := 1 ]
    by.cluster.equalized = by.cluster.equalized[ , lapply( .SD, sum, na.rm = TRUE ), by = .( cluster ), .SDcols = ( X1.cells:X2.cells ) ]
    by.cluster.equalized[ , x.total := X1.cells + X2.cells ]
    setorderv( by.cluster.equalized, c( "cluster" ) )
    
    by.cluster.merged = merge( by.cluster, by.cluster.equalized, by = "cluster" )
    
    xlim.use = range( unlist( by.cluster.merged[ , .( x.total.x, x.total.y ) ] ) )
    
    by.cluster[ , plot( x.total, X1.cells / x.total, col = "white", xlim = xlim.use, ylim = c( 0, 1 ), ylab = "Proportion X1", xlab = "# cells" ) ]
    abline( h = target.proportion )
    by.cluster.merged[ , arrows( x.total.x, X1.cells.x / x.total.x, x.total.y, X1.cells.y / x.total.y, length = 0, lty = 2 ) ]
    by.cluster[ , points( x.total, X1.cells / x.total, pch = 20 ) ]
    by.cluster.equalized[ , points( x.total, X1.cells / x.total ) ]
    by.cluster.equalized[ , points( x.total, X1.cells / x.total, pch = 20, col = "white" ) ]
    by.cluster[ X1.cells / x.total > target.proportion, text( x.total, X1.cells / x.total, labels = cluster, pos = 3 ) ]
    by.cluster[ X1.cells / x.total < target.proportion, text( x.total, X1.cells / x.total, labels = cluster, pos = 1 ) ]
    
    title( main = paste( input.table[ , .N ], "|", equalized.table[ , .N ], "(", round( equalized.table[ , .N ] / input.table[ , .N ], 2 ), ")" ) )
    
  }
  
  message( input.table[ , .N ], " starting" )
  message( equalized.table[ , .N ], " ending ( ", round( equalized.table[ , .N ] / input.table[ , .N ], 2 ), " )" )
  
  return( equalized.table$cell )
  
}


#'
#'
process.seurat.diff.exp.output = function( diff.exp, chr.table ) {
  
  diff.exp = data.table::as.data.table( diff.exp, keep.rownames = "gene" )
  diff.exp = merge( diff.exp, chr.table[ , .( chr, gene ) ], by = "gene", all.x = TRUE )
  colnames( diff.exp ) = c( "gene", "p.value", "log.2.fc", "pct.1", "pct.2", "p.value.adj", "chr" )
  setorderv( diff.exp, c( "p.value" ) )
  diff.exp[ , fold.change := 2 ^ log.2.fc ]
  diff.exp$q.value = qvalue( diff.exp$p.value )$qvalues
  diff.exp = diff.exp[ , .( chr, gene, pct.1, pct.2, fold.change, log.2.fc, p.value, q.value, bonferroni = p.value.adj ) ]
  diff.exp[ , chr := gsub( "chr", "", chr, fixed = TRUE ) ]
  diff.exp
  
  return( diff.exp )
  
}


#'
#'
liftover.x.positions = function( positions.to.liftover, chain.object ) {
  
  pos.dt = data.table::data.table( pos = positions.to.liftover, group = 1:length( positions.to.liftover ) )
  
  granges.object = pos.dt[ , GRanges( seqnames = "chrX", ranges = IRanges( start = pos, end = pos ) ) ]
  liftover.results = data.table::as.data.table( liftOver( granges.object, chain.object ) )
  pos.dt = merge( pos.dt, liftover.results[ , .( group, new.pos = start ) ], all.x = TRUE, by = "group" )
  pos.dt$group = NULL
  failed.liftover = pos.dt[ is.na( new.pos ), pos ]
  prop.failed = length( failed.liftover ) / length( positions.to.liftover )
  
  message( length( failed.liftover ), " of ", length( positions.to.liftover ), " failed liftover ( ", round( prop.failed, 3 ), " )" )
  
  return( pos.dt$new.pos )
  
}



# clustering ------------------------------------------------------------------------------------------------------

#'
#'
get.library.paths.and.metadata = function() {
  
  if( file.exists( "library_paths_and_metadata.txt" ) ) {
    
    message( "using existing library_paths_and_metadata.txt" )
    library.paths.and.metadata = fread( "library_paths_and_metadata.txt" )
    
  } else { 
    
    message( "copy my.name, cell.features, and sparse.dge columns ( press enter when ready )" )
    readline()
    library.paths.and.metadata = data.table::as.data.table( read.table( "clipboard", header = TRUE, colClasses = "character" ) )
    write.table( library.paths.and.metadata, "library_paths_and_metadata.txt", row.names = FALSE, quote = FALSE )
    message( "wrote new library_paths_and_metadata.txt" )
    
  }
  
  return( library.paths.and.metadata )
  
}


#'
#'
create.seurat.for.global.clustering = function( library.paths.and.metadata, clustering.dir, exclude.mito.genes = FALSE,
                                                n.features = 10000, cells.use = NULL ) {
  
  # determine if a new clustering should be started
  if ( file.exists( "seurat.RDS" ) ) {
    
    message( "start a new clustering? ( y / n )" )
    if ( readline() == "y" ) { new.clustering = TRUE } else { new.clustering = FALSE }
    
  } else { new.clustering = TRUE }
  
  # get seurat object
  if ( !new.clustering ) {
    
    message( "loading existing seurat" )
    seurat = readRDS( "seurat.RDS" )
    
  } else {
    
    message( "loading dge" )
    
    dge.list = as.list( unique( library.paths.and.metadata$sparse.dge ) )
    
    if( length( dge.list ) > 1 ) {
      
      dge.list = lapply( dge.list, readRDS )
      dge.list = lapply( dge.list, t )
      dge = rBind.fill( dge.list ) # Matrix.utils function ( no longer maintained on Cran )
      dge = as( dge, "sparseMatrix" )
      dge = t( dge )
      
    } else { 
      
      dge = readRDS( dge.list[[ 1 ]] )
      
    }
    
    gc()
    
    # keep only a subset of cells
    if ( !is.null( cells.use ) ) {
      
      message( "using only specified cells ( ", length( cells.use ), " )" )
      dge = dge[ , cells.use ]
      
    }
    
    # remove excluded barcodes before clustering?
    if( file.exists( paste0( root.dir, "/excluded_cell_barcodes.txt" ) ) ) {
      
      message( "remove excluded cell barcodes? ( y / n )" )
      if ( readline() == "y" ) { 
        
        excluded.barcodes = scan( paste0( root.dir, "/excluded_cell_barcodes.txt" ), what = "character" )
        dge = dge[ , !( colnames( dge ) %in% excluded.barcodes ) ]
        
      }
      
    }
    
    if ( exclude.mito.genes ) {
      
      message( "excluding mito genes" )
      mito.genes = grep( "^MT-", rownames( dge ), value = TRUE )
      dge = dge[ !( rownames( dge ) %in% mito.genes ), ]
      
    }
    
    # create basic seurat object ( i.e. through pca )
    message( "creating seurat object" )
    seurat = CreateSeuratObject( counts = dge )
    
    message( "running sctransform" )
    message( "using ", n.features, " features" ) # default is 3000 ( for SCTransform )
    if ( n.features == "all" ) {
      
      seurat = SCTransform( seurat, method = "glmGamPoi",
                            variable.features.n = nrow( seurat@assays$RNA@counts ), verbose = FALSE )
      
    } else {
      
      seurat = SCTransform( seurat, method = "glmGamPoi", 
                            variable.features.n = n.features, verbose = FALSE )

    }
    
    message( "running pca" )
    seurat = RunPCA( seurat, npcs = 100, verbose = FALSE )
    
    message( "saving seurat object" )
    saveRDS( seurat, "seurat.RDS" )
    
  }
  
  return( seurat )
  
}


#'
#'
update.cell.feature.columns = function( cell.features, dge ) {
  
  # EDIT: only used in generate_dges_and_cell_features.R as part of the cerebro pipeline ( just copy this code to script? )
  
  colnames( cell.features ) = gsub( "_", ".", colnames( cell.features ), fixed = TRUE )
  colnames( cell.features ) = gsub( "cell.barcode", "cbc", colnames( cell.features ), fixed = TRUE )
  cell.features[ , cell := paste0( library, "_", cbc ) ]
  
  cell.features = cell.features[ , .( cell, library, num.transcripts, num.reads, pct.ribosomal, pct.coding, pct.intronic, 
                                      pct.intergenic, pct.utr, pct.genic, pct.mt, ribomito ) ]
  
  cell.features[ , pct.intronic := pct.intronic / ( pct.intronic + pct.genic ) ]
  cell.features[ , pct.utr := pct.utr / ( pct.intronic + pct.genic ) ]
  cell.features[ , pct.mt := pct.mt / ( pct.intronic + pct.genic ) ]
  cell.features[ , efficiency := num.transcripts / num.reads ]
  
  cell.features = cell.features[ cell %in% colnames( dge ) ]
  
  binary.dge = dge > 0
  n.genes.dt = data.table::data.table( cell = colnames( binary.dge ), num.genes = colSums( binary.dge ) )
  cell.features = merge( cell.features, n.genes.dt, by = "cell" )
  cell.features[ , genes.per.transcript := num.genes / num.transcripts ]
  setorderv( cell.features, "num.transcripts", order = -1 )
  
}


#'
#'
load.cell.features = function( cell.feature.paths, names.use, dge, subset.paths = TRUE ) {
  
  cell.feature.paths = unique( cell.feature.paths )
  
  message( "loading cell features" )
  cell.features = data.table::data.table()
  
  for ( i in 1:length( cell.feature.paths ) ) {
    
    message( names.use[ i ] )
    path = cell.feature.paths[ i ]
    
    if ( subset.paths ) {
      
      message( "subsetting path" )
      
      # break path into shorter components
      path.components = unlist( strsplit( path, "/", fixed = TRUE ) )
      base.path = paste( path.components[ 1:( length( path.components ) - 1 ) ], collapse = "/" )
      file.name = path.components[ length( path.components ) ]
      
      system.command = paste0( "subst M: ", base.path )
      system( system.command )
      
      if ( !file.exists( paste0( "M://", file.name ) ) ) { message( path, " not found" ) }
      
      # load cell features ( or cell selection report ) using path subsets
      if ( grepl( ".cell_features.RDS", file.name, fixed = TRUE ) ) {
        
        cell.features.new = readRDS( paste0( "M://", file.name ) )
        
      } else if ( grepl( ".cell_selection_report.txt", file.name, fixed = TRUE ) ) { 
        
        cell.features.new = fread( paste0( "M://", file.name ) ) 
        
      } else { message( "ERROR: input format not recognized" ); return() }
      
      system( "subst M: /D" )
      
    } else {
      
      message( "using full path" )
      
      # load cell features ( or cell selection report ) using full path
      if ( grepl( ".cell_features.RDS", path, fixed = TRUE ) ) {
        
        cell.features.new = readRDS( path )
        
      } else if ( grepl( ".cell_selection_report.txt", path, fixed = TRUE ) ) { 
        
        cell.features.new = fread( path ) 
        
      } else { message( "ERROR: input format not recognized"); return() }
      
    }
    
    cell.features.new = data.table::as.data.table( cell.features.new )
    colnames( cell.features.new ) = gsub( "_", ".", colnames( cell.features.new ), fixed = TRUE )
    colnames( cell.features.new ) = gsub( "cell.barcode", "cbc", colnames( cell.features.new ), fixed = TRUE )
    cell.features.new = cell.features.new[ , cbc:pct.mt ]
    cell.features.new$replicate = names.use[ i ]
    
    cell.features.new[ , cell := paste0( names.use[ i ], "_", cbc ) ]
    cell.features = rbind( cell.features, cell.features.new )
    
  }
  
  cell.features[ , replicate := factor( replicate ) ]
  setorderv( cell.features, c( "replicate", "num.transcripts" ), order = c( 1, -1 ) )
  
  # update and add a few columns
  cell.features[ , pct.intronic := pct.intronic / ( pct.intronic + pct.genic ) ]
  cell.features[ , pct.utr := pct.utr / ( pct.intronic + pct.genic ) ]
  cell.features[ , pct.mt := pct.mt / ( pct.intronic + pct.genic ) ]
  cell.features[ , efficiency := num.transcripts / num.reads ]
  
  binary.dge = dge > 0
  n.genes.dt = data.table::data.table( cell = colnames( binary.dge ), num.genes = colSums( binary.dge ) )
  
  rep.table.a = cell.features[ , table( replicate ) ]
  cell.features = merge( cell.features, n.genes.dt, by = "cell" )
  rep.table.b = cell.features[ , table( replicate ) ]
  
  if ( any( rep.table.b / rep.table.a == 0 ) ) { message( "WARNING: replicates lost upon merging with dge" ) }
  
  cell.features[ , genes.per.transcript := num.genes / num.transcripts ]
  
  return( cell.features )
  
}


#'
#'
load.cell.features.plot.and.save = function( cell.feature.paths, names.use, dge, subset.paths = TRUE ) {
  
  if( file.exists( "cell_features.RDS" ) ) {
    
    message( "loading existing cell.features" )
    cell.features = readRDS( "cell_features.RDS" )
    
  } else {
    
    cell.features = load.cell.features( cell.feature.paths, names.use, dge, subset.paths )
    
  }
  
  if ( !dir.exists( "cell_feature_pdfs" ) ) { dir.create( "cell_feature_pdfs" ) }
  setwd( "cell_feature_pdfs" )
  
  options( warn = 0 ) 
  
  for( rep in unique( cell.features$replicate ) ) {
    
    pdf( paste0( rep, ".pdf" ) )
    message( "creating cell feature pdf for ", rep )
    
    cols.exclude = c( "cell", "cbc", "replicate" )
    
    for ( variable in colnames( cell.features )[ !( colnames( cell.features ) %in% cols.exclude ) ] ) {
      
      values = cell.features[ , ..variable ]
      variance = var( values )
      
      if ( variance > 1 ) {
        
        cell.features[ replicate == rep, smoothScatter( log10( num.transcripts ), log10( get( variable ) ), main = variable, 
                                                        xlab = "log10( num.transcripts )", 
                                                        ylab = paste0( "log10( ", variable, " )" ) ) ]
        
      } else {
        
        cell.features[ replicate == rep, smoothScatter( log10( num.transcripts ), get( variable ), main = variable,
                                                        xlab = "log10( num.transcripts )", 
                                                        ylab = variable ) ]
        
      }
      
    }
    
    dev.off()
    
  }  
  
  setwd( ".." )
  
  message( "saving cell features" )
  saveRDS( cell.features, "cell_features.RDS" )
  
  return( cell.features )
  
}


#'
#'
create.cell.metadata.table = function( seurat ) {
  
  cell.metadata = data.table::data.table( cell = colnames( seurat@assays$RNA@counts ) )
  cell.metadata[ , random := sample( .N ) ]
  cell.metadata[ , replicate := factor( get.cell.barcode.prefix( cell ) ) ]
  
  return( cell.metadata )
  
}


#'
#'
create.2D.embeddings = function( seurat, cell.metadata, pcs.use, reduction.use = "pca" ) {
  
  if ( length( cell.metadata$umap.x ) > 0 ) { cell.metadata$umap.x = NULL; cell.metadata$umap.y = NULL }
  
  seurat = RunUMAP( seurat, dims = pcs.use, reduction = reduction.use, verbose = FALSE )
  umap.dt = data.table::as.data.table( seurat@reductions$umap@cell.embeddings, keep.rownames = "cell" )
  colnames( umap.dt )[ 2:3 ] = c( "umap.x", "umap.y" )
  cell.metadata = merge( cell.metadata, umap.dt, by = "cell", all.x = TRUE )
  
  n.replicates = nlevels( cell.metadata$replicate )
  
  x.range = range( cell.metadata$umap.x )
  png( paste0( "umap_with_replicates_", max( pcs.use ), "_", reduction.use, ".png" ), width = 1200, height = 800 )
  cell.metadata[ ( random ), plot( umap.x, umap.y, pch = 20, col = rainbow( n.replicates, s = 0.5, v = 0.9 )[ replicate ], 
                                   axes = FALSE, xlab = "", ylab = "", 
                                   xlim = c( x.range[ 1 ] - ( x.range[ 2 ] - x.range[ 1 ] ) / 2, x.range[ 2 ]  ) ) ]
  legend( "topleft", legend = levels( cell.metadata$replicate ), pch = 15, col = rainbow( n.replicates, s = 0.5, v = 0.9 ), bty = "n" )
  dev.off()
  
  return( cell.metadata )
  
}


#'
#'
plot.cell.features = function( cell.features, cell.metadata ) {
  
  if ( !dir.exists( "feature_plots" ) ) { dir.create( "feature_plots" ) }
  setwd( "feature_plots" )
  if ( !dir.exists( "cell_features" ) ) { dir.create( "cell_features" ) }
  setwd( "cell_features" )
  
  # shared.cols = intersect( colnames( cell.metadata ), colnames( cell.features ) )
  tmp = merge( cell.metadata[ , . ( cell, umap.x, umap.y ) ], cell.features, by = "cell" )
  
  col.classes = sapply( cell.features, class )
  
  for( column.name in colnames( cell.features )[ col.classes %in% c( "integer", "numeric" ) ] ) {
    
    message( "plotting ", column.name )
    
    values = tmp[ , column.name, with = FALSE ]
    variance = var( values )
    
    if ( variance > 1 ) {
      
      message( "   using log transformation" )
      
      loading = create.loading.colors( log10( tmp[[ column.name ]] + 1e-4 ), 
                                       center.at = min( log10( tmp[[ column.name ]] + 1e-4 ) ),
                                       center.col = "lightgrey", neg.col = "royalblue", pos.col = "seagreen" )
      title = paste0( "log10( ", column.name, " )" )
      
    } else {
      
      loading = create.loading.colors( tmp[[ column.name ]], 
                                       center.col = "lightgrey", neg.col = "royalblue", pos.col = "seagreen" )
      title = column.name
      
    }
    
    png( paste0( gsub( ".", "_", column.name, fixed = TRUE ), ".png" ), width = 800, height = 800 )
    cell.metadata[ loading$order, plot( umap.x, umap.y, pch = 20, col = loading$colors[ loading$order ], axes = FALSE, 
                                        main = title, xlab = "", ylab = "" ) ]
    dev.off()
    
  }
  
  setwd( ".." )
  setwd( ".." )
  
}


#'
#'
annotate.cells.using.scpred.model = function( cell.metadata, dge, model.path, selection.threshold = 0.95, use.name = "", 
                                              plot = TRUE ) {
  
  reference = readRDS( model.path )
  
  # input dge should be seurat@assays$RNA@counts or equivalent
  query = CreateSeuratObject( dge )
  query = NormalizeData( query )
  query = scPredict( query, reference, threshold = selection.threshold, recompute_alignment = TRUE )
  
  scpred.data = query@meta.data[ , grepl( "scpred", colnames( query@meta.data ), fixed = TRUE ) ]
  n.models = ncol( scpred.data ) - 3
  next.best = apply( scpred.data[ , 1:n.models ], 1, function( x ) { sort( x, partial = n.models - 1 )[ n.models - 1 ] } )
  
  new.table = data.table::data.table( cell = row.names( scpred.data ), scpred.call = scpred.data$scpred_prediction, 
                          best.model = scpred.data$scpred_no_rejection,
                          max.score = scpred.data$scpred_max, next.best.score = next.best )
  
  new.table = merge( cell.metadata, new.table, by = "cell" )
  print( round( table( new.table$scpred.call ) / nrow( new.table ), 2 ) )
  
  annotation.levels = c( sort( new.table[ , unique( scpred.call[ scpred.call != "unassigned" ] ) ] ), "unassigned" )
  new.table$scpred.call = factor( new.table$scpred.call, levels = annotation.levels )
  
  if ( plot ) {
    
    annotation.palette = c( rainbow( nlevels( new.table$scpred.call ) - 1, s = 0.5, v = 0.9 ), "lightgrey" )
    
    # plot in session
    new.table[ , plot( umap.x, umap.y, pch = 20, col = annotation.palette[ scpred.call ], 
                       axes = FALSE, xlab = "", ylab = "" ) ] 
    legend( "topleft", levels( new.table$scpred.call ), fill = annotation.palette, bty = "n", border = NA )
    scpred.cluster.centers = new.table[ , lapply( .SD, median ), by = scpred.call, .SDcols = umap.x:umap.y ]
    text( scpred.cluster.centers[ , umap.x:umap.y ], labels = scpred.cluster.centers$scpred.call )
    
    # save plot
    png( paste0( "scpred_assignments_", use.name, "_umap.png" ), width = 800, height = 800 )
    new.table[ , plot( umap.x, umap.y, pch = 20, col = annotation.palette[ scpred.call ], 
                       axes = FALSE, xlab = "", ylab = "" ) ] 
    legend( "topleft", levels( new.table$scpred.call ), fill = annotation.palette, bty = "n", border = NA )
    scpred.cluster.centers = new.table[ , lapply( .SD, median ), by = scpred.call, .SDcols = umap.x:umap.y ]
    text( scpred.cluster.centers[ , umap.x:umap.y ], labels = scpred.cluster.centers$scpred.call )
    dev.off()
    
  }
  
  return( new.table[ , .( cell, scpred.call, best.model, max.score, next.best.score ) ] )
  
}


#'
#'
start.sub.clustering = function( n.pcs = 100, n.features = 3000, root.dir ) {
  
  clustering.dir = getwd()
  
  if ( !file.exists( "cell_metadata.RDS" ) ) {
    
    message( "run global clustering first ( run.subclustering = FALSE )" )
    
  } else {
    
    library.paths.and.metadata = fread( "library_paths_and_metadata.txt" )
    global.cell.metadata = readRDS( "cell_metadata.RDS" )
    cell.features = readRDS( "cell_features.RDS" )
    
    # plot global clusters
    global.annotation.centers = global.cell.metadata[ , lapply( .SD, median ), by = .( annotation ), .SDcols = umap.x:umap.y ]
    global.annotation.cell.counts = global.cell.metadata[ , .N, by = .( annotation ) ]
    global.n.annotations = nlevels( global.cell.metadata$annotation )
    global.cell.metadata[ , plot( umap.x, umap.y, pch = 20, 
                                  col = rainbow( global.n.annotations, s = 0.25, v = 0.9 )[ annotation ], 
                                  axes = FALSE, xlab = "", ylab = "" ) ]  
    text( global.annotation.centers[ , umap.x:umap.y ], 
          labels = paste0( global.annotation.centers$annotation, "\n", "n = ", global.annotation.cell.counts$N ) )
    
    message( "enter *unique* annotations to subcluster ( example: astro )" )
    
    print( levels( global.cell.metadata$annotation ) )
    annotation.use = readline()
    message( "using ", annotation.use )
    
    if ( !dir.exists( "sub_clustering" ) ) { dir.create( "sub_clustering" ) }
    setwd( "sub_clustering" )
    if ( !dir.exists( annotation.use ) ) { dir.create( annotation.use ) }
    subclustering.dir = paste0( clustering.dir, "/sub_clustering/", annotation.use )
    setwd( subclustering.dir )
    
    if ( file.exists( "seurat.RDS" ) ) {
      
      message( "start a new clustering? ( y / n )" )
      if ( readline() == "y" ) { new.clustering = TRUE } else { new.clustering = FALSE }
      
    } else { new.clustering = TRUE }
    
    if ( !new.clustering ) {
      
      message( "loading existing seurat" )
      seurat = readRDS( "seurat.RDS" )
      
    } else {
      
      dge = readRDS( paste0( clustering.dir, "/seurat.RDS"  ) )@assays$RNA@counts
      dge = dge[ , global.cell.metadata[ annotation %in% annotation.use, cell ] ]
      
      if( file.exists( paste0( root.dir, "/excluded_cell_barcodes.txt" ) ) ) {
        
        message( "remove excluded cell barcodes? ( y / n )" )
        if ( readline() == "y" ) { 
          
          excluded.barcodes = scan( paste0( root.dir, "/excluded_cell_barcodes.txt" ), what = "character" )
          dge = dge[ , !( colnames( dge ) %in% excluded.barcodes ) ]
          
        }
        
      }
      
      gene.sums = rowSums( dge )
      dge = dge[ gene.sums > 0, ]
      
      message( "creating seurat object" )
      seurat = CreateSeuratObject( counts = dge )
      
      message( "running sctransform" )
      message( "using ", n.features, " features" ) # default is 3000 ( for SCTransform )
      if ( n.features == "all" ) {
        
        seurat = SCTransform( seurat, method = "glmGamPoi",
                              variable.features.n = nrow( seurat@assays$RNA@counts ), verbose = FALSE )
        
      } else {
        
        seurat = SCTransform( seurat, method = "glmGamPoi", 
                              variable.features.n = n.features, verbose = FALSE )
        
      }
      
      message( "generating ", n.pcs, " pcs" )
      seurat = RunPCA( seurat, npcs = n.pcs, verbose = FALSE )
      
      message( "saving seurat object" )
      saveRDS( seurat, "seurat.RDS" )
      
      if ( !dir.exists( "dges" ) ) { dir.create( "dges" ) }
      setwd( "dges" )
      for ( i in 1:nrow( library.paths.and.metadata ) ) {
        
        library = library.paths.and.metadata[ i, my.name ]
        sub.dge = dge[ , grep( library, colnames( dge ), fixed = TRUE ) ]
        out.file = paste0( library, "_sparse_dge.RDS" )
        saveRDS( sub.dge, out.file )
        library.paths.and.metadata[ i, sparse.dge := paste0( getwd(), "/", out.file ) ]
        
      }
      setwd( ".." )
      
      write.table( library.paths.and.metadata, "library_paths_and_metadata.txt", row.names = FALSE, quote = FALSE )
      
      return( seurat )
      
    }
    
  }
  
}


#'
#'
update.excluded.cell.barcodes = function( cells.exclude, root.dir ) {
  
  if ( file.exists( paste0( root.dir, "/excluded_cell_barcodes.txt" ) ) ) {
    
    prev.excluded = scan( paste0( root.dir, "/excluded_cell_barcodes.txt" ), what = character() )
    cells.exclude = unique( c( prev.excluded, cells.exclude ) )
    
    message( "updating excluded_cell_barcodes.txt" )
    write.table( cells.exclude, paste0( root.dir, "/excluded_cell_barcodes.txt" ), 
                 row.names = FALSE, col.names = FALSE, quote = FALSE )
    
  } else {
    
    message( "writing new excluded_cell_barcodes.txt" )
    write.table( cells.exclude, paste0( root.dir, "/excluded_cell_barcodes.txt" ), 
                 row.names = FALSE, col.names = FALSE, quote = FALSE )
    
  }
  
}

#' Test if string ends with given suffix
#'
#' @param theString string to be tested
#' @param theExt suffix to be tested for
#' @return TRUE if theString ends with theExt
#' @export
strEndsWith<-function(theString,theExt) {
  return(substring(theString,1+nchar(theString)-nchar(theExt))==theExt)
}

#' Creates a file connection with the given open mode.
#' @param file If ends with ".gz", a gzfile() is created; else a regular file() connection.
#' @param open mode in which file is opened.  Default: "rb"
#' @export
open_conn = function(file, open="") {
  if (strEndsWith(file, ".gz")) {
    # work around bug in gzfile:
    # https://stackoverflow.com/questions/45665496/how-would-one-readlines-from-a-gzip-file-in-r
    if (nchar(open) == 0) {
      open = "rb"
    }
    if (!strEndsWith(open, "b")) {
      open=paste0(open, "b")
    }
    return(gzcon(file(file, open=open)))
  } else {
    return(file(file, open=open))
  }
}

#' Helper function for writing a possibly gzipped data.table/data.frame
#' @param outPath path to output file
#' @param x data.table/data.frame to write
#' @param col.names whether to write column names
#' @param row.names whether to write row names
#' @param quote whether to quote strings
#' @param sep column separator
#' @param ... additional arguments to write.table
#' @import utils
write_table_helper<-function(outPath, x, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t", ...) {
  outconn <- open_conn(outPath, "w")
  utils::write.table(x, outconn, col.names=col.names, row.names=row.names, quote=quote, sep=sep)
  close(outconn)
}