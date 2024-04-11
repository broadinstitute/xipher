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
get.corresponding.labels = function( labels, from, to ) {
  
  from = as.character( from )
  labels = as.character( labels )
  names( to ) = from
  new.labels = to[ labels ]
  names( new.labels ) = NULL
  
  return( new.labels )
  
}









# x analyses ------------------------------------------------------------------------------------------------------

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
#' @import data.table stringr
phase.x.variants = function( unphased.dac, gen.unit.use = "pos", xi.genes = NULL, genes.exclude = NULL, 
                             cells.exclude = NULL, min.cells = 10, max.xi.proportion = 0.2, time.limit = 48,
                             max.gen.units = 1000,
                             X_contig_name=default_X_contig_name ) {
  
  first.start = Sys.time()
  current.time = Sys.time()
  
  copy.unphased.dac = data.table::copy( unphased.dac )
  # copy.unphased.dac = copy( dac )
  copy.unphased.dac[ chr == X_contig_name ]
  
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
#' @import binom
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
#' @import data.table
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
#' @import data.table
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




# clustering ------------------------------------------------------------------------------------------------------

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

default_X_contig_name <- "X"