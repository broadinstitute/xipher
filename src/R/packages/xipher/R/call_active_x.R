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

#' phase x chromosome SNPs
#' @param outErrorRatePath path to write error rate table
#' @param outXCallsPath path to write active-X calls
#' @param preparedDacPath file produced by prepareDacClp()
#' @param phasePath file produced by phaseClp()
#' 
#' TODO: Steve: document the rest of the parameters
#' @export 
#' @import data.table
callActiveXClp <- function(outErrorRatePath,
                     outXCallsPath,
                     preparedDacPath,
                     phasePath,
                     max.error.rate=0.25,
                     min.cells.final=4,
                     calling.purity.min=0.85,
                     likelihood.min=0.9) {
  dac <- data.table::fread(preparedDacPath)
  flip.record <- data.table::fread(phasePath)


  # get error rates and x calls -------------------------------------------------------------------------------------
  
  # print libraries
  message( table( sapply( unique( dac$cell ), get.cell.barcode.prefix ) ) )
  
  phased.positions = unlist( strsplit( flip.record[ which.max( pct.total ), contains ], ",", fixed = TRUE ) )
  message( length( phased.positions ), " phased postions" )
  flip.positions = flip.record[ pos %in% phased.positions & flip, pos ]
  message( length( flip.positions ), " flip postions ( ", round( length( flip.positions ) / length( phased.positions ), 2 ), " )" )
  
  if ( length( phased.positions ) < 1 ) { 
    
    stop( "insufficient initial phasing ( quitting )" )
  }
  
  x.calls = call.active.x( dac, phased.positions, flip.positions )
  # saveRDS( x.calls, "x_calls_initial.RDS" )
  
  
  position.stats = generate.error.rates.and.update.phasing( dac, phased.positions, flip.positions, min.cell = min.cells.final )
  
  phased.positions = position.stats[ ( ( error.rate <= max.error.rate ) | ( error.rate >= ( 1 - max.error.rate ) ) ) & total >= min.cells.final, pos ]
  message( length( phased.positions ), " phased postions" )
  flip.positions = position.stats[ pos %in% phased.positions & error.rate > 0.5, pos ]
  message( length( flip.positions ), " flip postions ( ", round( length( flip.positions ) / length( phased.positions ), 2 ), " )" )
  position.uncertainties = position.stats[ , .( pos, gene, uncertainty ) ]
  
  x.calls = call.active.x( dac, phased.positions, flip.positions, position.uncertainties = position.uncertainties,
                           calling.purity.min = calling.purity.min, likelihood.min = likelihood.min )
  
  write_table_helper(outErrorRatePath, position.stats)
  write_table_helper(outXCallsPath, x.calls)
  
  message( "done!" )
}