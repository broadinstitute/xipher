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
#' @param outPhasePath path to write x SNP phase table
#' @param preparedDacPath file produced by prepareDacClp()
#' @param X_contig_name name of the X contig in the reference genome
#' 
#' TODO: Steve: document the rest of the parameters
#' @export 
#' @import data.table
phaseClp <- function(outPhasePath,
                     preparedDacPath, 
                     min.cells=2,
                     max.xi.proportion=0.1,
                     max.gen.units=1000,
                     X_contig_name=default_X_contig_name) {
  dac <- data.table::fread(preparedDacPath)
  flip.record = phase.x.variants( dac, min.cells = min.cells, max.xi.proportion = max.xi.proportion, 
                                  time.limit = 2, max.gen.units = max.gen.units, X_contig_name=X_contig_name )
  write_table_helper(outPhasePath, flip.record)
}
