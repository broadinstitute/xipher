#!/bin/bash
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

# Run from the package directory.
# roxygen2 must be installed in the R environment.

set -e

PROGNAME=`basename $0`
function usage () {
    cat >&2 <<EOF
$*

USAGE: $progname [-f]

Roxygenize, build and check package in current directory.

Options:

-f: Fast: Don't R CMD check
EOF
}


fast=0
while getopts ":f" options; do
  case $options in
    f ) fast=1;;
    h ) usage;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;

  esac
done
shift $(($OPTIND - 1))


if [ -f NAMESPACE ]
then rm NAMESPACE
fi
if [ -d man ]
then rm -r man
fi
Rscript -e 'roxygen2::roxygenise()' -e 'warnings()'
R CMD build .

if (( $fast != 1 ))
then R CMD check *.tar.gz
fi


