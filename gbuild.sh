#!/bin/sh
R CMD build Rgbp
R CMD check --as-cran Rgbp_*.tar.gz
