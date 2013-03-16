#!/bin/sh
R CMD build Rgbp
R CMD check --as-cran --no-examples Rgbp_*.tar.gz
