#!/bin/sh
R CMD build Rgbp
R CMD check --no-examples Rgbp_*.tar.gz
