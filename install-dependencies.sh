#!/bin/bash
# Script to install external dependencies for datelife R package

# install pathd8
# get the package
wget https://www2.math.su.se/PATHd8/PATHd8.zip -P /tmp
# unzip it
unzip /tmp/PATHd8 -d /usr/local/bin/
# compile
cc PATHd8.c -O3 -lm -o PATHd8
