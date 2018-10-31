#!/bin/bash

rm -rf macos
mkdir macos

LIB=/usr/local/lib
LIBGFORTRAN=${LIB}/libgfortran.a
LIBQUADMATH=${LIB}/libquadmath.a
LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

g++ *.cpp -o macos/rtm-gwas-gsc -O2 -std=c++11 -llapack -lrefblas $LIBGFORTRAN $LIBQUADMATH $LIBGCC
