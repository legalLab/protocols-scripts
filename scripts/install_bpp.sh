#!/usr/bin/env sh
# script to install Yang's BPP Software
# works on Ubuntu Linux 16.04 (April 2016)

# install deps
sudo apt-get install gcc

# download latest version
wget http://abacus.gene.ucl.ac.uk/software/bpp3.2a.tgz

# untar
tar -zxvf bpp3.2a.tgz

# cd 
cd bpp3.2a

# remove windows executables
rm *.exe

# compile BPP
gcc -o bpp -O3 bpp.c tools.c -lm