#!/usr/bin/env sh

# positional parameters:
# $0 = basename of program!
# $1 = name of your PyRAD .loci file

# clean up the .loci file by replacing the // lines with a single hash, and replacing whitespace with newlines (i.e. convert to fasta) 
sed -e 's/^\/\/.*/#/g' -e 's/\s\+/\n/g' $1 > cleaned.loci

# delete the last line (as there is no info below the split)
sed -i '$d' cleaned.loci

# split the file into component loci (fasta format)
csplit -f locus -n 4 cleaned.loci "/#/" "{*}"

# run the R script to add the missing data
Rscript pyrad2gphocs.R

# replace the hash/newline after the locus number
sed -i ':a;N;$!ba;s/#\n//g' gphocs.loci

# clean up
rm cleaned.loci
rm locus*