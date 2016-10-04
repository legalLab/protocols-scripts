# pyrad2gphocs
Scripts to convert PyRAD output into G-PhoCS input

### Description
These scripts convert PyRAD output files in .loci format into input files suitable for use in G-PhoCS (<http://compgen.cshl.edu/GPhoCS/>). It's a little bit hacky, and probably would be better written from scratch in Python.

Here an R script is called from within a shell script, and the R script needs to be modified to include the path to the "ape" package (mine is in a non-standard location)

Need to also provide a file with the names of all the individuals. This file needs to be named "samples_names.txt". Either change the name of your file to match this name, or modify the R script to change the file name there.

Don't need to change your input file to "example.loci", as any file name can be parsed.

Outputs a file named "gphocs.loci".

### To run
```
./pyrad2gphocs.sh example.loci
```

### Note
Probably redundant, as it looks like gphocs output functionality is available in the more recent versions of PyRAD (e.g. 3.0.6).