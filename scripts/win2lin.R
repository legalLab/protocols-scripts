#!/usr/bin/env Rscript

# function to convert a Windows file path to Linux if on a Linux machine, or keep as Windows if on Windows
# only works for relative file paths (e.g. 'C:\\path\\to' will not work, but '..\\path\\to' or 'path\\to' will work) 
# Rupert Collins August 2018

win2lin <- function(path){
    if(.Platform$OS.type=="windows"){
    return(path)}
    else{
    path.lin <- gsub("\\\\", "/", path)
    return(path.lin)}
}

# example 
# win2lin("..\\file\\path\\to")