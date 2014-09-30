---
title: "NEWS"
author: "Kevin Keenan"
date: "17/04/2014"
output: html_document
---

diveRsity update notes
======================

v1.9.5 up
---------------

Bug fix in internal function, `rgp`, allowing horizontal locus names in input file.


v1.9.6
--------------

Meirmans & Hedrick, (2011) $G''_{ST}$ added to `diffCalc`. The statistic is calculated for all available levels of genetic integration.

v1.9.7
-----------

Alcala et al, 2014 Nm estimation is added to `divMigrate` as a option to calculation relative migration. Gst is also added.

v1.9.8
------------
Bug fixed in 'rgp' function to correctly deal with leading whitespace in sample
names.

v1.9.8.1
------------
Bug fixed in check function within rgp to ensure allele names were recorded correctly.

### Original

```r
 check <- function(args){
    #args <- list(...)
    npops <- length(args)
    nchars <- mean(nchar(names(args[[1]])))
    pad <- paste("%0", nchars, "g", sep = "")
    rnames <- sprintf(pad, 
                      unique(sort(as.numeric(unlist(lapply(args,
                                                           names))))))
    
    out <- matrix(0, nrow = length(rnames), ncol = npops)
    rownames(out) <- as.character(rnames)
    for(i in 1:npops){
      out[match(names(args[[i]]), rownames(out)),i] <- as.numeric(args[[i]])
    }
    return(out)
  }
```

### Fixed code

```r
 check <- function(args, gp){
    #args <- list(...)
    npops <- length(args)
    pad <- paste("%0", gp, "g", sep = "")
    rnames <- sprintf(pad, 
                      unique(sort(as.numeric(unlist(lapply(args,
                                                           names))))))
    
    out <- matrix(0, nrow = length(rnames), ncol = npops)
    rownames(out) <- as.character(rnames)
    for(i in 1:npops){
      out[match(names(args[[i]]), rownames(out)),i] <- as.numeric(args[[i]])
    }
    return(out)
  }
```

v1.9.8.2
-----------
Multilocus gst and nm values for directional migration were not behaving as expected, when
dealing with 100% missing data. The below code was added to indicate loci with 100% missing data:

```r
  # fix allele frequencies
  dat$af <- lapply(dat$af, function(x){
    cs <- colSums(x)
    x[,cs == 0] <- NA
    return(x)
  })
```
These loci can now be identified in 'pwHt' and NA returned, instead of NaN, as was previously occuring. This previous behaviour resulted in NaN values derived from 100% missing data being treated in the same way as NaN values drived from two population samples not sharing any alleles, which was obviously incorrect.

v1.9.8.3
----------
Fixed a bug in the fastScan function which did not recognise "\r" as a line ending.

#### Original

```r
 fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
```

##### Fixed
```r
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    # replace Mac encoded line endings
    if(grep("\r", buf) == 1L){
      buf <- gsub("\r", "\n", buf)
    }
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
```

This new version of the function checks to see if "\r" is a character present in the file, then replaces these with "\n" so that the downstream functions work normally.