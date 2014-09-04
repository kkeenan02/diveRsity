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

```
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

```
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