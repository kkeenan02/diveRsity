# NEWS


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

v1.9.8.4
----------

Speed up gpSampler function by keeping data as charater strings, since this is the only data structure required.

Also adds a bug fix to `fastScan` function. Some files (e.g. EASYPOP output .gen files) have "\r\n" file endings. when replacing "\r" with "\n" in these files, lines ended with "\n\n" meaning blank lines between each real line. The code below fixes this:

```r
fastScan <- function(fname) {
  s <- file.info(fname)$size
  buf <- readChar(fname, s, useBytes = TRUE)
  # replace Mac encoded line endings
  if(grep("\r", buf) == 1L){
    buf <- gsub("\r", "\n", buf)
    buf <- gsub("\n\n", "\n", buf)
  }
  return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
}
```


v1.9.8.5
---------

Fixes a memory leak by managing parallel cluster in a more appropriate way.

v1.9.8.6
---------

Added a faster C++ verison of tabMerge function

v1.9.8.7
--------

Added a faster C++ version of tabMergePw for pairwise hsum calculation.

v1.9.8.8
--------

Fixes a bug in pairwise D (Jost, 2008) calculation.

v1.9.8.9
--------
Fixes a bug in pairwise in calculations in `inFunc`. The following line of code resulted in a vector for pairwise populations that were fixed for the same allele, when a matrix was expected by downstream code.

```r
y <- y[rowSums(y) != 0,]
```

By deleting this line, the code downstream generates `NAN` values, which are now summed using `na.rm = TRUE`.

v1.9.9.0
--------

Added updates for divMigrate-online.

v1.9.9.40
---------
Exact tests of HWE added in divBasic function

v1.9.9.43
---------
New chiCalc function with new output structure and pairwise calculations. P Values are now calculated using `fisher.test` rather than `chisq.test`.

v1.9.9.44
---------
Users can now choose whether to write the genepop object generated by snp2gen to file, or just return it to the workspace for use in other diveRsity functions (suggestion from Mark Ravinet).

v1.9.9.45
---------
`divBasic` now returns `mainTab` as a list of dataframes, rather than a character matrix. The function needs a lot of work interms of speed when calculating $F_{IS}$ and $A_{R}$.

v1.9.9.49
---------
`diffCalc` crashed when some loci had 100% missing data. The problem function was `glbWCcpp`. This version deals with this bug. Extensive testing of the fix is still required and will be carried out for v1.9.9.50

v1.9.9.51
---------
The `corplot` function as been revamped to make use of ggplot2 graphics.

v1.9.9.53
---------
Fixed a bug in diffCalc pairwise Fst calculations. This stemmed from the fix added in v1.9.9.49.

#### old code

```r
hsum <- lapply(hsum, function(x){
   x[!(names(x) == "NA")]
})
```

#### New code

```r
hsum <- lapply(hsum, function(x){
   x <- lapply(x, function(y){
      y <- y[!(names(y) == "NA")]
      return(y)
   })
})
```