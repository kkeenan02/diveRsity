################################################################################
# gpSampler                                                                    #
################################################################################
#' @export
gpSampler <- function(infile = NULL, samp_size = 10, outfile = NULL){
  fastScan <- function(infile) {
    s <- file.info(infile)$size
    buf <- readChar(infile, s, useBytes = TRUE)
    # replace Mac encoded line endings
    if(length(grep("\r", buf)) != 0L){
      buf <- gsub("\r", "\n", buf)
      buf <- gsub("\n\n", "\n", buf)
    }
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  dat <- fastScan(infile)
  # strip whitespace from the beginning an end of lines
  dat <- sapply(dat, function(x){
    sub("^\\s+", "", x)
  })
  dat <- sapply(dat, function(x){
    return(sub("\\s+$", "", x))
  })
  names(dat) <- NULL
  popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
  if(popLoc[1] == 3){
    if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
      locs <- strsplit(dat[2], split = "\\s+")[[1]]
      if(length(locs) == 1){
        locs <- strsplit(dat[2], split = ",")[[1]]
      }
      locs <- as.character(sapply(locs, function(x){
        x <- strsplit(x, split = "")[[1]]
        if(is.element(",", x)){
          x <- x[-(which(x == ","))]
        }
        return(paste(x, collapse = ""))
      }))
      dat <- c(dat[1], locs, dat[-(1:2)])
    }
  } else {
    locs <- as.character(dat[2:(popLoc[1]-1)])
  }
  # strip whitespace from locus names
  locs <- as.character(sapply(locs, function(x){
    return(strsplit(x, split = "\\s+")[[1]][1])
  }))
  # npops
  popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
  npops <- length(popLoc)
  no_col <- length(locs)+1
  nloci <- length(locs)
  # get genotypes
  strt <- popLoc + 1
  ends <- c(popLoc[-1] - 1, length(dat))
  # generate a list of pop data
  popDat <- function(dat, strt, nd){
    return(dat[strt:nd])
  }
  pop_dat <- mapply(popDat, strt = strt, nd = ends, MoreArgs = list(dat = dat),
                    SIMPLIFY = FALSE)
  # function for sampling strings
  smplR <- function(x, n){
    id <- sample(length(x), n, replace = FALSE)
    return(x[id])
  }
  # resample data
  if(length(samp_size) == 1L){
    out <- lapply(pop_dat, smplR, n = samp_size)
  } else {
    out <- mapply(smplR, x = pop_dat, n = samp_size, SIMPLIFY = FALSE)
  }
  # construct the new output
  out <- lapply(out, function(x){
    c("POP", x)
  })
  out <- do.call("c", out)
  out <- c(paste("gpSampler n=", samp_size[1], sep = ""), locs, out)
  out <- paste(out, collapse = "\n")
  writeLines(out, outfile)
}
################################################################################
# gpSampler                                                                    #
################################################################################