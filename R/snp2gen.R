#' snp2gen function
#'
#' Convert genotypes in a snp matrix to genepop genotype files.
#'
#' Kevin Keenan, QUB, 2014
snp2gen <- function(infile = NULL, prefix_length = 2){
  if(is.null(infile)){
    stop("Please provide an input file!")
  }
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    # replace Mac encoded line endings
    if(length(grep("\r", buf)) != 0L){
      buf <- gsub("\r", "\n", buf)
      buf <- gsub("\n\n", "\n", buf)
    }
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  if(is.list(infile)){
    infile <- as.matrix(infile)
    dat <- sapply(1:nrow(infile), function(i){
      out <- paste(infile[i,], collapse = "\t")
      return(out)
    })
    dat <- c(paste(colnames(infile), collapse = "\t"), dat)
    pre <- "snp2gen"
  } else {
    # deal with relative paths
    pre <- strsplit(infile, split = "\\.")[[1]]
    if(length(pre) == 3L){
      pre <- paste(getwd(), pre[2], sep = "")
    } else {
      pre <- paste(getwd(), paste(pre[-length(pre)], collapse = ""), sep = "")
    }
    dat <- fastScan(infile)
    # deal with empty last lines in files
    if(dat[length(dat)] == ""){
      dat <- dat[-length(dat)]
    }
  }
  inds <- strsplit(dat[1], split = "\\s+")[[1]][-1]
  splitNames <- lapply(inds, function(x){
    return(strsplit(x, split = "")[[1]])
  })
  pop_pre <- sapply(splitNames, function(x){
    return(paste(x[1:prefix_length], collapse = ""))
  })
  prefixes <- unique(pop_pre)
  pop_idx <- lapply(prefixes, function(x){
    which(x == pop_pre)
  })
  npops <- length(prefixes)
  pop_sizes <- sapply(pop_idx, length)
  # organise genotypes into matrix
  genos <- t(sapply(dat[-1], function(x){
    return(strsplit(x, split = "\\s+")[[1]])
  }))
  dimnames(genos) <- list(NULL, NULL)
  # extract snp names
  locs <- as.vector(genos[,1])
  nloci <- length(locs)
  genos <- genos[,-1]
  genos <- gsub("A", "01", toupper(genos))
  genos <- gsub("C", "02", toupper(genos))
  genos <- gsub("G", "03", toupper(genos))
  genos <- gsub("T", "04", toupper(genos))
  genos <- gsub("-", "00", toupper(genos))
  # extract populations
  genos <- lapply(pop_idx, function(x){
    if(length(x) == 1L){
      paste(genos[,x], collapse = "\t")
    } else {
      apply(genos[,x], 2, paste, collapse = "\t")
    }
  })
  # get ind names
  indNames <- lapply(pop_idx, function(x){
    return(paste(inds[x], " ,", "\t", sep = ""))
  })
  genos <- mapply(paste0, indNames, genos, SIMPLIFY = FALSE)
  genos <- unlist(sapply(genos, function(x){
    c("POP", x)
  }))
  genos <- c("snp2gen-converted", paste(locs, collapse = ", "),
             genos)
  genos <- paste0(genos, collapse = "\n")
  writeLines(genos, "snp2gen_converted.gen")
}