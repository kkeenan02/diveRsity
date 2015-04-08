#' snp2gen function
#' 
#' Convert genotypes in a snp matrix to genepop genotype files.
#' 
#' Kevin Keenan, QUB, 2014
snp2gen <- function(infile = NULL, prefix_length = 2, write = FALSE){
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
  genos <- strsplit(genos, split = "")
  al1 <- matrix(sapply(genos, "[", 1), ncol = sum(pop_sizes),
                nrow = nloci, byrow = FALSE)
  al2 <- matrix(sapply(genos, "[", 2), ncol = sum(pop_sizes),
                nrow = nloci, byrow = FALSE)
  # geno array [ind, loc, allele]
  genos <- array(NA, dim = c(nrow(al1), ncol(al1), 2))
  genos[,,1] <- al1
  genos[,,2] <- al2
  # replace nucleotides with genepop ids
  genos[toupper(genos) == "A"] <- "01"
  genos[toupper(genos) == "C"] <- "02"
  genos[toupper(genos) == "G"] <- "03"
  genos[toupper(genos) == "T"] <- "04"
  genos[genos == "-"] <- "00"
  # extract populations
  pop_list <- lapply(pop_idx, function(x){
    return(genos[,x,])
  })
  # convert pop_list to genepop format
  pop_list <- lapply(pop_list, function(x){
    out <- apply(x, 2, function(y){
      return(paste(y[,1], y[,2], sep = ""))
    })
    return(t(out))
  })
  pop_list <- lapply(pop_list, function(x){
    return(rbind(rep(NA, ncol(x)), x))
  })
  # get ind names
  indNames <- lapply(pop_idx, function(x){
    return(c("pop", paste(inds[x], " ,", sep = "")))
  })
  pop_list <- cbind(do.call("c", indNames),
                    do.call("rbind", pop_list))
  pop_list <- rbind(c(paste(pre, "-converted", sep = ""), 
                      rep(NA, nloci)),
                    c(c(paste(locs[1:(nloci-1)], ",", sep = ""), 
                        locs[nloci]), NA), pop_list)
  pop_list[is.na(pop_list)] <- ""
  if(write){
    # write the results
    fl <- file(paste(pre, "_converted.gen", sep = ""), "w")
    for(i in 1:nrow(pop_list)){
      out <- pop_list[i,]
      if(all(out[-1] == "\t")){
        out <- out[1]
      }
      cat(out, "\n", file = fl, sep = "\t")
    }
    close(fl)
    z <- gc()
  }
  return(as.data.frame(pop_list))
}