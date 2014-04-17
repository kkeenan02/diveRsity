#' @import ggplot2 shiny qgraph plotrix
NULL
################################################################################
################################################################################
##                              diveRsity v1.9.5                              ##  
##                            by Kevin Keenan QUB                             ##  
##            An R package for the calculation of differentiation             ##
##              statistics and locus informativeness statistics               ##  
##                V 1.2.0 and up allows parallel computations                 ##  
##                            GPL3 Kevin Keenan 2014                          ##  
################################################################################
################################################################################

#' @export
divPart <- function(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
                    WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE, 
                    bootstraps = 0, plot = FALSE, parallel = FALSE){
  ############################ Argument definitions ############################
  .Deprecated(new = "fastDivPart", msg = "This function name is no longer in use. Please use 'fastDivPart' instead, \nSee ?fastDivPart for usage details.", 
              old = "divPart")
}
################################################################################
# div.part: deprecated
################################################################################
#' @export
# div.part, a wrapper function for the calculation of differentiation stats.
div.part<-function(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
                   WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE, 
                   bootstraps = 0, plot = FALSE, parallel = FALSE){
  ############################ Argument definitions ############################
  .Deprecated(new = "divPart", msg = "This function name is no longer in use. Please use 'divPart' instead, \nSee ?divPart for usage details.", 
              old = "div.part")
}
################################################################################
# End div.part
################################################################################
#
#
#
#
#
#
#
################################################################################
# readGenepopX, a function for the generation of basic population parameters   #
################################################################################
readGenepopX <- function (x) {
  #gp=x$gp
  infile=x$infile
  bootstrap=x$bootstrap
  locs=x$locs
  data1 <- fileReader(infile)
  if(is.null(x$gp)){
    rownames(data1) <- NULL
    data1 <- as.matrix(data1)
    # determine genepop format
    p1 <- which(toupper(data1[,1]) == "POP")[1] + 1
    gp <- as.numeric(names(sort(-table(sapply(data1[p1,(-1)], nchar)/2)))[1])
    data1 <- as.data.frame(data1)
  } else {
    gp <- x$gp
  }
  
  if(gp == 3){
    data1[data1==0]<-NA
    data1[data1=="999999"]<-NA
    data1[data1=="000000"]<-NA
    data1[data1=="NANA"]<-NA
  } else if(gp == 2){
    data1[data1==0]<-NA
    data1[data1=="9999"]<-NA
    data1[data1=="0000"]<-NA
    data1[data1=="NA"]<-NA
  }
  raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
  pop_weights<- 1/pop_sizes
  
  n_harmonic<-npops/sum(pop_weights)
  
  N<-pop_sizes
  
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list<-list()
  for (i in 1:npops){
    pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                   2:(nloci+1)])
  }
  # check if all populations have at least some data at loci
  extCheck <- sapply(1:length(pop_list), function(i){
    sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
  })
  if (sum(extCheck) > 0){
    npops <- npops - sum(extCheck)
    pop_list <- pop_list[-(which(extCheck == TRUE))]
    pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
    pop_names <- pop_names[-(which(extCheck == TRUE))]
    pop_weights <- pop_weights[-(which(extCheck == TRUE))]
    N <- N[-(which(extCheck == TRUE))]
    #raw_data fix
    noPop <- which(extCheck == TRUE)
    indexer <- lapply(noPop, function(i){
      (pop_pos[i] + 1):(pop_pos[(i+1)])
    })
    indexer <- unlist(indexer)
    raw_data <- raw_data[-(indexer), ]    
  }  
  if (gp==3) {
    plMake<-function(x){
      out <- matrix(sprintf("%06g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "0000NA"] <- "    NA"
      }
      return(out)
    }
  } else if (gp==2) {
    plMake<-function(x){
      out <- matrix(sprintf("%04g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "00NA"] <- "  NA"
      }
      return(out)
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  
  
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  
  
  if(bootstrap == T){
    bs<-function(x){
      return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
    }
    pop_list<-lapply(pop_list, bs)
  }  
  
  ###vectorize loci_pop_sizes#####################################################
  
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }#
  pre_loci_pop_sizes<-lapply(pop_list,lps)#
  pls<-matrix(ncol=nloci,nrow=npops)#
  for(i in 1:length(pre_loci_pop_sizes)){#
    pls[i,]<-pre_loci_pop_sizes[[i]]#
  }#
  #convert pls to loci_pop_sizes format
  loci_pop_sizes<-split(pls,col(pls))
  
  
  #vectorized loci_pop_weights##################################################
  
  pre_loc_weights<- 1/pls
  loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
  loci_harm_N<-npops/colSums(pre_loc_weights)
  
  #end vectorized loci_pop_weights##############################################
  
  ###vectorize pop_alleles########################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles####################################################
  
  #vectorize allele_names#######################################################
  
  alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res<-list()
    for(i in 1:ncol(x[[1]])){
      res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    }
    return(res)
  }
  
  allele_names<-lapply(pop_alleles,alln)
  
  
  loci_combi<-allele_names[[1]]
  for(j in 1:nloci){
    for(i in 2:npops){
      loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
    }
  }
  
  #all_alleles vectorized#######################################################
  
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  
  #end all_alleles vectorized###################################################
  
  aa<-all_alleles
  aa<-lapply(aa, FUN=`list`, npops)
  afMatrix<-function(x){
    np<-x[[2]]
    z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
    rownames(z)<-x[[1]]
    return(z)
  }
  allele_freq<-lapply(aa,afMatrix)
  
  
  #combine pop_alleles
  parbind<-function(x){
    rbind(x[[1]],x[[2]])
  }
  pa1<-lapply(pop_alleles, parbind)
  #create a function to tabulate the occurance of each allele
  afTab<-function(x){
    lapply(1:ncol(x), function(i){
      return(table(x[,i]))
    })
  }
  actab<-lapply(pa1, afTab)
  
  afs<-function(x){
    afsint<-function(y){
      length(na.omit(y))/2
    }
    apply(x,2,afsint)
  }
  indtyppop<-lapply(pa1,afs)
  #calculate allele frequencies
  afCalcpop<-lapply(1:length(actab), function(x){
    lapply(1:length(actab[[x]]),function(y){
      actab[[x]][[y]]/(indtyppop[[x]][y]*2)
    })
  })
  #assign allele freqs to frequency matrices
  obs_count<-allele_freq
  for(i in 1:npops){
    for(j in 1:nloci){
      allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
      obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
    }
  }
  
  
  
  indtyp<-list()
  for(i in 1:nloci){
    indtyp[[i]]<-vector()
  }
  for(i in 1:npops){
    for(j in 1:nloci){
      indtyp[[j]][i]<-indtyppop[[i]][j]
    }
  }
  
  if(bootstrap==T){
    ind_vectors<-list()
    for(i in 1:npops){
      ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
    }
    
    
    pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                     ncol=(nloci+1))
    pre_data[1,]<-c("Title",rep("\t",nloci))
    for(i in 2:(nloci+1)){
      pre_data[i,1]<-loci_names[(i-1)]
    }
    pop_data<-list()
    for(i in 1:npops){
      pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                  cbind(ind_vectors[[i]],pop_list[[i]])),
                            ncol=(nloci+1))
    }
    bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
    for(i in 2:npops){
      bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
    }
    bs_data_file<-data.frame(bs_data_file)
  }
  nalleles<-vector()
  for(i in 1:nloci){
    nalleles[i]<- nrow(allele_freq[[i]])
  }
  ##############################################################################
  if(bootstrap==T){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         locs=locs,
         bs_file=bs_data_file,
         obs_allele_num=obs_count)
  } else if(bootstrap==F){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         locs=locs,
         obs_allele_num=obs_count)
  }
}
################################################################################
# readGenepopX end                                                             #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# plotter, a function to create interactive plots of results from divPart      #
################################################################################
plotter<-function(x,img="1200x600"){
  x=x
  spot.radius=5
  jjj<-x
  require("sendplot")
  fl_ext<-c(".tif","Dot.png","Dot.tif")
  bs_res<-list()
  lso123<-list()
  accDat<-list()
  sp.header<-list()
  pw_res<-list()
  pwso<-list()
  pw<-list()
  plot_data321<-list()
  if(jjj$plot_loci==TRUE && jjj$plot_pw==FALSE){
    bs_res<<-jjj$bs_res
    lso123<<-jjj$lso123
    accDat<<-jjj$accDat
    if(length(lso123) > 150){
      image.size <- "2400x1200"
    } else {
      image.size=img
    }
    #Gst_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[1]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[1]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[1]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    #clean up
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[1]],"_locus_stat_",fl_ext,sep=""))
    #G'st_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[2]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[2]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[2]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[2]],"_locus_stat_",fl_ext,sep=""))
    #Djost_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[3]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[3]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[3]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[3]],"_locus_stat_",fl_ext,sep=""))
    #Fst_WC_loci
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[4]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[4]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[4]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[4]],"_locus_stat_",fl_ext,sep=""))
    }
    if(exists("jjj", where=".GlobalEnv")==TRUE){
      rm(jjj, pos=".GlobalEnv")
    }
    if(exists("accDat", where=".GlobalEnv")==TRUE){
      rm(accDat, pos=".GlobalEnv")
    }
    if(exists("bs_res", where=".GlobalEnv")==TRUE){
      rm(bs_res, pos=".GlobalEnv")
    }
    if(exists("lso123", where=".GlobalEnv")==TRUE){
      rm(lso123, pos=".GlobalEnv")
    }
    if(exists("sp.header", where=".GlobalEnv")==TRUE){
      rm(sp.header, pos=".GlobalEnv")
    }
    if(exists("plot_data321", where=".GlobalEnv")==TRUE){
      rm(plot_data321, pos=".GlobalEnv")
    }
    #rm(jjj,accDat,bs_res,lso123,sp.header,pos=".GlobalEnv")
    
  } else if(jjj$plot_loci==FALSE && jjj$plot_pw==TRUE){
    accDat<<-jjj$accDat
    pw_res<<-jjj$pw_res
    pwso<<-jjj$pwso
    pw<<-jjj$pw
    plot_data321<<-jjj$plot_data321
    if(length(pwso) > 150){
      image.size <- "2400x1200"
    } else {
      image.size=img
    }
    #Gst_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[1]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[1]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[1]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[1]],"_pairwise_stats_",fl_ext,sep=""))
    #G'st_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[2]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[2]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[2]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[2]],"_pairwise_stats_",fl_ext,sep=""))
    #Djost_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[3]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[3]],
                               image.size=image.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[3]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",fl_ext,sep=""))
    #Fst_WC_pw
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[4]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[4]],
                                 image.size=image.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[4]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[4]],"_pairwise_stats_",
                   fl_ext,sep=""))
    }
    
    if(exists("jjj", where=".GlobalEnv")==TRUE){
      rm(jjj, pos=".GlobalEnv")
    }
    if(exists("accDat", where=".GlobalEnv")==TRUE){
      rm(accDat, pos=".GlobalEnv")
    }
    if(exists("pw_res", where=".GlobalEnv")==TRUE){
      rm(pw_res, pos=".GlobalEnv")
    }
    if(exists("pwso", where=".GlobalEnv")==TRUE){
      rm(pwso, pos=".GlobalEnv")
    }
    if(exists("sp.header", where=".GlobalEnv")==TRUE){
      rm(sp.header, pos=".GlobalEnv")
    }
    if(exists("plot_data321", where=".GlobalEnv")==TRUE){
      rm(plot_data321, pos=".GlobalEnv")
    }
    if(exists("pw", where=".GlobalEnv")==TRUE){
      rm(pw, pos=".GlobalEnv")
    }
    #rm(jjj,accDat,plot_data,pw,pw_res,pwso,sp.header,pos=".GlobalEnv")
    
  } else if(jjj$plot_loci==TRUE && jjj$plot_pw==TRUE){
    bs_res<<-jjj$bs_res
    lso123<<-jjj$lso123
    accDat<<-jjj$accDat
    pw_res<<-jjj$pw_res
    pwso<<-jjj$pwso
    pw<<-jjj$pw
    plot_data321<<-jjj$plot_data321
    if(length(lso123) > 150 && length(pwso) > 150){
      pwimage.size <- "2400x1200"
      locimage.size <- "2400x1200"
    } else if(length(lso123) > 150 && length(pwso) <= 150){
      pwimage.size <- img
      locimage.size <- "2400x1200"
    } else if(length(lso123) <= 150 && length(pwso) > 150){
      pwimage.size <- "2400x1200"
      locimage.size <- img
    } else {
      locimage.size <- img
      pwimage.size <- img
    }
    #Gst_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[1]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[1]],
                               image.size=locimage.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[1]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[1]],"_locus_stat_",fl_ext,sep=""))
    #G'st_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[2]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[2]],
                               image.size=locimage.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[2]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[2]],"_locus_stat_",fl_ext,sep=""))
    #Djost_loci
    suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[3]],
                               x.pos=jjj$x.pos_loci,
                               y.pos=jjj$y.pos_loci[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_loci[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_loci[[3]],
                               image.size=locimage.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_loci[[3]],
                                                "_locus_stat_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_loci[[3]],"_locus_stat_",fl_ext,sep=""))
    #Fst_WC_loci
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_loci[[4]],
                                 x.pos=jjj$x.pos_loci,
                                 y.pos=jjj$y.pos_loci[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_loci[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_loci[[4]],
                                 image.size=locimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_loci[[4]],
                                                  "_locus_stat_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_loci[[4]],"_locus_stat_",fl_ext,sep=""))
    }
    #Gst_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[1]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[1]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[1]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[1]],
                               image.size=pwimage.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[1]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[1]],"_pairwise_stats_",
                 fl_ext,sep=""))
    #G'st_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[2]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[2]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[2]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[2]],
                               image.size=pwimage.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[2]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[2]],"_pairwise_stats_",fl_ext,sep=""))
    #Djost_pw
    suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[3]],
                               x.pos=jjj$x.pos_pw,
                               y.pos=jjj$y.pos_pw[[3]],
                               xy.type="points",
                               plot.extras=jjj$plot.extras_pw[[3]],
                               mai.mat=NA,
                               mai.prc=FALSE,
                               xy.labels=jjj$xy.labels_pw[[3]],
                               image.size=pwimage.size,
                               spot.radius=5,
                               fname.root=paste(jjj$fn_pre_pw[[3]],
                                                "_pairwise_stats_",sep=""),
                               dir=jjj$direct,
                               window.size="2100x1000"))
    unlink(paste(jjj$direct,jjj$fn_pre_pw[[3]],"_pairwise_stats_",
                 fl_ext,sep=""))
    #Fst_WC_pw
    if(jjj$fst==TRUE){
      suppressWarnings(imagesend(plot.call=jjj$plot.call_pw[[4]],
                                 x.pos=jjj$x.pos_pw,
                                 y.pos=jjj$y.pos_pw[[4]],
                                 xy.type="points",
                                 plot.extras=jjj$plot.extras_pw[[4]],
                                 mai.mat=NA,
                                 mai.prc=FALSE,
                                 xy.labels=jjj$xy.labels_pw[[4]],
                                 image.size=pwimage.size,
                                 spot.radius=5,
                                 fname.root=paste(jjj$fn_pre_pw[[4]],
                                                  "_pairwise_stats_",sep=""),
                                 dir=jjj$direct,
                                 window.size="2100x1000"))
      unlink(paste(jjj$direct,jjj$fn_pre_pw[[4]],"_pairwise_stats_",
                   fl_ext,sep=""))
    }
    if(exists("jjj", where=".GlobalEnv")==TRUE){
      rm(jjj, pos=".GlobalEnv")
    }
    if(exists("accDat", where=".GlobalEnv")==TRUE){
      rm(accDat, pos=".GlobalEnv")
    }
    if(exists("pw_res", where=".GlobalEnv")==TRUE){
      rm(pw_res, pos=".GlobalEnv")
    }
    if(exists("pwso", where=".GlobalEnv")==TRUE){
      rm(pwso, pos=".GlobalEnv")
    }
    if(exists("sp.header", where=".GlobalEnv")==TRUE){
      rm(sp.header, pos=".GlobalEnv")
    }
    if(exists("plot_data321", where=".GlobalEnv")==TRUE){
      rm(plot_data321, pos=".GlobalEnv")
    }
    if(exists("pw", where=".GlobalEnv")==TRUE){
      rm(pw, pos=".GlobalEnv")
    }
    if(exists("bs_res", where=".GlobalEnv")==TRUE){
      rm(bs_res, pos=".GlobalEnv")
    }
    if(exists("lso123", where=".GlobalEnv")==TRUE){
      rm(lso123, pos=".GlobalEnv")
    }
  }
}
################################################################################
# plotter end                                                                  #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# in.calc, a wrapper function for the calculation of locus informativeness     #
################################################################################
#' @export
in.calc<-function(infile, outfile = NULL, gp = 3, bs_locus = FALSE,
                  bs_pairwise = FALSE, bootstraps = 0, plot = FALSE,
                  parallel = FALSE){
  .Deprecated(new = "inCalc", msg = "This function name is no longer in use. Please use 'inCalc' instead. \nSee ?inCalc for usage details.", 
              old = "in.calc")
}
################################################################################
# in.calc end                                                                  #
################################################################################
#
#
#
#
#
################################################################################
# readGenepop.user, a usable function for basic population parameters          #
################################################################################
#' @export
readGenepop.user<- function (infile = NULL, gp = 3, bootstrap = FALSE) {
  .Deprecated(new = "readGenepop", msg = "This function name is no longer in use. Please use 'readGenepop' instead. \nSee ?readGenepop for usage details.", 
              old = "readGenepop.user")
}
################################################################################
# End readGenepop.user: deprecated
################################################################################
#
#
#
#
#
################################################################################
# readGenepop, a usable function for basic population parameters               #
################################################################################
#' @export
readGenepop <- function (infile=NULL, gp=3, bootstrap=FALSE) {
  gp=gp
  infile=infile
  bootstrap=bootstrap
  data1 <- fileReader(infile)
  if(gp == 3){
    data1[data1==0]<-NA
    data1[data1=="999999"]<-NA
    data1[data1=="000000"]<-NA
    data1[data1=="NANA"]<-NA
  } else if(gp == 2){
    data1[data1==0]<-NA
    data1[data1=="9999"]<-NA
    data1[data1=="0000"]<-NA
    data1[data1=="NA"]<-NA
  }
  raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
  pop_weights<- 1/pop_sizes
  
  n_harmonic<-npops/sum(pop_weights)
  
  N<-pop_sizes
  
  nloci<- (pop_pos[1]-2)
  if(nloci != (ncol(raw_data)-1)){
    stop("Check your input file for formatting errors!")
  }
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list<-list()
  for (i in 1:npops){
    pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                   2:(nloci+1)])
  }
  # check if all populations have at least some data at loci
  extCheck <- sapply(1:length(pop_list), function(i){
    sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
  })
  # remove pops with no data (for Erin Landguth 12/12)
  if (sum(extCheck) > 0){
    npops <- npops - sum(extCheck)
    pop_list <- pop_list[-(which(extCheck == TRUE))]
    pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
    pop_names <- pop_names[-(which(extCheck == TRUE))]
    pop_weights <- pop_weights[-(which(extCheck == TRUE))]
    N <- N[-(which(extCheck == TRUE))]
    #raw_data fix
    noPop <- which(extCheck == TRUE)
    indexer <- lapply(noPop, function(i){
      (pop_pos[i] + 1):(pop_pos[(i+1)])
    })
    indexer <- unlist(indexer)
    raw_data <- raw_data[-(indexer), ]    
  }
  
  
  if (gp==3) {
    plMake<-function(x){
      out <- matrix(sprintf("%06g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "0000NA"] <- "    NA"
      }
      return(out)
    }
  } else if (gp==2) {
    plMake<-function(x){
      out <- matrix(sprintf("%04g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "00NA"] <- "  NA"
      }
      return(out)
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  
  
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  
  
  if(bootstrap == T){
    bs<-function(x){
      return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
    }
    pop_list<-lapply(pop_list, bs)
  }  
  
  ###vectorize loci_pop_sizes#####################################################
  
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }#
  pre_loci_pop_sizes<-lapply(pop_list,lps)#
  pls<-matrix(ncol=nloci,nrow=npops)#
  for(i in 1:length(pre_loci_pop_sizes)){#
    pls[i,]<-pre_loci_pop_sizes[[i]]#
  }#
  #convert pls to loci_pop_sizes format
  loci_pop_sizes<-split(pls,col(pls))
  
  
  #vectorized loci_pop_weights##################################################
  
  pre_loc_weights<- 1/pls
  loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
  loci_harm_N<-npops/colSums(pre_loc_weights)
  
  #end vectorized loci_pop_weights##############################################
  
  ###vectorize pop_alleles########################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles####################################################
  
  #vectorize allele_names#######################################################
  
  alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res<-list()
    for(i in 1:ncol(x[[1]])){
      res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    }
    return(res)
  }
  
  allele_names<-lapply(pop_alleles,alln)
  
  
  loci_combi<-allele_names[[1]]
  for(j in 1:nloci){
    for(i in 2:npops){
      loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
    }
  }
  
  #all_alleles vectorized#######################################################
  
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  
  #end all_alleles vectorized###################################################
  
  aa<-all_alleles
  aa<-lapply(aa, FUN=`list`, npops)
  afMatrix<-function(x){
    np<-x[[2]]
    z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
    rownames(z)<-x[[1]]
    return(z)
  }
  allele_freq<-lapply(aa,afMatrix)
  
  
  #combine pop_alleles
  parbind<-function(x){
    rbind(x[[1]],x[[2]])
  }
  pa1<-lapply(pop_alleles, parbind)
  #create a function to tabulate the occurance of each allele
  afTab<-function(x){
    lapply(1:ncol(x), function(i){
      return(table(x[,i]))
    })
  }
  actab<-lapply(pa1, afTab)
  
  afs<-function(x){
    afsint<-function(y){
      length(na.omit(y))/2
    }
    apply(x,2,afsint)
  }
  indtyppop<-lapply(pa1,afs)
  #calculate allele frequencies
  afCalcpop<-lapply(1:length(actab), function(x){
    lapply(1:length(actab[[x]]),function(y){
      actab[[x]][[y]]/(indtyppop[[x]][y]*2)
    })
  })
  #assign allele freqs to frequency matrices
  obs_count<-allele_freq
  for(i in 1:npops){
    for(j in 1:nloci){
      allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
      obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
    }
  }
  
  
  
  indtyp<-list()
  for(i in 1:nloci){
    indtyp[[i]]<-vector()
  }
  for(i in 1:npops){
    for(j in 1:nloci){
      indtyp[[j]][i]<-indtyppop[[i]][j]
    }
  }
  
  if(bootstrap==T){
    ind_vectors<-list()
    for(i in 1:npops){
      ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
    }
    
    
    pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                     ncol=(nloci+1))
    pre_data[1,]<-c("Title",rep("\t",nloci))
    for(i in 2:(nloci+1)){
      pre_data[i,1]<-loci_names[(i-1)]
    }
    pop_data<-list()
    for(i in 1:npops){
      pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                  cbind(ind_vectors[[i]],pop_list[[i]])),
                            ncol=(nloci+1))
    }
    bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
    for(i in 2:npops){
      bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
    }
    bs_data_file<-data.frame(bs_data_file)
  }
  nalleles<-vector()
  for(i in 1:nloci){
    nalleles[i]<- nrow(allele_freq[[i]])
  }
  ##############################################################################
  if(bootstrap==T){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         #locs=locs,
         bs_file=bs_data_file,
         obs_allele_num=obs_count)
  } else if(bootstrap==F){
    list(npops=npops, 
         nloci=nloci, 
         pop_alleles=pop_alleles, 
         pop_list=pop_list,
         loci_names=loci_names, 
         pop_pos=pop_pos, 
         pop_sizes=pop_sizes,
         allele_names=allele_names,
         all_alleles=all_alleles,
         allele_freq=allele_freq,
         raw_data=raw_data,
         loci_harm_N=loci_harm_N,
         n_harmonic=n_harmonic,
         pop_names=pop_names,
         indtyp=indtyp,
         nalleles=nalleles,
         #locs=locs,
         obs_allele_num=obs_count)
  }
}
################################################################################
# readGenepop end                                                              #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# pre.divLowMemory, a low memory consumption function for locus bootstrapping  #
################################################################################
pre.divLowMemory <- function(y){
  locs <- y$locs
  fst <- y$fst
  min <- y$min
  if(is.null(min)){
    min = TRUE
  }
  # define all functions first
  # define readGenepopX function
  #############################################################################
  # readGenepopX, a function for the generation of basic population parameters#
  #############################################################################
  readGenepopX <- function (x) {
    infile=x$infile
    #gp=x$gp
    bootstrap=x$bootstrap
    # define file reader
    ###########################################################################
    # Master file reader
    ###########################################################################
    fileReader <- function (infile) {
      if (typeof(infile) == "list") {
        return(infile)
      } else if (typeof(infile) == "character") {
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if (ext == "arp") {
          convRes <- arp2gen(infile)
          if (!is.null(convRes)) {
            cat("Arlequin file converted to genepop format! \n")
            infile <- paste(flForm[1], ".gen", sep = "")
          } else {
            infile <- paste(flForm[1], ".gen", sep = "")
          }
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
          locs <- strsplit(dat[2], split = "\\s+")[[1]]
          if(length(locs != 1)){
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
        
        
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if (popLoc[1] == 3) {
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:length(dat)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x) {
          x <- unlist(strsplit(x, split = "\\s+"))
          if (is.element("", x)) {
            x <- x[-(which(x == ""))]
          }
          if (is.element(",", x)) {
            x <- x[-(which(x == ","))]
          }
          if (length(x) != 1 && length(x) != no_col) {
            x <- paste(x, collapse = "")
          }
          if (length(x) < no_col) {
            tabs <- paste(rep(NA, (no_col - length(x))), 
                          sep = "\t", collapse = "\t")
            line <- paste(x, tabs, sep = "\t")
            line <- unlist(strsplit(line, split = "\t"))
            return(line)
          } else {
            return(x)
          }
        })
      }
      out <- as.data.frame(t(dat1))
      rownames(out) <- NULL
      return(out)
    }
    data1 <- fileReader(infile)
    if(is.null(x$gp)){
      rownames(data1) <- NULL
      data1 <- as.matrix(data1)
      # determine genepop format
      p1 <- which(toupper(data1[,1]) == "POP")[1] + 1
      gp <- as.numeric(names(sort(-table(sapply(data1[p1,(-1)], nchar)/2)))[1])
      data1 <- as.data.frame(data1)
    } else {
      gp <- x$gp
    }
    if(gp == 3){
      data1[data1==0]<-NA
      data1[data1=="999999"]<-NA
      data1[data1=="000000"]<-NA
      data1[data1=="NANA"]<-NA
    } else if(gp == 2){
      data1[data1==0]<-NA
      data1[data1=="9999"]<-NA
      data1[data1=="0000"]<-NA
      data1[data1=="NA"]<-NA
    }    
    raw_data<-data1
    npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                    which(data1[,1]=="pop")))
    pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
                which(data1[,1]=="pop"),(nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
    pop_weights<- 1/pop_sizes
    
    n_harmonic<-npops/sum(pop_weights)
    
    N<-pop_sizes
    
    nloci<- (pop_pos[1]-2)
    loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
    pop_list<-list()
    for (i in 1:npops){
      pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                     2:(nloci+1)])
    }
    # check if all populations have at least some data at loci
    extCheck <- sapply(1:length(pop_list), function(i){
      sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
    })
    if (sum(extCheck) > 0){
      npops <- npops - sum(extCheck)
      pop_list <- pop_list[-(which(extCheck == TRUE))]
      pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
      pop_names <- pop_names[-(which(extCheck == TRUE))]
      pop_weights <- pop_weights[-(which(extCheck == TRUE))]
      N <- N[-(which(extCheck == TRUE))]
      #raw_data fix
      noPop <- which(extCheck == TRUE)
      indexer <- lapply(noPop, function(i){
        (pop_pos[i] + 1):(pop_pos[(i+1)])
      })
      indexer <- unlist(indexer)
      raw_data <- raw_data[-(indexer), ]    
    }  
    if (gp==3) {
      plMake<-function(x){
        out <- matrix(sprintf("%06g",as.numeric(x)),
                      nrow = nrow(x), ncol = ncol(x))
        if (Sys.info()["sysname"] == "Darwin"){
          out[out == "0000NA"] <- "    NA"
        }
        return(out)
      }
    } else if (gp==2) {
      plMake<-function(x){
        out <- matrix(sprintf("%04g",as.numeric(x)),
                      nrow = nrow(x), ncol = ncol(x))
        if (Sys.info()["sysname"] == "Darwin"){
          out[out == "00NA"] <- "  NA"
        }
        return(out)
      }
    }
    suppressWarnings(pop_list<-lapply(pop_list, plMake))
    
    if (gp == 3){
      for(i in 1:npops){
        pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
      }
    } else if (gp == 2){
      for(i in 1:npops){
        pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
      }
    }
    
    if(bootstrap == T){
      bs<-function(x){
        return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
      }
      pop_list<-lapply(pop_list, bs)
    }  
    
    ###vectorize loci_pop_sizes###############################################
    
    lps<-function(x){#
      lsp_count<-as.vector(colSums(!is.na(x)))#
      return(lsp_count)#
    }#
    pre_loci_pop_sizes<-lapply(pop_list,lps)#
    pls<-matrix(ncol=nloci,nrow=npops)#
    for(i in 1:length(pre_loci_pop_sizes)){#
      pls[i,]<-pre_loci_pop_sizes[[i]]#
    }#
    #convert pls to loci_pop_sizes format
    loci_pop_sizes<-split(pls,col(pls))
    
    
    #vectorized loci_pop_weights##############################################
    
    pre_loc_weights<- 1/pls
    loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
    loci_harm_N<-npops/colSums(pre_loc_weights)
    
    #end vectorized loci_pop_weights##########################################
    
    ###vectorize pop_alleles##################################################
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    #end vectorize pop_alleles################################################
    
    #vectorize allele_names###################################################
    
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    
    allele_names<-lapply(pop_alleles,alln)
    
    
    loci_combi<-allele_names[[1]]
    for(j in 1:nloci){
      for(i in 2:npops){
        loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
      }
    }
    
    #all_alleles vectorized###################################################
    
    aaList<-function(x){
      return(sort(unique(x,decreasing=FALSE)))
    }
    all_alleles<-lapply(loci_combi,aaList)
    
    #end all_alleles vectorized###############################################
    
    aa<-all_alleles
    aa<-lapply(aa, FUN=`list`, npops)
    afMatrix<-function(x){
      np<-x[[2]]
      z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
      rownames(z)<-x[[1]]
      return(z)
    }
    allele_freq<-lapply(aa,afMatrix)
    
    
    #combine pop_alleles
    parbind<-function(x){
      rbind(x[[1]],x[[2]])
    }
    pa1<-lapply(pop_alleles, parbind)
    #create a function to tabulate the occurance of each allele
    afTab<-function(x){
      lapply(1:ncol(x), function(i){
        return(table(x[,i]))
      })
    }
    actab<-lapply(pa1, afTab)
    
    afs<-function(x){
      afsint<-function(y){
        length(na.omit(y))/2
      }
      apply(x,2,afsint)
    }
    indtyppop<-lapply(pa1,afs)
    #calculate allele frequencies
    afCalcpop<-lapply(1:length(actab), function(x){
      lapply(1:length(actab[[x]]),function(y){
        actab[[x]][[y]]/(indtyppop[[x]][y]*2)
      })
    })
    #assign allele freqs to frequency matrices
    obs_count<-allele_freq
    for(i in 1:npops){
      for(j in 1:nloci){
        allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
        obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
      }
    }
    
    indtyp<-list()
    for(i in 1:nloci){
      indtyp[[i]]<-vector()
    }
    for(i in 1:npops){
      for(j in 1:nloci){
        indtyp[[j]][i]<-indtyppop[[i]][j]
      }
    }
    
    if(bootstrap==T){
      ind_vectors<-list()
      for(i in 1:npops){
        ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
      }
      pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                       ncol=(nloci+1))
      pre_data[1,]<-c("Title",rep("\t",nloci))
      for(i in 2:(nloci+1)){
        pre_data[i,1]<-loci_names[(i-1)]
      }
      pop_data<-list()
      for(i in 1:npops){
        pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                    cbind(ind_vectors[[i]],pop_list[[i]])),
                              ncol=(nloci+1))
      }
      bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
      for(i in 2:npops){
        bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
      }
      bs_data_file<-data.frame(bs_data_file)
    }
    nalleles<-vector()
    for(i in 1:nloci){
      nalleles[i]<- nrow(allele_freq[[i]])
    }
    ##########################################################################
    list(pop_list = pop_list,
         npops = npops,
         nloci = nloci,
         pop_sizes = pop_sizes,
         pop_alleles = pop_alleles,
         all_alleles = all_alleles,
         allele_freq = allele_freq,
         loci_harm_N = loci_harm_N,
         loci_names = loci_names,
         pop_names = pop_names,
         indtyp = indtyp,
         gp = gp)
  }
  ############################################################################
  # readGenepopX end                                                          #
  ############################################################################
  #
  #
  ############################################################################
  data1 <- readGenepopX(y)
  ############################################################################
  if(fst){
    # define the fst function
    ##########################################################################
    # fstWC: a function co calculate weir and cockerhams fis, fit, and fst
    ##########################################################################
    fstWC<-function(x){
      badData <- sapply(x$indtyp, function(y){
        is.element(0, y)
      })
      if(sum(badData) > 0){
        nl <- x$nloci - (sum(badData))
      } else{
        nl <- x$nloci
      }
      gdData<-which(!badData)
      badData<-which(badData)
      if (nl == 1) {
        all_genot<-x$pop_list[[1]][,gdData]
        if(x$npops > 1){
          for(i in 2:x$npops){
            all_genot <- c(all_genot, x$pop_list[[i]][,gdData])
          }
        }
        all_genot <- matrix(all_genot, ncol = 1)
      } else {
        all_genot<-matrix(x$pop_list[[1]][,gdData], ncol = length(gdData))
        if(x$npops > 1){
          for(i in 2:x$npops){
            all_genot<-rbind(all_genot, x$pop_list[[i]][,gdData])
          }
        }
      }
      genot<-apply(all_genot,2,unique)
      genot<-lapply(genot, function(x){
        if (sum(is.na(x))>0){
          y<-which(is.na(x)==TRUE)
          x_new<-x[-y]
          return(x_new)
        } else {
          return(x)
        }
      })
      #count genotypes
      
      genoCount<-list()
      for(i in 1:ncol(all_genot)){
        genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
        for(j in 1:length(genot[[i]])){
          genoCount[[i]][,j]<-length(which(all_genot[,i] == genot[[i]][j]))
        }
        if (x$gp==3){
          colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
                                          substr(genot[[i]],4,6),sep="")
        } else if (x$gp==2){
          colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
                                          substr(genot[[i]],3,4),sep="")
        }
      }
      
      h_sum<-list()
      for(i in 1:ncol(all_genot)){
        h_sum[[i]]<-vector()
        cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
        for(j in 1:length(x$all_alleles[[gdData[i]]])){
          het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
          het_id2<-lapply(het_id1, sum)
          het_id2<-as.vector(het_id2)
          het_id3<-which(het_id2==1)
          h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
        }
      }
      indtyp_tot<-lapply(x$indtyp, sum)
      kk_hsum <- lapply(1:ncol(all_genot), function(i){
        list(h_sum[[i]], indtyp_tot[[gdData[i]]])
      })
      kk_hbar<-lapply(kk_hsum, function(x){
        return(x[[1]]/x[[2]])
      })
      
      pdat <- lapply(1:ncol(all_genot), function(i){
        list(x$allele_freq[[gdData[i]]], x$indtyp[[gdData[i]]])
      })
      
      kk_p<-lapply(pdat, function(x){
        if(is.null(x[[1]])==FALSE){
          apply(x[[1]], 1, function(y){
            y*(2*x[[2]])
          })
        }
      })
      res<-matrix(0,(x$nloci+1),3)
      colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
      rownames(res)<-c(x$loci_names, "All")
      A<-vector()
      a<-vector()
      b<-vector()
      c<-vector()
      for(i in 1:ncol(all_genot)){
        kknbar<-indtyp_tot[[gdData[i]]]/x$npops
        kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
                 indtyp_tot[[gdData[i]]])/(x$npops-1)
        kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
        kkptild[kkptild=="NaN"]<-NA
        kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
        kks2<-colSums(x$indtyp[[gdData[i]]] * (kkptild-rep(kkpbar, each = x$npops))^2)/((x$npops-1)*kknbar)
        kkA <- kkpbar * (1-kkpbar)-(x$npops-1)*kks2/x$npops
        kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
        kkb <- (kknbar/(kknbar - 1))*(kkA-((2*kknbar-1)/(4*kknbar))*kk_hbar[[i]])
        #kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
        kkc<-kk_hbar[[i]]/2
        A[i]<-sum(kkA)
        a[i]<-sum(kka)
        b[i]<-sum(kkb)
        c[i]<-sum(kkc)
        res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
        res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
        res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
      }
      res[res=="NaN"]<-NA
      res[res==0.000]<-NA
      sumA<-sum(na.omit(A))
      suma<-sum(na.omit(a))
      sumb<-sum(na.omit(b))
      sumc<-sum(na.omit(c))
      res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
      res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
      res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
      #res[is.na(res)]<-NaN
      list(Fstats=res,
           multiLoc<-res[(x$nloci+1),])
    }
    ##########################################################################
    # end fstWC
    ##########################################################################
    
    if(locs == TRUE){
      fstats <- fstWC(data1)[[1]]
    }else if (locs == FALSE){
      fstats <- fstWC(data1)[[2]]
    }
  }
  ##############################################################################
  # create 'easy use' objects from data1 (readGenepopX output)
  # pl = pop_list
  pl<-data1$pop_list
  # np = npops
  np<-data1$npops
  # nl = nloci
  nl<-data1$nloci
  # ps = pop sizes
  ps<-data1$pop_sizes
  # pa = pop alleles
  pa<-data1$pop_alleles
  # ant = allele names total
  ant<-data1$all_alleles
  # af = allele frequencies
  af<-data1$allele_freq
  # lnharm = locus harmonic sample size
  lnharm<-round(as.numeric(data1$loci_harm_N), 4)
  # ln = locus names
  ln<-data1$loci_names
  # pn = population names
  pn<-data1$pop_names
  # ntpl = number (of individuals) typed per locus
  nt<-data1$indtyp
  
  # remove data1 to save ram
  rm(data1)
  # garbage collect
  zz <- gc(reset = TRUE)
  rm(zz)
  ##############################################################################
  #observed heterozygosity count vectorize######################################
  
  ohcFUN<-function(x){
    lapply(1:ncol(x[[1]]), function(y){
      (x[[1]][,y]!=x[[2]][,y])*1 #multiply by 1 to conver logical to numeric
    })
  }
  ohc_data<-lapply(pa, ohcFUN)
  ohcConvert<-function(x){
    matrix(unlist(x),nrow=length(x[[1]]))
  }
  ohc<-lapply(ohc_data,ohcConvert)
  
  
  #end observed heterozygosity count vectorize##################################
  #exhmf & exhtf vectorize######################################################
  
  # calculate Heterozygosity exp
  tsapply <- function(...){t(sapply(...))}
  exhtf <- tsapply(af, function(x){
    apply(x, 2, function(y){
      1 - (sum(y^2))
    })
  })
  
  #end exhmf & exhtf vectorize##################################################
  #mean frequency vectorize#####################################################
  
  mf<-lapply(af,function(x){
    rowSums(x)/np  
  })
  ht<-sapply(mf, function(x){
    1- sum(x^2)
  })
  ht[ht=="NaN"]<-NA
  
  #end mean frequency vectorize#################################################
  
  
  ###end locus stats legacy code
  #locus stats vectorize########################################################
  
  hs<-round(rowSums(exhtf)/np,4)
  hs_est<-round(hs*((2*lnharm)/((2*lnharm)-1)),4)
  ht_est<-round((ht + (hs_est/(2*lnharm*np))),4)
  ht_est[ht_est=="NaN"]<-NA
  hst<-(ht-hs)/(1-hs)
  dst<-ht-hs
  gst<-dst/ht
  djost<-((ht-hs)/(1-hs))*(np/(np-1))
  hst_est<-(ht_est-hs_est)/(1-hs_est)
  dst_est<-ht_est- hs_est
  gst_est<-(ht_est-hs_est)/ht_est
  gst_max<-((np-1)*(1-hs))/(np-1+hs)
  gst_est_max<-(((np-1)*(1-hs_est))/(np-1+hs_est))
  gst_hedrick<-gst/gst_max
  gst_est_hedrick<-gst_est/gst_est_max
  gst_est_hedrick[gst_est_hedrick > 1] <- 1
  djost_est<-(np/(np-1))*((ht_est-hs_est)/(1 - hs_est))
  
  #end locus stats vectorize####################################################
  # Across all loci stats #
  ht_mean<-round(mean(na.omit(ht)),4)
  hs_mean<-round(mean(hs),4)
  gst_all<-round((ht_mean-hs_mean)/ht_mean,4)
  gst_all_max<-round(((np-1)*(1-hs_mean))/(np-1+hs_mean),4)
  gst_all_hedrick<-round(gst_all/gst_all_max,4)
  djost_all<-round(((ht_mean-hs_mean)/(1-hs_mean))*(np/(np-1)),4)
  ##############################################################################
  # Across all loci estimated stats #
  hs_est_mean<-round(mean(hs_est),4)
  ht_est_mean<-round(mean(na.omit(ht_est)),4)
  gst_est_all<-round((ht_est_mean-hs_est_mean)/ht_est_mean,4)
  gst_est_all_max<-round((((np-1)*(1-hs_est_mean))/(np-1+hs_est_mean)),4)
  gst_est_all_hedrick<-round(gst_est_all/gst_est_all_max,4)
  gst_est_all_hedrick[gst_est_all_hedrick > 1] <- 1
  #djost_est_all<-round((np/(np-1))*((ht_est_mean-hs_est_mean)/
  #(1 - hs_est_mean)),4)
  if (nl == 1){
    djost_est_all <- round(djost_est,4)
  } else {
    djost_est_all<-round(1/((1/mean(na.omit(djost_est))+(var(na.omit(djost_est))*
                                                           ((1/mean(na.omit(djost_est)))^3)))),4)
  }
  djost_est[djost_est==0]<-NaN
  djost[djost==0]<-NaN
  ##############################################################################
  if(fst == TRUE){
    if(locs == TRUE && min == FALSE){  
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else if (locs == TRUE && min == TRUE){
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else if (locs == FALSE && min == FALSE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    } else if(locs == FALSE && min == TRUE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn,
           fstats=fstats)
    }
  } else {
    if(locs==T && min == FALSE){  
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn)
    } else if (locs == TRUE && min == TRUE){
      list(hs=hs,
           hst=hst,
           dst=dst,
           gst=gst,
           djost=djost,
           hs_est=hs_est,
           ht_est=ht_est,
           hst_est=hst_est,
           dst_est=dst_est,
           gst_est=gst_est,
           djost_est=djost_est,
           gst_max=gst_max,
           gst_est_max=gst_est_max,
           gst_hedrick=gst_hedrick,
           gst_est_hedrick=gst_est_hedrick,
           ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn)
    } else if (locs == FALSE && min == FALSE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           pop_list=pl,
           pop_names=pn)
    } else if(locs == FALSE && min == TRUE){
      list(ht_mean=ht_mean,
           hs_mean=hs_mean,
           gst_all=gst_all,
           gst_all_max=gst_all_max,
           gst_all_hedrick=gst_all_hedrick,
           djost_all=djost_all,
           hs_est_mean=hs_est_mean,
           ht_est_mean=ht_est_mean,
           gst_est_all=gst_est_all,
           gst_est_all_max=gst_est_all_max,
           pop_sizes=ps,
           gst_est_all_hedrick=gst_est_all_hedrick,
           djost_est_all=djost_est_all,
           locus_names=ln,
           locus_harmonic_N=lnharm,
           npops=np,
           nloci=nl,
           #pop_list=pl,
           pop_names=pn)
    }
  }
}
################################################################################
# end pre.divLowMemory                                                         #
################################################################################

#
#
#
#
#
#
#
#
#
################################################################################
# corPlot, plot the relationship between divPart stats and number of alleles   #
################################################################################
#' @export
corPlot<-function(x,y){
  x=x
  y=y
  par(mfrow=c(2,2))
  par(mar=c(4,5,2,2))
  sigStar <- function(x){
      if(x$p.value < 0.001) {
          return("***")
      } else if (x$p.value < 0.01) {
          return("**")
      } else if (x$p.value < 0.05) {
          return("*")
      } else {
          return("ns")
      }
  }
  #Fst
  plot(y[[2]][1:(nrow(y[[2]])-1),8]~x[[16]],
       pch=16,xlab=expression(N[alleles]),ylab=expression(hat(theta)),
       ylim=c(0,1),las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),8]~x[[16]]),col="red",lwd=2)
  cor1<-cor.test(y[[2]][1:(nrow(y[[2]])-1),8],x[[16]])
  sig <- sigStar(cor1)
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor1$estimate[[1]],3),
       " ", sig, sep=""),cex=2)
  #gst
  plot(y[[2]][1:(nrow(y[[2]])-1),4]~x[[16]],pch=16,
       xlab=expression(N[alleles]),ylab=expression(G[st]),ylim=c(0,1),
       las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),4]~x[[16]]),col="red",lwd=2)
  cor2<-cor.test(y[[2]][1:(nrow(y[[2]])-1),4],x[[16]])
  sig <- sigStar(cor2)
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor2$estimate[[1]],3),
       " ", sig, sep=""),cex=2)
  #g'st
  plot(y[[2]][1:(nrow(y[[2]])-1),5]~x[[16]],pch=16,
       xlab=expression(N[alleles]),ylab=expression("G'"[st]),ylim=c(0,1),
       las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),5]~x[[16]]),col="red",lwd=2)
  cor3<-cor.test(y[[2]][1:(nrow(y[[2]])-1),5],x[[16]])
  sig <- sigStar(cor3)
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor3$estimate[[1]],3),
       " ", sig, sep=""),cex=2)
  #D
  plot(y[[2]][1:(nrow(y[[2]])-1),6]~x[[16]],pch=16,
       xlab=expression(N[alleles]),ylab=expression(D[est]),ylim=c(0,1),
       las=1)
  abline(lm(y[[2]][1:(nrow(y[[2]])-1),6]~x[[16]]),col="red",lwd=2)
  cor4<-cor.test(y[[2]][1:(nrow(y[[2]])-1),6],x[[16]])
  sig <- sigStar(cor4)
  text(x=max(x[[16]])/1.5,y=0.8,
       labels=paste("r = ",round(cor4$estimate[[1]],3),
       " ", sig, sep=""),cex=2)
}
################################################################################
# end corPlot                                                                  #
################################################################################
#
#
#
#
#
#
#
#
#
################################################################################
# difPlot, plot all pairwise population pairs                                  #
################################################################################
#' @export
difPlot <- function (x, outfile= NULL, interactive = FALSE) {
  on=outfile
  inta <- interactive
  #output from divPart
  require("plotrix")
  if(!is.element("pairwise", names(x))){
    stop(paste("To plot pairwise differentiation the argument 'pairwise' in ",
               "either 'divPart' or fastDivPart' must be TRUE!", sep = "\n"))
  }
  if(!is.element("gstEst", names(x$pairwise))){
    idx <- c(4, 7, 5, 6)
    fun <- "divPart"
  } else {
    idx <- c(1, 4, 2, 3)
    fun <- "fastDivPart"
  }
  
  if(is.null(on) == TRUE && inta == TRUE){
    of = paste(getwd(),"/", sep = "")
  } else {
    suppressWarnings(dir.create(paste(getwd(), "/", 
                                      on, "-[diveRsity]", "/", sep="")))
    of=paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
  }
  
  if(!exists("inta",-1)){
    inta <- FALSE
  }
  if(inta == TRUE) {
    sp.header<-list()
    colleer<-list()
    colleer <<- colorRampPalette(c("blue","white"))
    require(sendplot)
    direct<-of
    pwc<-combn(ncol(x[[3]][[1]]),2)
    pwNames<-paste(colnames(x[[3]][[1]])[pwc[1,]],
                   colnames(x[[3]][[1]])[pwc[2,]],
                   sep=' vs ')
    
    gst_lab <- round(as.vector(x[[3]][[idx[1]]]), 4)
    gst_lab <- na.omit(gst_lab)
    collab111<-list()
    #
    if(length(x[[3]]) > 6 && fun == "divPart"){
      fst_lab <- round(as.vector(x[[3]][[idx[2]]]), 4)
      fst_lab<-na.omit(fst_lab)
    } else if(length(x[[3]]) == 4 && fun == "fastDivPart"){
      fst_lab <- round(as.vector(x[[3]][[idx[2]]]), 4)
      fst_lab<-na.omit(fst_lab)
    }
    #
    gpst_lab <- round(as.vector(x[[3]][[idx[3]]]), 4)
    gpst_lab<-na.omit(gpst_lab)
    #
    Dest_lab <- round(as.vector(x[[3]][[idx[4]]]), 4)
    Dest_lab<-na.omit(Dest_lab)
    #
    
    fl_ext<-c(".tif","Dot.png","Dot.tif")
    if (length(x[[3]]) > 6 && fun == "divPart"){
      xy.labels <-  data.frame(pops = pwNames,
                               Nei_Gst = gst_lab,
                               Weir_Theta = fst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    } else if(length(x[[3]]) == 4 && fun == "fastDivPart"){
      xy.labels <-  data.frame(pops = pwNames,
                               Nei_Gst = gst_lab,
                               Weir_Theta = fst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    } else {
      xy.labels <-  data.frame(pop = pwNames,
                               Nei_Gst = gst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    }
    #Nei Gst
    abx<-list()
    abx<<-x
    collab111 <<- c(round(min(gst_lab),3),
                    round(mean(gst_lab),3),
                    round(max(gst_lab),3))
    
    plot.call <- paste("image(1:nrow(abx[[3]][[", 
                       idx[1],
                       "]]), 1:nrow(abx[[3]][[", 
                       idx[1], 
                       "]]), abx[[3]][[",
                       idx[1], 
                       "]], ylab='', xlab='', main='Pairwise Gst', ",
                       "xaxt='n', yaxt='n', col = colleer(50), ",
                       "las=1, cex.main=3)", sep = "")
    ##
    plot.extras <- paste("color.legend(nrow(abx[[3]][[",
                         idx[1], 
                        "]])/5, nrow(abx[[3]][[",
                         idx[1],
                         "]])/3, nrow(abx[[3]][[",
                         idx[1],
                         "]])/4, nrow(abx[[3]][[",
                         idx[1],
                         "]])/1.2, collab111, rect.col=colleer(50), ",
                         "gradient='y', cex=3)", sep = "")
    ##
    suppressWarnings(imagesend(plot.call=plot.call,
                               x.pos=pwc[2,],
                               y.pos=pwc[1,],
                               xy.type="points",
                               xy.labels = xy.labels,
                               plot.extras=plot.extras,
                               fname.root="Gst_matrix",
                               dir=of,
                               image.size="1050x800",
                               font.size=18,
                               spot.radius = 10,
                               font.type = "Arial",
                               window.size="1100x800"))
    #clean up
    unlink(paste(of,"Gst_matrix",fl_ext,sep=""))
    #
    #
    #
    #
    #
    #Fst
    if(length(x[[3]]) > 6 && fun == "divPart"){
      collab111 <<- c(round(min(fst_lab),3),
                      round(mean(fst_lab),3),
                      round(max(fst_lab),3))
      plot.call <- paste("image(1:nrow(abx[[3]][[",
                         idx[2],
                         "]]),1:nrow(abx[[3]][[",
                         idx[2],
                         "]]), abx[[3]][[",
                         idx[2],
                         "]],ylab = '',xlab = '',xaxt = 'n',yaxt = 'n', ",
                         "main = 'Pairwise Fst', col = colleer(50), ",
                         "las = 1,cex.main = 3)", sep = "")
      ##
      plot.extras <- paste("color.legend(nrow(abx[[3]][[",
                           idx[2],
                           "]])/5, nrow(abx[[3]][[",
                           idx[2],
                           "]])/3, nrow(abx[[3]][[",
                           idx[2],
                           "]])/4, nrow(abx[[3]][[",
                           idx[2],
                           "]])/1.2, collab111, rect.col=colleer(50), ",
                           "gradient='y', cex=3)", sep = "")
      #
      suppressWarnings(imagesend(plot.call=plot.call,
                                 x.pos=pwc[2,],
                                 y.pos=pwc[1,],
                                 xy.type="points",
                                 xy.labels = xy.labels,
                                 plot.extras=plot.extras,
                                 fname.root="Fst_matrix",
                                 dir=of,
                                 image.size="1050x800",
                                 font.size=18,
                                 spot.radius = 10,
                                 font.type ="Arial",
                                 window.size="1100x800"))
      #clean up
      unlink(paste(of,"Fst_matrix",fl_ext,sep=""))
    } else if(length(x[[3]]) == 4 && fun == "fastDivPart"){
      collab111 <<- c(round(min(fst_lab),3),
                      round(mean(fst_lab),3),
                      round(max(fst_lab),3))
      plot.call <- paste("image(1:nrow(abx[[3]][[",
                         idx[2],
                         "]]),1:nrow(abx[[3]][[",
                         idx[2],
                         "]]), abx[[3]][[",
                         idx[2],
                         "]],ylab = '',xlab = '',xaxt = 'n',yaxt = 'n', ",
                         "main = 'Pairwise Fst', col = colleer(50), ",
                         "las = 1,cex.main = 3)", sep = "")
      ##
      plot.extras <- paste("color.legend(nrow(abx[[3]][[",
                           idx[2],
                           "]])/5, nrow(abx[[3]][[",
                           idx[2],
                           "]])/3, nrow(abx[[3]][[",
                           idx[2],
                           "]])/4, nrow(abx[[3]][[",
                           idx[2],
                           "]])/1.2, collab111, rect.col=colleer(50), ",
                           "gradient='y', cex=3)", sep = "")
      #
      suppressWarnings(imagesend(plot.call=plot.call,
                                 x.pos=pwc[2,],
                                 y.pos=pwc[1,],
                                 xy.type="points",
                                 xy.labels = xy.labels,
                                 plot.extras=plot.extras,
                                 fname.root="Fst_matrix",
                                 dir=of,
                                 image.size="1050x800",
                                 font.size=18,
                                 spot.radius = 10,
                                 font.type ="Arial",
                                 window.size="1100x800"))
      #clean up
      unlink(paste(of,"Fst_matrix",fl_ext,sep=""))
    }
    #
    #
    #
    #
    #
    #G'st
    collab111 <<- c(round(min(gpst_lab),3),
                    round(mean(gpst_lab),3),
                    round(max(gpst_lab),3))
    plot.call <- paste("image(1:nrow(abx[[3]][[",
                       idx[3],
                       "]]),1:nrow(abx[[3]][[",
                       idx[3],
                       "]]), abx[[3]][[",
                       idx[3],
                       "]], ylab='', xlab='', xaxt='n', yaxt='n', ",
                       "main='Pairwise Gst (Hedrick)',col = colleer(50), ",
                       "las=1,cex.main=3)", sep = "")
    
    plot.extras <- paste("color.legend(nrow(abx[[3]][[",
                         idx[3],
                         "]])/5,nrow(abx[[3]][[",
                         idx[3],
                         "]])/3, nrow(abx[[3]][[",
                         idx[3],
                         "]])/4,nrow(abx[[3]][[",
                         idx[3],
                         "]])/1.2,collab111, rect.col=colleer(50), ",
                         "gradient='y',cex=3)", sep = "")
    ##
    suppressWarnings(imagesend(plot.call=plot.call,
                               x.pos=pwc[2,],
                               y.pos=pwc[1,],
                               xy.type="points",
                               xy.labels = xy.labels,
                               plot.extras=plot.extras,
                               fname.root="G_prime_st_matrix",
                               dir=of,
                               image.size="1050x800",
                               font.size=18,
                               spot.radius = 10,
                               font.type = "Arial",
                               window.size="1100x800"))
    #clean up
    unlink(paste(of,"G_prime_st_matrix",fl_ext,sep=""))
    #
    #
    #
    #
    #
    #
    #Dest
    collab111 <<- c(round(min(Dest_lab),3),
                    round(mean(Dest_lab),3),
                    round(max(Dest_lab),3))
    plot.call <- paste("image(1:nrow(abx[[3]][[",
                       idx[4],
                       "]]),1:nrow(abx[[3]][[",
                       idx[4],
                       "]]), abx[[3]][[", 
                       idx[4],
                       "]], ylab = '',xlab = '',xaxt = 'n', ",
                       "yaxt = 'n',main = 'Pairwise D (Jost)', ",
                       "col = colleer(50),las=1,cex.main=3)", sep = "")
    
    plot.extras <- paste("color.legend(nrow(abx[[3]][[",
                         idx[4],
                         "]])/5,nrow(abx[[3]][[",
                         idx[4],
                         "]])/3, nrow(abx[[3]][[",
                         idx[4],
                         "]])/4,nrow(abx[[3]][[",
                         idx[4],
                         "]])/1.2,collab111, rect.col = colleer(50), ",
                         "gradient='y',cex=3)", sep = "")
    ##
    suppressWarnings(imagesend(plot.call=plot.call,
                               x.pos=pwc[2,],
                               y.pos=pwc[1,],
                               xy.type="points",
                               xy.labels = xy.labels,
                               plot.extras=plot.extras,
                               fname.root="D_matrix_",
                               dir=of,
                               image.size="1050x800",
                               font.size=18,
                               spot.radius = 10,
                               font.type = "Arial",
                               window.size="1100x800"))
    #lean up
    
    unlink(paste(of,"D_matrix_",fl_ext,sep=""))
    
  } else {
    
    #
    if(length(x[[3]]) > 6 && fun == "divPart"){
      par(mfrow=c(2,2))
    } else if(length(x[[3]]) == 4 && fun == "fastDivPart"){
      par(mfrow=c(2,2))
    } else {
      par(mfrow=c(3,1))
    }
    colleer<-colorRampPalette(c("blue","white"))
    cols<-colleer(50)
    #Gst
    image(1:nrow(x[[3]][[idx[1]]]),
          1:nrow(x[[3]][[idx[1]]]),
          x[[3]][[idx[1]]],
          ylab="population",
          xlab="population",
          main="Pairwise Gst",
          col=cols,
          las=1)
    gst<-as.vector(x[[3]][[idx[1]]])
    gst<-as.vector(na.omit(gst))
    collab111<-c(round(min(gst),3),
                 round(mean(gst),3),
                 round(max(gst),3))
    
    color.legend(nrow(x[[3]][[idx[1]]])/5,
                 nrow(x[[3]][[idx[1]]])/3,
                 nrow(x[[3]][[idx[1]]])/4,
                 nrow(x[[3]][[idx[1]]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
    if(length(x[[3]]) > 6 && fun == "divPart"){
      #Fst
      image(1:nrow(x[[3]][[idx[2]]]),
            1:nrow(x[[3]][[idx[2]]]),
            x[[3]][[idx[2]]],
            ylab="population",
            xlab="population",
            main="Pairwise Theta",
            col = cols,
            las=1)
      fst<-as.vector(x[[3]][[idx[2]]])
      fst<-as.vector(na.omit(fst))
      collab111<-c(round(min(fst),3),round(mean(fst),3),round(max(fst),3))
      
      color.legend(nrow(x[[3]][[idx[2]]])/5,
                   nrow(x[[3]][[idx[2]]])/3,
                   nrow(x[[3]][[idx[2]]])/4,
                   nrow(x[[3]][[idx[2]]])/1.2,
                   collab111,
                   cols,
                   gradient="y")
    } else if(length(x[[3]]) == 4 && fun == "fastDivPart"){
      #Fst
      image(1:nrow(x[[3]][[idx[2]]]),
            1:nrow(x[[3]][[idx[2]]]),
            x[[3]][[idx[2]]],
            ylab="population",
            xlab="population",
            main="Pairwise Theta",
            col = cols,
            las=1)
      fst<-as.vector(x[[3]][[idx[2]]])
      fst<-as.vector(na.omit(fst))
      collab111<-c(round(min(fst),3),round(mean(fst),3),round(max(fst),3))
      
      color.legend(nrow(x[[3]][[idx[2]]])/5,
                   nrow(x[[3]][[idx[2]]])/3,
                   nrow(x[[3]][[idx[2]]])/4,
                   nrow(x[[3]][[idx[2]]])/1.2,
                   collab111,
                   cols,
                   gradient="y")
    }
    #Hedrick's Gst
    image(1:nrow(x[[3]][[idx[3]]]),
          1:nrow(x[[3]][[idx[3]]]),
          x[[3]][[idx[3]]],
          ylab="population",
          xlab="population",
          main="Pairwise G'st",
          col = cols)
    gprimest<-as.vector(x[[3]][[idx[3]]])
    gprimest<-as.vector(na.omit(gprimest))
    collab111<-c(round(min(gprimest),3),
                 round(mean(gprimest),3),
                 round(max(gprimest),3))
    
    color.legend(nrow(x[[3]][[idx[3]]])/5,
                 nrow(x[[3]][[idx[3]]])/3,
                 nrow(x[[3]][[idx[3]]])/4,
                 nrow(x[[3]][[idx[3]]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
    #Jost's D
    image(1:nrow(x[[3]][[idx[4]]]),
          1:nrow(x[[3]][[idx[4]]]),
          x[[3]][[idx[4]]],
          ylab="population",
          xlab="population",
          main="Pairwise Jost's D",
          col = cols,
          las=1)
    D<-as.vector(x[[3]][[idx[4]]])
    D<-as.vector(na.omit(D))
    collab111<-c(round(min(D),3),
                 round(mean(D),3),
                 round(max(D),3))
    
    color.legend(nrow(x[[3]][[idx[4]]])/5,
                 nrow(x[[3]][[idx[4]]])/3,
                 nrow(x[[3]][[idx[4]]])/4,
                 nrow(x[[3]][[idx[4]]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
  }
  if(exists("abx", where=".GlobalEnv")==TRUE){
    rm(abx, pos=".GlobalEnv")
  }
  if(exists("collab111", where=".GlobalEnv")==TRUE){
    rm(collab111, pos=".GlobalEnv")
  }
  if(exists("colleer", where=".GlobalEnv")==TRUE){
    rm(colleer, pos=".GlobalEnv")
  }
  if(exists("sp.header", where=".GlobalEnv")==TRUE){
    rm(sp.header, pos=".GlobalEnv")
  }
  par(mfrow = c(1,1))
}
################################################################################
# end dif.Plot                                                                 #
################################################################################
#
#
#
#
#
#
#
###############################################################################
#                 chiCalc, a function for the assessment of                   #
#             population heterogeniety from microsatellite data.              #
#        Input data should be given in the 2 or 3 digit genepop format        #
#                       By Kevin Keenan, QUB, 2013                            #
###############################################################################
#' @export
chiCalc <- function(infile = NULL, outfile = NULL, gp = 3, minFreq = NULL){
  inputs <- list(infile = infile, gp = gp, bootstrap = FALSE)
  minFreq <- minFreq
  
  # define read genepop function
  #############################################################################
  # readGenepopX, a function for the generation of basic population parameters #
  #############################################################################
  readGenepopX <- function (x) {
    #gp=x$gp
    infile=x$infile
    bootstrap=x$bootstrap
    data1 <- fileReader(infile)
    if(is.null(x$gp)){
      rownames(data1) <- NULL
      data1 <- as.matrix(data1)
      # determine genepop format
      p1 <- which(toupper(data1[,1]) == "POP")[1] + 1
      gp <- as.numeric(names(sort(-table(sapply(data1[p1,(-1)], nchar)/2)))[1])
      data1 <- as.data.frame(data1)
    } else {
      gp <- x$gp
    }
    if(gp == 3){
      data1[data1==0]<-NA
      data1[data1=="999999"]<-NA
      data1[data1=="000000"]<-NA
      data1[data1=="NANA"]<-NA
    } else if(gp == 2){
      data1[data1==0]<-NA
      data1[data1=="9999"]<-NA
      data1[data1=="0000"]<-NA
      data1[data1=="NA"]<-NA
    }    
    raw_data<-data1
    npops<-length(which(toupper(data1[,1]) == "POP"))
    pop_pos<- c(which(toupper(data1[,1]) == "POP"),(nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
    pop_weights<- 1/pop_sizes
    
    n_harmonic<-npops/sum(pop_weights)
    
    N<-pop_sizes
    
    nloci<- (pop_pos[1]-2)
    loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
    pop_list<-list()
    for (i in 1:npops){
      pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                     2:(nloci+1)])
    }
    # check if all populations have at least some data at loci
    extCheck <- sapply(1:length(pop_list), function(i){
      sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
    })
    if (sum(extCheck) > 0){
      npops <- npops - sum(extCheck)
      pop_list <- pop_list[-(which(extCheck == TRUE))]
      pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
      pop_names <- pop_names[-(which(extCheck == TRUE))]
      pop_weights <- pop_weights[-(which(extCheck == TRUE))]
      N <- N[-(which(extCheck == TRUE))]
      #raw_data fix
      noPop <- which(extCheck == TRUE)
      indexer <- lapply(noPop, function(i){
        (pop_pos[i] + 1):(pop_pos[(i+1)])
      })
      indexer <- unlist(indexer)
      raw_data <- raw_data[-(indexer), ]    
    }  
    if (gp==3) {
      plMake<-function(x){
        out <- matrix(sprintf("%06g",as.numeric(x)),
                      nrow = nrow(x), ncol = ncol(x))
        if (Sys.info()["sysname"] == "Darwin"){
          out[out == "0000NA"] <- "    NA"
        }
        return(out)
      }
    } else if (gp==2) {
      plMake<-function(x){
        out <- matrix(sprintf("%04g",as.numeric(x)),
                      nrow = nrow(x), ncol = ncol(x))
        if (Sys.info()["sysname"] == "Darwin"){
          out[out == "00NA"] <- "  NA"
        }
        return(out)
      }
    }
    suppressWarnings(pop_list<-lapply(pop_list, plMake))
    
    if (gp == 3){
      for(i in 1:npops){
        pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
      }
    } else if (gp == 2){
      for(i in 1:npops){
        pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
      }
    }
    
    if(bootstrap == T){
      bs<-function(x){
        return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
      }
      pop_list<-lapply(pop_list, bs)
    }  
    
    ###vectorize loci_pop_sizes###############################################
    
    lps<-function(x){#
      lsp_count<-as.vector(colSums(!is.na(x)))#
      return(lsp_count)#
    }#
    pre_loci_pop_sizes<-lapply(pop_list,lps)#
    pls<-matrix(ncol=nloci,nrow=npops)#
    for(i in 1:length(pre_loci_pop_sizes)){#
      pls[i,]<-pre_loci_pop_sizes[[i]]#
    }#
    #convert pls to loci_pop_sizes format
    loci_pop_sizes<-split(pls,col(pls))
    
    
    #vectorized loci_pop_weights##############################################
    
    pre_loc_weights<- 1/pls
    loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
    loci_harm_N<-npops/colSums(pre_loc_weights)
    
    #end vectorized loci_pop_weights##########################################
    
    ###vectorize pop_alleles##################################################
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    #end vectorize pop_alleles################################################
    
    #vectorize allele_names###################################################
    
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    
    allele_names<-lapply(pop_alleles,alln)
    
    
    loci_combi<-allele_names[[1]]
    for(j in 1:nloci){
      for(i in 2:npops){
        loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
      }
    }
    
    #all_alleles vectorized###################################################
    
    aaList<-function(x){
      return(sort(unique(x,decreasing=FALSE)))
    }
    all_alleles<-lapply(loci_combi,aaList)
    
    #end all_alleles vectorized###############################################
    
    aa<-all_alleles
    aa<-lapply(aa, FUN=`list`, npops)
    afMatrix<-function(x){
      np<-x[[2]]
      z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
      rownames(z)<-x[[1]]
      return(z)
    }
    allele_freq<-lapply(aa,afMatrix)
    
    
    #combine pop_alleles
    parbind<-function(x){
      rbind(x[[1]],x[[2]])
    }
    pa1<-lapply(pop_alleles, parbind)
    #create a function to tabulate the occurance of each allele
    afTab<-function(x){
      lapply(1:ncol(x), function(i){
        return(table(x[,i]))
      })
    }
    actab<-lapply(pa1, afTab)
    
    afs<-function(x){
      afsint<-function(y){
        length(na.omit(y))/2
      }
      apply(x,2,afsint)
    }
    indtyppop<-lapply(pa1,afs)
    #calculate allele frequencies
    afCalcpop<-lapply(1:length(actab), function(x){
      lapply(1:length(actab[[x]]),function(y){
        actab[[x]][[y]]/(indtyppop[[x]][y]*2)
      })
    })
    #assign allele freqs to frequency matrices
    obs_count<-allele_freq
    for(i in 1:npops){
      for(j in 1:nloci){
        allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
        obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
      }
    }
    
    
    
    indtyp<-list()
    for(i in 1:nloci){
      indtyp[[i]]<-vector()
    }
    for(i in 1:npops){
      for(j in 1:nloci){
        indtyp[[j]][i]<-indtyppop[[i]][j]
      }
    }
    
    if(bootstrap==T){
      ind_vectors<-list()
      for(i in 1:npops){
        ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
      }
      
      
      pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                       ncol=(nloci+1))
      pre_data[1,]<-c("Title",rep("\t",nloci))
      for(i in 2:(nloci+1)){
        pre_data[i,1]<-loci_names[(i-1)]
      }
      pop_data<-list()
      for(i in 1:npops){
        pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                    cbind(ind_vectors[[i]],pop_list[[i]])),
                              ncol=(nloci+1))
      }
      bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
      for(i in 2:npops){
        bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
      }
      bs_data_file<-data.frame(bs_data_file)
    }
    nalleles<-vector()
    for(i in 1:nloci){
      nalleles[i]<- nrow(allele_freq[[i]])
    }
    ##########################################################################
    if(bootstrap==T){
      list(npops=npops, 
           nloci=nloci, 
           pop_alleles=pop_alleles, 
           pop_list=pop_list,
           loci_names=loci_names, 
           pop_pos=pop_pos, 
           pop_sizes=pop_sizes,
           allele_names=allele_names,
           all_alleles=all_alleles,
           allele_freq=allele_freq,
           raw_data=raw_data,
           loci_harm_N=loci_harm_N,
           n_harmonic=n_harmonic,
           pop_names=pop_names,
           indtyp=indtyp,
           nalleles=nalleles,
           #locs=locs,
           bs_file=bs_data_file,
           obs_allele_num=obs_count)
    } else if(bootstrap==F){
      list(npops=npops, 
           nloci=nloci, 
           pop_alleles=pop_alleles, 
           pop_list=pop_list,
           loci_names=loci_names, 
           pop_pos=pop_pos, 
           pop_sizes=pop_sizes,
           allele_names=allele_names,
           all_alleles=all_alleles,
           allele_freq=allele_freq,
           raw_data=raw_data,
           loci_harm_N=loci_harm_N,
           n_harmonic=n_harmonic,
           pop_names=pop_names,
           indtyp=indtyp,
           nalleles=nalleles,
           #locs=locs,
           obs_allele_num=obs_count)
    }
  }
  ############################################################################
  # readGenepopX end                                                          #
  ############################################################################
  #
  ############################################################################
  # Main body of chiCalc function                                            #
  ############################################################################
  # extract observed allele numbers using readGenepopX
  dat <- readGenepopX(inputs)
  allNum <- dat$obs_allele_num
  # Calculate column sums 
  csum <- lapply(allNum, function(x){
    apply(x, 2, function(y){
      return(sum(y))
    })
  })
  # Calculate row sums
  rsum <- lapply(allNum, function(x){
    apply(x, 1, function(y){
      return(sum(y))
    })
  })
  # Calculate expected numbers
  expNum <- lapply(1:dat$nloci, function(i){
    mat <- sapply(1:length(rsum[[i]]), function(j){
      cols <- sapply(csum[[i]], function(x){
        (x * rsum[[i]][j])/sum(rsum[[i]])
      })
    })
    return(t(mat))
  })
  # Calculate chi values per allele per population
  chisq <- lapply(1:dat$nloci, function(i){
    out <- sapply(1:ncol(allNum[[i]]), function(j){
      return(((allNum[[i]][,j] - expNum[[i]][,j])^2)/expNum[[i]][,j])
    })
    out <- matrix(out, ncol = dat$npops)
    rownames(out) <- dat$all_alleles[[i]]
    return(out)
  })
  # Calculate chi values across populations
  alleleChi <- lapply(chisq, function(x){
    apply(x, 1, FUN = 'sum')
  })
  # Assign loci names to alleleChi
  names(alleleChi) <- dat$loci_names
  # Link mean allele frequency and Chisq values
  chiFreq <- lapply(1:dat$nloci, function(i){
    out <- matrix(rbind(alleleChi[[i]], rowMeans(dat$allele_freq[[i]])), 
                  nrow = 2)
    dimnames(out) <- list(c("Chi", "Freq"), dat$all_alleles[[i]])
    return(round(out, 4))
  })
  ############################################################################
  # Calculate chi values for all allele data
  ############################################################################
  # Calculate locus sums (chi)
  locsumAll <- round(sapply(alleleChi, FUN = 'sum'),4)
  # Calculate locus degrees of freedom (Basic formula = k-1)
  locDf <- sapply(dat$all_alleles, function(x){
    return(length(x) - 1)
  })
  # Calulcate p values using allele alleles
  pAll <- round(pchisq(q = locsumAll, df = locDf, lower.tail = FALSE),4)
  # Create visual significance indicator 
  sigStar <- function(x){
    if(is.na(x)){
      return(NA)
    } else if(x < 0.001) {
      return("***")
    } else if (x < 0.01) {
      return("**")
    } else if (x < 0.05) {
      return("*")
    } else {
      return("ns")
    }
  }
  sig <- sapply(pAll, sigStar)
  # Compile locus results into dataframe
  resOut <- matrix(cbind(dat$loci_names, locsumAll, locDf, pAll, sig), ncol = 5)
  
  
  # Calculate overall statistics
  
  chiTotal <- round(sum(locsumAll), 4)
  dfTotal <- round(sum(locDf), 4)
  pTotal <- round(pchisq(q = chiTotal, df = dfTotal, lower.tail = FALSE), 4)
  sigTot <- sigStar(pTotal)
  if(pTotal == 0){
    pTotal <- '0.0000'
  }
  # Add overall chisq stats to results dataframe
  resOut <- rbind(resOut, c("Overall", chiTotal, dfTotal, pTotal, sigTot))
  colnames(resOut) <- c("locus", "chisq", "df", "p.value", "signif")
  
  ############################################################################
  # Calculate stats for all alleles above nominal level(s)
  ############################################################################
  # check if minFreq is a single value
  if(length(minFreq) == 1){
    chifreqMin <- sapply(chiFreq, function(x){
      x[, which(x["Freq",] >= minFreq)]
    })
    locsumMin <- sapply(chifreqMin, function(x){
      sum(x[1,])
    })
    locdfMin <- sapply(chifreqMin, function(x){
      ncol(x) - 1
    })
    pMin <- round(pchisq(q = locsumMin, df = locdfMin, lower.tail = FALSE), 4)
    sigMin <- sapply(pMin, sigStar)
    chitotMin <- round(sum(locsumAll), 4)
    dftotMin <- round(sum(locDf), 4)
    ptotMin <- round(pchisq(q = chiTotal, df = dfTotal, lower.tail = FALSE), 4)
    sigtotMin <- sigStar(ptotMin)
    if(ptotMin == 0){
      ptotMin <- '0.0000'
    }
    resMin <- matrix(cbind(locsumMin, locdfMin, pMin, sigMin), ncol = 4)
    colnames(resMin) <- c(paste("chisq(", minFreq, ")", sep = ""),
                          paste("df(", minFreq, ")", sep = ""),
                          paste("p.value(", minFreq, ")", sep = ""),
                          paste("signif(", minFreq, ")", sep = ""))
    resMin <- rbind(resMin, c(chitotMin, dftotMin, ptotMin, sigtotMin))
    resOut <- cbind(resOut, resMin)
  } else if (length(minFreq) > 1){
    res <- lapply(minFreq, function(z){
      chifreqMin <- sapply(chiFreq, function(x){
        x[, which(x["Freq",] >= z)]
      })
      locsumMin <- sapply(chifreqMin, function(x){
        if(ncol(x) == 0 || is.null(ncol(x))){
          return(NA)
        } else {
          return(sum(x[1,]))
        }
      })
      locdfMin <- sapply(chifreqMin, function(x){
        if(ncol(x) == 0 || is.null(ncol(x))){
          return(NA)
        } else {
          ncol(x) - 1
        }
      })
      pMin <- round(pchisq(q = locsumMin, df = locdfMin, 
                           lower.tail = FALSE), 4)
      sigMin <- sapply(pMin, sigStar)
      chitotMin <- round(sum(locsumMin, na.rm = TRUE), 4)
      dftotMin <- round(sum(locdfMin, na.rm = TRUE), 4)
      ptotMin <- round(pchisq(q = chitotMin, df = dftotMin, 
                              lower.tail = FALSE), 4)
      sigtotMin <- sigStar(ptotMin)
      if(ptotMin == 0){
        ptotMin <- '0.0000'
      }
      resMin <- matrix(cbind(locsumMin, locdfMin, pMin, sigMin), ncol = 4)
      colnames(resMin) <- c(paste("chisq(", z, ")", sep = ""),
                            paste("df(", z, ")", sep = ""),
                            paste("p.value(", z, ")", sep = ""),
                            paste("signif(", z, ")", sep = ""))
      resMin <- rbind(resMin, c(chitotMin, dftotMin, ptotMin, sigtotMin))
      return(resMin)
    })
    for(i in 1:length(minFreq)){
      resOut <- cbind(resOut, res[[i]])
    }
  }
  on <- outfile
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[diveRsity]","/",sep="")))
    of <- paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
    write.table(resOut, append = FALSE,
                 file = paste(of, "[chi].txt", sep = ""), 
                 sep = "\t", eol = "\n", quote = FALSE, col.names = TRUE,
                row.names = FALSE)
  }
  return(resOut)
}

###############################################################################
# END
###############################################################################
#
#
#
#
###############################################################################
# try to include diveRsity online
#' @export
divOnline <- function(){
    runApp(system.file('diveRsity-online', package = 'diveRsity'))
}
################################################################################
# END
################################################################################
#
#
#
#
#
#
#
# try to include microPlexer app
#' @export
microPlexer <- function(){
  runApp(system.file('microPlexer', package = 'diveRsity'))
}
################################################################################
# END
################################################################################
#
#
#
#
#
#
#
################################################################################
# Calculate basic stats
################################################################################
#' @export
divBasic <- function (infile = NULL, outfile = NULL, gp = 3, 
                      bootstraps = NULL) {
  
  on = outfile
  # create a results dir
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[diveRsity]","/",sep="")))
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
  }
  
  data1 <- fileReader(infile)
  data1[data1 == 0] <- NA
  data1[data1 == "999999"] <- NA
  data1[data1 == "000000"] <- NA
  data1[data1 == "9999"] <- NA
  data1[data1 == "0000"] <- NA
  #raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"),(nrow(data1)+1))
  pop_sizes <- sapply(1:npops, function(i){
    pop_pos[(i+1)] - pop_pos[i]-1
  })
  minSize <- min(pop_sizes) 
  pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list <- lapply(1:npops, function(i){
    return(as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                           2:(nloci+1)]))
  })
  if (gp==3) {
    plMake<-function(x){
      out <- matrix(sprintf("%06g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "0000NA"] <- "    NA"
      }
      return(out)
    }
  } else if (gp==2) {
    plMake<-function(x){
      out <- matrix(sprintf("%04g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "00NA"] <- "  NA"
      }
      return(out)
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  
  
  # define a function for calculating allelic richness 
  ARfun <- function(pop_list){
    
    bser<-function(x){
      return(matrix(x[sample(nrow(x), minSize, replace = TRUE), ],ncol=ncol(x)))
    }
    
    pop_list<-lapply(pop_list, bser)
    
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    allele_names<-lapply(pop_alleles,alln)
    Alls <- lapply(allele_names, function(x){
      sapply(x, function(y){
        length(y)
      })
    })
    nAlls <- matrix(unlist(Alls), ncol = npops)
    
    return(nAlls)
  }
  ############################# END AR function ###############################
  # Calculate allelic richness
  ARdata <- replicate(1000, ARfun(pop_list))
  
  AR <- apply(ARdata, 2, function(x){
    round(rowMeans(x), 2)
  })
  meanAR <- apply(ARdata, 3, colMeans, na.rm = TRUE)
  arLCI <- apply(meanAR, 1, quantile, probs = 0.025, na.rm = TRUE)
  arUCI <- apply(meanAR, 1, quantile, probs = 0.975, na.rm = TRUE)
  #locSD <- apply(ARdata, c(1,2), sd, na.rm = TRUE)
  locSD <- rbind(arLCI, arUCI)
  colnames(locSD) <- pop_names
  rownames(locSD) <- c("Lower_CI", "Upper_CI")
  ###vectorize loci_pop_sizes#################################################
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }
  locPopSize <- sapply(pop_list,lps)
  ###vectorize pop_alleles####################################################
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
      return(pl)
    }
  }
  pop_alleles<-lapply(pop_list,pl_ss)
  #end vectorize pop_alleles##################################################
  #vectorize allele_names#####################################################
  pop_alleles<-lapply(pop_list,pl_ss)
  # calcluate the observed heterozygosity
  ohcFUN<-function(x){
    lapply(1:ncol(x[[1]]), function(y){
      (x[[1]][,y]!=x[[2]][,y])*1 #multiply by 1 to conver logical to numeric
    })
  }
  ohc_data<-lapply(pop_alleles, ohcFUN)
  ohcConvert<-function(x){
    matrix(unlist(x),nrow=length(x[[1]]))
  }
  ohc<-lapply(ohc_data,ohcConvert)
  rm(ohc_data)
  hetObs <- sapply(ohc, function(x){
    apply(x, 2, function(y){
      sum(na.omit(y))/length(na.omit(y))
    })
  })
  # End
  alln <- function(x){ # where x is the object pop_alleles (returned by pl_ss())
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    })
  }
  allele_names<-sapply(pop_alleles,alln)
  # Count the number of alleles observed in each population sample per locus
  obsAlls <- apply(allele_names, 2, function(x){
    sapply(x, function(y){
      length(y)
    })
  })
  # Calculate expected He
  if(npops == 1){
    loci_combi <- allele_names[,1]
  } else {
    loci_combi <- apply(allele_names, 1, FUN = 'unlist')
  }
  # fix loci_combi for SNP format
  if(is.matrix(loci_combi)){
    loci_combi <- lapply(1:ncol(loci_combi), function(i){
      return(loci_combi[,i])
    })
  }
  aaList<-function(x){
    return(sort(unique(x,decreasing=FALSE)))
  }
  all_alleles<-lapply(loci_combi,aaList)
  # Create allele frequency holders
  allele_freq <- lapply(1:ncol(pop_list[[1]]), function(i){
    Nrow <- length(all_alleles[[i]])
    Ncol <- length(pop_list)
    mat <- matrix(rep(0,(Ncol * Nrow)), ncol = Ncol)
    rownames(mat) <- all_alleles[[i]]
    return(mat)
  })
  # rbind pop_alleles
  pa1 <- lapply(pop_alleles, function(x){
    rbind(x[[1]],x[[2]])
  })
  
  # Count alleles
  actab <- lapply(pa1, function(x){
    lapply(1:ncol(x), function(i){
      table(x[,i])
    })
  })
  # Count the number of individuals typed per locus per pop
  indtyppop <- lapply(pa1, function(x){
    apply(x, 2, function(y){
      length(na.omit(y))/2
    })
  })
  #calculate allele frequencies
  afCalcpop<-sapply(1:length(actab), function(x){
    sapply(1:length(actab[[x]]),function(y){
      list(actab[[x]][[y]]/(indtyppop[[x]][y]*2))
    })
  })
  preFreq <- lapply(1:nrow(afCalcpop), function(i){
    lapply(1:ncol(afCalcpop), function(j){
      afCalcpop[i,j][[1]]
    })
  })
  rm(afCalcpop)  # remove afCalcpop
  # Assign allele frequencies per locus
  for(i in 1:nloci){
    for(j in 1:npops){
      allele_freq[[i]][names(preFreq[[i]][[j]]), j] <- preFreq[[i]][[j]]
    }
  }
  # calculate Heterozygosity exp
  if(npops > 1){
    tsapply <- function(...){t(sapply(...))}
    hetExp <- tsapply(allele_freq, function(x){
      apply(x, 2, function(y){
        1 - (sum(y^2))
      })
    })
  } else {
    hetExp <- as.matrix(sapply(allele_freq, function(x){
      apply(x, 2, function(y){
        1 - (sum(y^2))
      })
    }), ncol = 1)
  }
  
  
  totAlls <- sapply(allele_freq, FUN = "nrow")
  # Calculate the proportion of alleles per sample
  propAlls <- apply(obsAlls, 2, function(x){
    round((x/totAlls)*100, 2)
  })
  # R function to calculate expected and observed genetype 
  # numbers for HWE testing
  
  # generate all possible genotypes for each locus per population
  posGeno <- apply(allele_names, 2, function(x) {
    lapply(x, function(y) {
      if (length(y) == 0) {
        return(NA)
      } else {
        genos <- expand.grid(y, y)
        genos.sort <- t(apply(genos, 1, sort))
        genos <- unique(genos.sort)
        geno <- paste(genos[, 1], genos[, 2], sep = "")
        return(geno)
      }
    })
  })
  
  # Count the number of each genotype observed
  # define a genotype counting function
  obsGeno <- lapply(1:npops, function(i){
    lapply(1:nloci, function(j){
      sapply(posGeno[[i]][[j]], function(x){
        if(is.na(x)){
          return(NA)
        } else {
          length(which(pop_list[[i]][,j] == x))
        }
      })
    })
  })
  
  
  expGeno <- lapply(1:npops, function(i){
    lapply(1:nloci, function(j){
      sapply(posGeno[[i]][[j]], function(x){
        if(is.na(x)){
          return(NA)
        } else {
          if(gp == 3){
            allele1 <- substr(x, 1, 3)
            allele2 <- substr(x, 4, 6)
          } else {
            allele1 <- substr(x, 1, 2)
            allele2 <- substr(x, 3, 4)
          }
          Freq1 <- allele_freq[[j]][which(rownames(allele_freq[[j]]) == allele1), i]
          Freq2 <- allele_freq[[j]][which(rownames(allele_freq[[j]]) == allele2), i]
          if(allele1 != allele2){
            expFreq <- 2 * (Freq1 * Freq2)
            return(as.vector(expNum <- expFreq * locPopSize[j, i]))
          } else {
            expFreq <- Freq1^2
            return(expNum <- as.vector(expFreq * locPopSize[j, i]))
          }
        }
      })
    })
  })
  
  # Calculate chi-sq
  chiDif <- sapply(1:npops, function(i){
    sapply(1:nloci, function(j){
      if(length(obsGeno[[i]][[j]]) == 1){
        return(NA)
      } else {
        top <- (obsGeno[[i]][[j]] - expGeno[[i]][[j]])^2
        chi <- top/expGeno[[i]][[j]]
        return(round(sum(chi), 2))
      }      
    })
  })
  
  # Calculate degrees of freedom
  df <- apply(allele_names, 2,   function(x){
    sapply(x, function(y){
      k <- length(y)
      if(k == 1){
        return(NA)
      } else {
        return((k*(k-1))/2)
      }
    })
  })
  
  # Calculate HWE significance
  HWE <- sapply(1:npops, function(i){
    round(pchisq(q = chiDif[,i], df = df[,i], lower.tail = FALSE), 4)
  })
  if(!is.null(bootstraps)){
    # write a function to calculate allele freq and  obsHet from pop_alleles object
    # convert pop_alleles into a list of arrays
    pa <- lapply(pop_alleles, function(x){
      return(array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), 2)))
    })
    # fit boot function
    bootPA <- function(pa){
      if(!is.list(pa)){
        idx <- sample(dim(pa)[1], dim(pa)[1], replace = TRUE)
        pa <- pa[idx,,]
        obsHet <- apply(pa, 2, function(x){
          (x[,1] != x[, 2])*1
        })
        obsHet <- apply(obsHet, 2, function(x){
          sum(na.omit(x))/length(na.omit(x))
        })
        # calculate expected
        htExp <- apply(pa, 2, function(x){
          af <- as.vector((table(c(x[,1], x[,2]))/(length(na.omit(x[,1]))*2))^2)
          return(1 - sum(af, na.rm = TRUE))
        })
        ht <- sum(htExp, na.rm=TRUE)
        ho <- sum(obsHet, na.rm=TRUE)
        overall <- (ht-ho)/ht
        return(c((htExp-obsHet)/htExp, overall))
      } else {
        out <- lapply(pa, function(pasub){
          idx <- sample(dim(pasub)[1], dim(pasub)[1], replace = TRUE)
          pasub <- pasub[idx,,]
          obsHet <- apply(pasub, 2, function(x){
            (x[,1] != x[, 2])*1
          })
          obsHet <- apply(obsHet, 2, function(x){
            sum(na.omit(x))/length(na.omit(x))
          })
          # calculate expected
          htExp <- apply(pasub, 2, function(x){
            af <- as.vector((table(c(x[,1], x[,2]))/(length(na.omit(x[,1]))*2))^2)
            return(1 - sum(af, na.rm = TRUE))
          })
          ht <- htExp
          ho <- obsHet
          fis <- (ht-ho)/ht
          fis[is.nan(fis)] <- NA
          overall <- (sum(ht,na.rm=TRUE)-sum(ho,na.rm=TRUE))/sum(ht,na.rm=TRUE)
          fis <- c(fis, overall)
          return(fis)
        })
        return(do.call("rbind", out))
      }
    }
    # calculate base fis
    fisCalc <- function(x, y){
      ho <- colSums(x, na.rm = TRUE)
      he <- colSums(y, na.rm = TRUE)
      return((he-ho)/he)
    }
    fisLoc <- round((hetExp[-(nloci+1),] - hetObs[-(nloci+1),])/hetExp[-(nloci+1),], 4)
    fisLoc[is.nan(fisLoc)] <- NA
    fisAct <- fisCalc(hetObs, hetExp)
    # calculate fis CIs
    fisBS <- replicate(bootstraps, bootPA(pa))
    # convert fisBS into list format
    fisBSloc <- lapply(1:npops, function(i){
      return(t(fisBS[i,-(nloci+1),]))
    })
    fisBSOverall <- lapply(1:npops, function(i){
      return(as.vector((fisBS[i,nloci+1,])))
    })
    # fix the bias
    biasCor <- lapply(1:npops, function(i){
      bs <- fisBSloc[[i]]
      mnBA <- colMeans(bs, na.rm = TRUE) - fisLoc[,i]
      mnBA[is.nan(mnBA)] <- NA
      bcCor <- t(apply(bs, 1, function(x){
        return(x - mnBA)
      }))
    })
    biasCorall <- lapply(1:npops, function(i){
      bs <- fisBSOverall[[i]]
      mnBA <- mean(bs, na.rm = TRUE) - fisAct[i]
      mnBA[is.nan(mnBA)] <- NA
      return(as.vector(bs - mnBA))
    })
    
    # bias cor CIs
    bcCILoc <- lapply(biasCor, function(x){
      apply(x, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    })
    bcCIall <- lapply(biasCorall, quantile, 
                      probs = c(0.025, 0.975), na.rm = TRUE)
    nbcCILoc <- lapply(fisBSloc, function(x){
      apply(x, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    })
    nbcCIall <- lapply(fisBSOverall, quantile, 
                       probs = c(0.025, 0.975), na.rm = TRUE)
    fisLoc <- lapply(1:npops, function(i){
      return(fisLoc[,i])
    })
    # arrange data for outpup
    # define function
    tableMake <- function(fisLoc, fisAct, bcCIall, bcCILoc, nbcCIall, nbcCILoc,
                          loci_names){
      out <- data.frame(fis = c(fisLoc, fisAct),
                        lower_CI = c(nbcCILoc[1,], nbcCIall[1]),
                        upper_CI = c(nbcCILoc[2,], nbcCIall[2]),
                        BC_lower_CI = c(bcCILoc[1,], bcCIall[1]),
                        BC_upper_CI = c(bcCILoc[2,], bcCIall[2]))
      rownames(out) <- c(loci_names, "overall")
      return(round(out, 4))
    }
    
    output <- mapply(tableMake, fisLoc = fisLoc,
                     fisAct = fisAct, bcCIall = bcCIall,
                     bcCILoc = bcCILoc, nbcCIall = nbcCIall,
                     nbcCILoc = nbcCILoc,  SIMPLIFY = FALSE,
                     MoreArgs = list(loci_names = loci_names))
    
  }
  
  # Calculate over all HWE significance
  HWEall <- round(pchisq(q = colSums(chiDif), df = colSums(df), 
                         lower.tail = FALSE), 4)
  # Add pop means or totals to each stat object
  # allelic richness
  AR <- round(rbind(AR, colMeans(AR)), 2)
  dimnames(AR) <- list(c(loci_names, "overall"), pop_names)
  # Number of individuals typed per locus per pop
  locPopSize <- rbind(locPopSize, round(colMeans(locPopSize), 2))
  dimnames(locPopSize) <- list(c(loci_names, "overall"), pop_names)
  # proportion of alleles per pop
  propAlls <- round(rbind(propAlls, colMeans(propAlls)), 2)
  dimnames(propAlls) <- list(c(loci_names, "overall"), pop_names)
  # Number of alleles observed per pop
  obsAlls <- rbind(obsAlls, colSums(obsAlls))
  dimnames(obsAlls) <- list(c(loci_names, "overall"), pop_names)
  # Observed heterozygosity
  hetObs <- round(rbind(hetObs, colMeans(hetObs)), 2)
  dimnames(hetObs) <- list(c(loci_names, "overall"), pop_names)
  # Expected heterozygosity
  hetExp <- round(rbind(hetExp, colMeans(hetExp)), 2)
  dimnames(hetExp) <- list(c(loci_names, "overall"), pop_names)
  # HWE
  HWE <- rbind(HWE, HWEall)
  # Compile information into writable format
  if(!is.null(bootstraps)){
    statComp <- lapply(1:npops, function(i){
      pop <- rbind(locPopSize[,i], obsAlls[,i],
                   propAlls[,i], AR[,i], hetObs[,i],
                   hetExp[,i], HWE[,i], output[[i]][,1],
                   output[[i]][,4], output[[i]][,5])
      return(pop)
    })
    if(npops > 1){
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "Fis",
                          "Fis_Low", "Fis_High", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
      for(i in 2:npops){
        writeOut <- rbind(writeOut, 
                          cbind(c(pop_names[i], 
                                  "N", "A", "%", "Ar", "Ho", "He", "HWE", "Fis",
                                  "Fis_Low", "Fis_High", "\t"),
                                rbind(c(loci_names, "Overall"), statComp[[i]], 
                                      rep("\t", nloci+1))))
      }
    } else {
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "Fis",
                          "Fis_Low", "Fis_High", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
    }
  } else {
    statComp <- lapply(1:npops, function(i){
      pop <- rbind(locPopSize[,i], obsAlls[,i],
                   propAlls[,i], AR[,i], hetObs[,i],
                   hetExp[,i], HWE[,i])
      return(pop)
    })
    if(npops > 1){
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
      for(i in 2:npops){
        writeOut <- rbind(writeOut, 
                          cbind(c(pop_names[i], 
                                  "N", "A", "%", "Ar", "Ho", "He", "HWE", "\t"),
                                rbind(c(loci_names, "Overall"), statComp[[i]], 
                                      rep("\t", nloci+1))))
      }
    } else {
      writeOut <- cbind(c(pop_names[1], 
                          "N", "A", "%", "Ar", "Ho", "He", "HWE", "\t"), 
                        rbind(c(loci_names, "Overall"), 
                              statComp[[1]], rep("\t", nloci+1)))
    }
  }
  
  
  if (!is.null(outfile)){
    write_res<-is.element("xlsx",installed.packages()[,1])
    if(write_res){
      library("xlsx")
      write.xlsx(writeOut, file = paste(of, "[divBasic].xlsx", sep = ""),
                 sheetName = "Basic stats", col.names = FALSE,
                 row.names = FALSE, append=FALSE)
    } else {
      out<-file(paste(of, "[divBasic].txt", sep = ""), "w")
      #cat(paste(colnames(pw_bs_out),sep=""),"\n",sep="\t",file=pw_bts)
      for(i in 1:nrow(writeOut)){
        cat(writeOut[i,], "\n", file = out, sep="\t")
      }
      close(out)
    }
  }
  if(!is.null(bootstraps)){
    list(locus_pop_size = locPopSize,
         Allele_number = obsAlls,
         proportion_Alleles = propAlls,
         Allelic_richness = AR,
         Ho = hetObs,
         He = hetExp,
         HWE = HWE,
         fis = output,
         arCIs = round(locSD, 4),
         mainTab = writeOut)
  } else {
    list(locus_pop_size = locPopSize,
         Allele_number = obsAlls,
         proportion_Alleles = propAlls,
         Allelic_richness = AR,
         Ho = hetObs,
         He = hetExp,
         HWE = HWE,
         arCIs = round(locSD, 4),
         mainTab = writeOut)
  }
  
}
################################################################################
# END
################################################################################
#
#
#
#
################################################################################
# Master file reader
################################################################################
fileReader <- function(infile){
  if (typeof(infile) == "list") {
    return(infile)
  } else if (typeof(infile) == "character") {
    flForm <- strsplit(infile, split = "\\.")[[1]]
    ext <- flForm[[length(flForm)]]
    if (ext == "arp") {
      convRes <- arp2gen(infile)
      if (!is.null(convRes)) {
        cat("Arlequin file converted to genepop format! \n")
        infile <- paste(flForm[1], ".gen", sep = "")
      } else {
        infile <- paste(flForm[1], ".gen", sep = "")
      }
    }
    dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
    if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
      locs <- strsplit(dat[2], split = "\\s+")[[1]]
      if(length(locs != 1)){
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
    
    
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    no_col <- popLoc[1] - 1
    if (popLoc[1] == 3) {
      locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
      dat <- c(dat[1], locs, dat[3:length(dat)])
    }
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    no_col <- popLoc[1] - 1
    dat1 <- sapply(dat, function(x) {
      x <- unlist(strsplit(x, split = "\\s+"))
      if (is.element("", x)) {
        x <- x[-(which(x == ""))]
      }
      if (is.element(",", x)) {
        x <- x[-(which(x == ","))]
      }
      if (length(x) != 1 && length(x) != no_col) {
        x <- paste(x, collapse = "")
      }
      if (length(x) < no_col) {
        tabs <- paste(rep(NA, (no_col - length(x))), 
                      sep = "\t", collapse = "\t")
        line <- paste(x, tabs, sep = "\t")
        line <- unlist(strsplit(line, split = "\t"))
        return(line)
      } else {
        return(x)
      }
    })
  }
  out <- as.data.frame(t(dat1))
  rownames(out) <- NULL
  return(out)
}
################################################################################
# END
################################################################################
#
#
#
#
#
#
#
################################################################################
# fstWC: a function co calculate weir and cockerhams fis, fit, and fst
################################################################################
fstWC<-function(x){
  badData <- sapply(x$indtyp, function(y){
    is.element(0, y)
  })
  if(sum(badData) > 0){
    nl <- x$nloci - (sum(badData))
  } else{
    nl <- x$nloci
  }
  gdData<-which(!badData)
  badData<-which(badData)
  if (nl == 1) {
    all_genot<-x$pop_list[[1]][,gdData]
    if(x$npops > 1){
      for(i in 2:x$npops){
        all_genot <- c(all_genot, x$pop_list[[i]][,gdData])
      }
    }
    all_genot <- matrix(all_genot, ncol = 1)
  } else {
    all_genot<-matrix(x$pop_list[[1]][,gdData], ncol = length(gdData))
    if(x$npops > 1){
      for(i in 2:x$npops){
        all_genot<-rbind(all_genot, x$pop_list[[i]][,gdData])
      }
    }
  }
  genot<-apply(all_genot,2,unique)
  genot<-lapply(genot, function(x){
    if (sum(is.na(x))>0){
      y<-which(is.na(x)==TRUE)
      x_new<-x[-y]
      return(x_new)
    } else {
      return(x)
    }
  })
  #count genotypes
  
  genoCount<-list()
  for(i in 1:ncol(all_genot)){
    genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
    for(j in 1:length(genot[[i]])){
      genoCount[[i]][,j]<-length(which(all_genot[,i] == genot[[i]][j]))
    }
    if (x$gp==3){
      colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
                                      substr(genot[[i]],4,6),sep="")
    } else if (x$gp==2){
      colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
                                      substr(genot[[i]],3,4),sep="")
    }
  }
  
  h_sum<-list()
  for(i in 1:ncol(all_genot)){
    h_sum[[i]]<-vector()
    cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
    for(j in 1:length(x$all_alleles[[gdData[i]]])){
      het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
      het_id2<-lapply(het_id1, sum)
      het_id2<-as.vector(het_id2)
      het_id3<-which(het_id2==1)
      h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
    }
  }
  indtyp_tot<-lapply(x$indtyp, sum)
  kk_hsum <- lapply(1:ncol(all_genot), function(i){
    list(h_sum[[i]], indtyp_tot[[gdData[i]]])
  })
  kk_hbar<-lapply(kk_hsum, function(x){
    return(x[[1]]/x[[2]])
  })
  
  pdat <- lapply(1:ncol(all_genot), function(i){
    list(x$allele_freq[[gdData[i]]], x$indtyp[[gdData[i]]])
  })
  
  kk_p<-lapply(pdat, function(x){
    if(is.null(x[[1]])==FALSE){
      apply(x[[1]], 1, function(y){
        y*(2*x[[2]])
      })
    }
  })
  res<-matrix(0,(x$nloci+1),3)
  colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
  rownames(res)<-c(x$loci_names, "All")
  A<-vector()
  a<-vector()
  b<-vector()
  c<-vector()
  for(i in 1:ncol(all_genot)){
    kknbar<-indtyp_tot[[gdData[i]]]/x$npops
    kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
             indtyp_tot[[gdData[i]]])/(x$npops-1)
    kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
    kkptild[kkptild=="NaN"]<-NA
    kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
    kks2<-colSums(x$indtyp[[gdData[i]]]*
                    (kkptild-rep(kkpbar,each = x$npops))^2)/((x$npops-1)*kknbar)
    kkA<-kkpbar*(1-kkpbar)-(x$npops-1)*kks2/x$npops
    kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
    kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
    kkc<-kk_hbar[[i]]/2
    A[i]<-sum(kkA)
    a[i]<-sum(kka)
    b[i]<-sum(kkb)
    c[i]<-sum(kkc)
    res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
    res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
    res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
  }
  res[res=="NaN"]<-NA
  res[res==0.000]<-NA
  sumA<-sum(na.omit(A))
  suma<-sum(na.omit(a))
  sumb<-sum(na.omit(b))
  sumc<-sum(na.omit(c))
  res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
  res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
  res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
  #res[is.na(res)]<-NaN
  list(Fstats=res,
       multiLoc<-res[(x$nloci+1),])
}
################################################################################
# end fstWC
################################################################################
#
#
#
#
#
#
#
################################################################################
# fstOnly: a memory efficient function to calculate WC Fst and Fit
################################################################################
#' @export
fstOnly <- function(infile = NULL, outfile = NULL, gp = 3, 
                    bs_locus = FALSE, bs_pairwise = FALSE, 
                    bootstraps = 0, parallel = FALSE){
  # create a directory
  if(!is.null(outfile)){
    suppressWarnings(dir.create(path=paste(getwd(), "/", outfile,
                                           "-[fstWC]", "/",sep="")))
  }
  
  # define the fstWC function
  #############################################################################
  # fstWC: a function co calculate weir and cockerhams fis, fit, and fst
  #############################################################################
  fstWC<-function(x){
    badData <- sapply(x$indtyp, function(y){
      is.element(0, y)
    })
    if(sum(badData) > 0){
      nl <- x$nloci - (sum(badData))
    } else{
      nl <- x$nloci
    }
    gdData<-which(!badData)
    badData<-which(badData)
    if (nl == 1) {
      all_genot<-x$pop_list[[1]][,gdData]
      if(x$npops > 1){
        for(i in 2:x$npops){
          all_genot <- c(all_genot, x$pop_list[[i]][,gdData])
        }
      }
      all_genot <- matrix(all_genot, ncol = 1)
    } else {
      all_genot<-matrix(x$pop_list[[1]][,gdData], ncol = length(gdData))
      if(x$npops > 1){
        for(i in 2:x$npops){
          all_genot<-rbind(all_genot, x$pop_list[[i]][,gdData])
        }
      }
    }
    genot<-apply(all_genot,2,unique)
    genot<-lapply(genot, function(x){
      if (sum(is.na(x))>0){
        y<-which(is.na(x)==TRUE)
        x_new<-x[-y]
        return(x_new)
      } else {
        return(x)
      }
    })
    #count genotypes
    
    genoCount<-list()
    for(i in 1:ncol(all_genot)){
      genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
      for(j in 1:length(genot[[i]])){
        genoCount[[i]][,j]<-length(which(all_genot[,i] == genot[[i]][j]))
      }
      if (x$gp==3){
        colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
                                        substr(genot[[i]],4,6),sep="")
      } else if (x$gp==2){
        colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
                                        substr(genot[[i]],3,4),sep="")
      }
    }
    
    h_sum<-list()
    for(i in 1:ncol(all_genot)){
      h_sum[[i]]<-vector()
      cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
      for(j in 1:length(x$all_alleles[[gdData[i]]])){
        het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
        het_id2<-lapply(het_id1, sum)
        het_id2<-as.vector(het_id2)
        het_id3<-which(het_id2==1)
        h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
      }
    }
    indtyp_tot<-lapply(x$indtyp, sum)
    kk_hsum <- lapply(1:ncol(all_genot), function(i){
      list(h_sum[[i]], indtyp_tot[[gdData[i]]])
    })
    kk_hbar<-lapply(kk_hsum, function(x){
      return(x[[1]]/x[[2]])
    })
    
    pdat <- lapply(1:ncol(all_genot), function(i){
      list(x$allele_freq[[gdData[i]]], x$indtyp[[gdData[i]]])
    })
    
    kk_p<-lapply(pdat, function(x){
      if(is.null(x[[1]])==FALSE){
        apply(x[[1]], 1, function(y){
          y*(2*x[[2]])
        })
      }
    })
    res<-matrix(0,(x$nloci+1),3)
    colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
    rownames(res)<-c(x$loci_names, "All")
    A<-vector()
    a<-vector()
    b<-vector()
    c<-vector()
    for(i in 1:ncol(all_genot)){
      kknbar<-indtyp_tot[[gdData[i]]]/x$npops
      kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
               indtyp_tot[[gdData[i]]])/(x$npops-1)
      kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
      kkptild[kkptild=="NaN"]<-NA
      kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
      kks2<-colSums(x$indtyp[[gdData[i]]]*
                      (kkptild-rep(kkpbar,each = x$npops))^2)/((x$npops-1)*
                                                                 kknbar)
      kkA<-kkpbar*(1-kkpbar)-(x$npops-1)*kks2/x$npops
      kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
      kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
      kkc<-kk_hbar[[i]]/2
      A[i]<-sum(kkA)
      a[i]<-sum(kka)
      b[i]<-sum(kkb)
      c[i]<-sum(kkc)
      res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
      res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
      res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
    }
    res[res=="NaN"]<-NA
    res[res==0.000]<-NA
    sumA<-sum(na.omit(A))
    suma<-sum(na.omit(a))
    sumb<-sum(na.omit(b))
    sumc<-sum(na.omit(c))
    res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
    res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
    res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
    #res[is.na(res)]<-NaN
    list(Fstats=res,
         multiLoc<-res[(x$nloci+1),])
  }
  #############################################################################
  # end fstWC
  #############################################################################
  #
  #
  #
  # define the readGenepopX function
  readGenepopX <- function (x) {
    infile=x$infile
    #gp=x$gp
    bootstrap=x$bootstrap
    # define file reader
    ###########################################################################
    # Master file reader
    ###########################################################################
    fileReader <- function(infile){
      if (typeof(infile) == "list") {
        return(infile)
      } else if (typeof(infile) == "character") {
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if (ext == "arp") {
          convRes <- arp2gen(infile)
          if (!is.null(convRes)) {
            cat("Arlequin file converted to genepop format! \n")
            infile <- paste(flForm[1], ".gen", sep = "")
          } else {
            infile <- paste(flForm[1], ".gen", sep = "")
          }
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
          locs <- strsplit(dat[2], split = "\\s+")[[1]]
          if(length(locs != 1)){
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
        
        
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if (popLoc[1] == 3) {
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:length(dat)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x) {
          x <- unlist(strsplit(x, split = "\\s+"))
          if (is.element("", x)) {
            x <- x[-(which(x == ""))]
          }
          if (is.element(",", x)) {
            x <- x[-(which(x == ","))]
          }
          if (length(x) != 1 && length(x) != no_col) {
            x <- paste(x, collapse = "")
          }
          if (length(x) < no_col) {
            tabs <- paste(rep(NA, (no_col - length(x))), 
                          sep = "\t", collapse = "\t")
            line <- paste(x, tabs, sep = "\t")
            line <- unlist(strsplit(line, split = "\t"))
            return(line)
          } else {
            return(x)
          }
        })
      }
      out <- as.data.frame(t(dat1))
      rownames(out) <- NULL
      return(out)
    }
    data1 <- fileReader(infile)
    if(is.null(x$gp)){
      rownames(data1) <- NULL
      data1 <- as.matrix(data1)
      # determine genepop format
      p1 <- which(toupper(data1[,1]) == "POP")[1] + 1
      gp <- as.numeric(names(sort(-table(sapply(data1[p1,(-1)], nchar)/2)))[1])
      data1 <- as.data.frame(data1)
    } else {
      gp <- x$gp
    }
    if(gp == 3){
      data1[data1==0]<-NA
      data1[data1=="999999"]<-NA
      data1[data1=="000000"]<-NA
      data1[data1=="NANA"]<-NA
    } else if(gp == 2){
      data1[data1==0]<-NA
      data1[data1=="9999"]<-NA
      data1[data1=="0000"]<-NA
      data1[data1=="NA"]<-NA
    }    
    raw_data<-data1
    npops<-length(which(toupper(data1[,1]) == "POP"))
    pop_pos<- c(which(toupper(data1[,1]) == "POP"), (nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
    pop_weights<- 1/pop_sizes
    
    n_harmonic<-npops/sum(pop_weights)
    
    N<-pop_sizes
    
    nloci<- (pop_pos[1]-2)
    loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
    pop_list<-list()
    for (i in 1:npops){
      pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                                     2:(nloci+1)])
    }
    # check if all populations have at least some data at loci
    extCheck <- sapply(1:length(pop_list), function(i){
      sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
    })
    if (sum(extCheck) > 0){
      npops <- npops - sum(extCheck)
      pop_list <- pop_list[-(which(extCheck == TRUE))]
      pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
      pop_names <- pop_names[-(which(extCheck == TRUE))]
      pop_weights <- pop_weights[-(which(extCheck == TRUE))]
      N <- N[-(which(extCheck == TRUE))]
      #raw_data fix
      noPop <- which(extCheck == TRUE)
      indexer <- lapply(noPop, function(i){
        (pop_pos[i] + 1):(pop_pos[(i+1)])
      })
      indexer <- unlist(indexer)
      raw_data <- raw_data[-(indexer), ]    
    }  
    if (gp==3) {
      plMake<-function(x){
        out <- matrix(sprintf("%06g",as.numeric(x)),
                      nrow = nrow(x), ncol = ncol(x))
        if (Sys.info()["sysname"] == "Darwin"){
          out[out == "0000NA"] <- "    NA"
        }
        return(out)
      }
    } else if (gp==2) {
      plMake<-function(x){
        out <- matrix(sprintf("%04g",as.numeric(x)),
                      nrow = nrow(x), ncol = ncol(x))
        if (Sys.info()["sysname"] == "Darwin"){
          out[out == "00NA"] <- "  NA"
        }
        return(out)
      }
    }
    suppressWarnings(pop_list<-lapply(pop_list, plMake))
    
    if (gp == 3){
      for(i in 1:npops){
        pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
      }
    } else if (gp == 2){
      for(i in 1:npops){
        pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
      }
    }
    
    if(bootstrap == T){
      bs<-function(x){
        return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
      }
      pop_list<-lapply(pop_list, bs)
    }  
    
    ###vectorize loci_pop_sizes###############################################
    
    lps<-function(x){#
      lsp_count<-as.vector(colSums(!is.na(x)))#
      return(lsp_count)#
    }#
    pre_loci_pop_sizes<-lapply(pop_list,lps)#
    pls<-matrix(ncol=nloci,nrow=npops)#
    for(i in 1:length(pre_loci_pop_sizes)){#
      pls[i,]<-pre_loci_pop_sizes[[i]]#
    }#
    #convert pls to loci_pop_sizes format
    loci_pop_sizes<-split(pls,col(pls))
    
    
    #vectorized loci_pop_weights##############################################
    
    pre_loc_weights<- 1/pls
    loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
    loci_harm_N<-npops/colSums(pre_loc_weights)
    
    #end vectorized loci_pop_weights##########################################
    
    ###vectorize pop_alleles##################################################
    if (gp==3){
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
        pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
        return(pl)
      }
    } else {
      pl_ss<-function(x){  # where x is object pop_list
        pl<-list()
        pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
        pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
        return(pl)
      }
    }
    pop_alleles<-lapply(pop_list,pl_ss)
    #end vectorize pop_alleles################################################
    
    #vectorize allele_names###################################################
    
    alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
      res<-list()
      for(i in 1:ncol(x[[1]])){
        res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
      }
      return(res)
    }
    
    allele_names<-lapply(pop_alleles,alln)
    
    
    loci_combi<-allele_names[[1]]
    for(j in 1:nloci){
      for(i in 2:npops){
        loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
      }
    }
    
    #all_alleles vectorized###################################################
    
    aaList<-function(x){
      return(sort(unique(x,decreasing=FALSE)))
    }
    all_alleles<-lapply(loci_combi,aaList)
    
    #end all_alleles vectorized###############################################
    
    aa<-all_alleles
    aa<-lapply(aa, FUN=`list`, npops)
    afMatrix<-function(x){
      np<-x[[2]]
      z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
      rownames(z)<-x[[1]]
      return(z)
    }
    allele_freq<-lapply(aa,afMatrix)
    
    
    #combine pop_alleles
    parbind<-function(x){
      rbind(x[[1]],x[[2]])
    }
    pa1<-lapply(pop_alleles, parbind)
    #create a function to tabulate the occurance of each allele
    afTab<-function(x){
      lapply(1:ncol(x), function(i){
        return(table(x[,i]))
      })
    }
    actab<-lapply(pa1, afTab)
    
    afs<-function(x){
      afsint<-function(y){
        length(na.omit(y))/2
      }
      apply(x,2,afsint)
    }
    indtyppop<-lapply(pa1,afs)
    #calculate allele frequencies
    afCalcpop<-lapply(1:length(actab), function(x){
      lapply(1:length(actab[[x]]),function(y){
        actab[[x]][[y]]/(indtyppop[[x]][y]*2)
      })
    })
    #assign allele freqs to frequency matrices
    obs_count<-allele_freq
    for(i in 1:npops){
      for(j in 1:nloci){
        allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
        obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
      }
    }
    
    indtyp<-list()
    for(i in 1:nloci){
      indtyp[[i]]<-vector()
    }
    for(i in 1:npops){
      for(j in 1:nloci){
        indtyp[[j]][i]<-indtyppop[[i]][j]
      }
    }
    
    if(bootstrap==T){
      ind_vectors<-list()
      for(i in 1:npops){
        ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
      }
      pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
                       ncol=(nloci+1))
      pre_data[1,]<-c("Title",rep("\t",nloci))
      for(i in 2:(nloci+1)){
        pre_data[i,1]<-loci_names[(i-1)]
      }
      pop_data<-list()
      for(i in 1:npops){
        pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
                                    cbind(ind_vectors[[i]],pop_list[[i]])),
                              ncol=(nloci+1))
      }
      bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
      for(i in 2:npops){
        bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
      }
      bs_data_file<-data.frame(bs_data_file)
    }
    nalleles<-vector()
    for(i in 1:nloci){
      nalleles[i]<- nrow(allele_freq[[i]])
    }
    ##########################################################################
    list(pop_list = pop_list,
         npops = npops,
         nloci = nloci,
         pop_sizes = pop_sizes,
         pop_alleles = pop_alleles,
         all_alleles = all_alleles,
         allele_freq = allele_freq,
         loci_harm_N = loci_harm_N,
         loci_names = loci_names,
         pop_names = pop_names,
         indtyp = indtyp,
         gp = gp)
  }
  ############################################################################
  # readGenepopX end                                                        #
  ############################################################################
  
  #setup a parallel cluster if parallel == TRUE
  if(parallel){
    para_pack <- is.element(c("parallel","doParallel","foreach","iterators"),
                            installed.packages()[,1])
    if(sum(para_pack) != 4){
      stop("Please install all required parallel packages")
    } else {
      library("doParallel")
      cores <- detectCores()
      cl <- makeCluster(cores)
      registerDoParallel(cl)
    }
  }
  
  
  # create the baseline
  accData <- readGenepopX(list(infile = infile, gp = gp, bootstrap = FALSE))
  glbFst <- fstWC(accData)
  if(bs_locus){
    # calculate locus bootstraps
    if(bs_locus && parallel){
      input <- list(infile = infile, gp = gp, bootstrap = TRUE)
      clusterExport(cl, c("fstWC", "readGenepopX", "input"), 
                    envir = environment())
      loc_stats <- parLapply(cl, 1:bootstraps, function(...){
        gps <- readGenepopX(input)
        fst <- fstWC(gps)$Fstats
        return(fst)
      })
    } else if(bs_locus && !parallel){
      input <- list(infile = infile, gp = gp, bootstrap = TRUE)
      loc_stats <- lapply(1:bootstraps, function(...){
        gps <- readGenepopX(input)
        fst <- fstWC(gps)$Fstats
        return(fst)
      })
    }
    # compile bs_locus results
    loc_fst <- sapply(loc_stats, function(x){
      return(x[,"Fst_WC"])
    })
    loc_fit <- sapply(loc_stats, function(x){
      return(x[,"Fit_WC"])
    })
    locBS <- list(loc_fst = loc_fst,
                  loc_fit = loc_fit)
    locBSci <- lapply(locBS, function(x){
      apply(x, 1, function(y){
        return(as.vector(((sd(y)/sqrt(length(y))) * 1.96)))
      })
    })
    locBSout <- lapply(1:2, function(i){
      if(i == 1){
        cbind(actual = round(glbFst[[1]][,"Fst_WC"], 4),
              lower_CI = round(glbFst[[1]][,"Fst_WC"] - locBSci[[i]], 4),
              upper_CI = round(glbFst[[1]][,"Fst_WC"] + locBSci[[i]], 4))
        
      } else if(i == 2){
        cbind(actual = round(glbFst[[1]][,"Fit_WC"], 4),
              lower_CI = round(glbFst[[1]][,"Fit_WC"] - locBSci[[i]], 4),
              upper_CI = round(glbFst[[1]][,"Fit_WC"] + locBSci[[i]], 4))
      }
    })
    # write the locus results
    locOut <- rbind(c("actual", "lower_CI", "upper_CI"),
                    locBSout[[1]],
                    c("--", "--", "--"),
                    c("","",""),
                    c("actual", "lower_CI", "upper_CI"),
                    locBSout[[2]])
    locNames <- c("Fst", rownames(locBSout[[1]]),
                  "--", "", "Fit", rownames(locBSout[[1]]))
    locOut <- cbind(locNames, locOut)
    dimnames(locOut) <- NULL
    
    if (!is.null(outfile)){
      of = paste(getwd(), "/", outfile, "-[fstWC]", "/", sep = "")
      if(is.element("xlsx", installed.packages()[, 1])){
        # write data to excel
        # Load dependencies
        library("xlsx")
        # standard stats
        write.xlsx(locOut, file = paste(of, "[fstWC].xlsx", sep = ""),
                   sheetName = "Locus_stats", col.names = FALSE,
                   row.names = FALSE, append = FALSE)
      } else {
        # text file alternatives
        std <- file(paste(of, "Locus_stats-[fstWC].txt", sep = ""), "w")
        for(i in 1:nrow(locOut)){
          cat(locOut[i,], "\n", file = std, sep = "\t")
        }
        close(std)
      }
    }
    names(locBSout) <- c("Fst", "Fit")
  }
  
  
  ##########################################################################
  ##                              PAIRWISE                                ##
  ##########################################################################
  
  
  
  # calculate pairwise bootstraps
  if(bs_pairwise){
    pw <- combn(accData$npops,2)
    pwmat <- pw + 1
    ind_vectors <- lapply(1:accData$npops, function(x){
      rep(x, accData$pop_sizes[[x]])}
    )
    #      
    pre_data <- matrix(rep("", ((accData$nloci + 1) * (accData$nloci + 1))),
                       ncol = (accData$nloci + 1))
    pre_data[1,] <- rep("", (accData$nloci + 1))
    #
    for(i in 2:(accData$nloci + 1)){
      pre_data[i, 1] <- accData$loci_names[(i-1)]
    }
    #
    pw_data<-list()
    for (i in 1:ncol(pw)){
      pw_data[[i]]<-data.frame(rbind(pre_data,
                                     c("POP",as.vector(rep("",accData$nloci))),
                                     cbind(ind_vectors[[pw[1,i]]],
                                           matrix(noquote(accData$pop_list
                                                          [[pw[1,i]]]),
                                                  ncol=accData$nloci)),
                                     c("POP",as.vector(rep("",accData$nloci))),
                                     cbind(ind_vectors[[pw[2,i]]],
                                           matrix(noquote(accData$pop_list
                                                          [[pw[2,i]]]),
                                                  ncol=accData$nloci))))
    }
    # define true stat res obj
    pw_glb <- matrix(rep(0, (2 * (ncol(pw)))), ncol = 2)
    
    # true stat input
    trueStatIn <- lapply(pw_data, function(x){
      list(infile = x, gp = gp, bootstrap = FALSE)
    })
    
    # calculate true stats
    if(parallel){
      clusterExport(cl, c("readGenepopX", "fstWC"), envir = environment())
      tparSapply <- function(...) t(parSapply(...))
      trueStat <- tparSapply(cl, trueStatIn, function(x){
        input <- readGenepopX(x)
        return(fstWC(input)[[2]][2:3])
      })
    } else {
      tsapply <- function(...) t(sapply(...))
      trueStat <- tsapply(cl, trueStatIn, function(x){
        input <- readGenepopX(x)
        return(fstWC(input)[[2]][2:3])
      })
    }
    
    fstBS <- function(x){
      fstIn <- readGenepopX(x)
      fstOut <- fstWC(fstIn)
      return(fstOut[[2]][2:3])
    }
    
    if(parallel){
      # bootstrap pairwise populations
      clusterExport(cl, c("pw_data", "gp", "fstWC", "bootstraps", 
                          "readGenepopX", "fstBS"), envir = environment())
      pwRES <- parLapply(cl, 1:ncol(pw), function(i){
        input <- list(infile = pw_data[[i]], gp = gp, bootstrap = TRUE)
        out <- replicate(bootstraps, fstBS(input))
        fstCI <- (sd(out[1, ])/sqrt(length(out[1,]))) * 1.96
        fitCI <- (sd(out[2, ])/sqrt(length(out[2,]))) * 1.96
        return(c(fst = fstCI, fit = fitCI))
      })
      stopCluster(cl)
    } else {
      # bootstrap pairwise populations
      pwRES <- lapply(1:ncol(pw), function(i){
        input <- list(infile = pw_data[[i]], gp = gp, bootstrap = TRUE)
        out <- replicate(bootstraps, fstBS(input))
        fstCI <- (sd(out[1, ])/sqrt(length(out[1,]))) * 1.96
        fitCI <- (sd(out[2, ])/sqrt(length(out[2,]))) * 1.96
        return(c(fst = fstCI, fit = fitCI))
      })
    }
    pwOut <- lapply(1:2, function(i){      
      lapply(1:ncol(pw), function(j){
        if(i == 1){
          return(round(c(trueStat[j,"Fst_WC"], 
                         trueStat[j,"Fst_WC"] - pwRES[[j]]["fst"],
                         trueStat[j,"Fst_WC"] + pwRES[[j]]["fst"]), 4))
        } else if(i == 2){
          return(round(c(trueStat[j,"Fit_WC"], 
                         trueStat[j,"Fit_WC"] - pwRES[[j]]["fit"],
                         trueStat[j,"Fit_WC"] + pwRES[[j]]["fit"]),4))
        }
      })
    })
    pwOut1 <- lapply(pwOut, function(x){
      out <- as.data.frame(do.call("rbind", x))
      colnames(out) <- c("actual", "lower_CI", "upper_CI")
      rownames(out) <- paste(accData$pop_names[pw[1,]], " v ",
                             accData$pop_names[pw[2,]], sep = "")
      return(out)
    })
    
    # write pw bootstrap results
    pwWriteOut <- rbind(c("actual", "lower_CI", "upper_CI"),
                        pwOut1[[1]],
                        c("--", "--", "--"),
                        c("","",""),
                        c("actual", "lower_CI", "upper_CI"),
                        pwOut1[[2]])
    pwNames <- c("Fst", rownames(pwOut1[[1]]),
                 "--", "", "Fit", rownames(pwOut1[[1]]))
    pwOut <- as.matrix(cbind(pwNames, pwWriteOut))
    dimnames(pwOut) <- NULL
    
    if (!is.null(outfile)){
      of = paste(getwd(), "/", outfile, "-[fstWC]", "/", sep = "")
      if(is.element("xlsx", installed.packages()[, 1])){
        # write data to excel
        # Load dependencies
        library("xlsx")
        # standard stats
        if(bs_locus){
          write.xlsx(pwOut, file = paste(of, "[fstWC].xlsx", sep = ""),
                     sheetName = "Pairwise_stats", col.names = FALSE,
                     row.names = FALSE, append = TRUE)
        } else {
          write.xlsx(pwOut, file = paste(of, "[fstWC].xlsx", sep = ""),
                     sheetName = "Pairwise_stats", col.names = FALSE,
                     row.names = FALSE, append = FALSE)
        }                   
      } else {
        # text file alternatives
        std <- file(paste(of, "Pairwise_stats-[fstWC].txt", sep = ""), "w")
        for(i in 1:nrow(pwOut)){
          cat(pwOut[i,], "\n", file = std, sep = "\t")
        }
        close(std)
      }
    }
    names(pwOut1) <- c("Fst", "Fit")
  }
  # return results to the R enviroment
  if(bs_locus && bs_pairwise){
    list(locus = locBSout,
         pairwise = pwOut1)
  } else if (bs_locus && !bs_pairwise){
    return(locus = locBSout)
  } else if (bs_pairwise && !bs_locus){
    return(pairwise = pwOut1)
  }
}
################################################################################
# END
################################################################################
#
#
#
#
#
################################################################################
# divRatio: calculates diversity standardised to yardstick popukation
################################################################################
#' @export
divRatio <- function(infile = NULL, outfile = NULL, gp = 3, pop_stats =  NULL, 
                     refPos = NULL, bootstraps = 1000,  parallel = FALSE) {
  popStats = pop_stats
  NBS = bootstraps
  # create a directory for output
  if(!is.null(outfile)){
    suppressWarnings(dir.create(path=paste(getwd(),"/", outfile,
                                           "-[diveRsity]","/",sep="")))
    of = paste(getwd(), "/", outfile, "-[diveRsity]", "/", sep = "")
    write_res <- is.element("xlsx", installed.packages()[, 1])
  }
  # read the allelic richness and heterozygosity functions
  data1 <- fileReader(infile)
  data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
  #raw_data<-data1
  npops<-length(which(toupper(data1[,1]) == "POP"))
  pop_pos<- c(which(toupper(data1[,1]) == "POP"), (nrow(data1)+1))
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  # Calculate the minimum sample size
  pop_sizes <- sapply(1:npops, function(i){
    pop_pos[(i+1)] - pop_pos[i]-1
  })
  #minSize <- min(pop_sizes) 
  pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
  nloci<- (pop_pos[1]-2)
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  pop_list <- lapply(1:npops, function(i){
    return(as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
                           2:(nloci+1)]))
  })  
  if (gp==3) {
    plMake<-function(x){
      out <- matrix(sprintf("%06g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "0000NA"] <- "    NA"
      }
      return(out)
    }
  } else if (gp==2) {
    plMake<-function(x){
      out <- matrix(sprintf("%04g",as.numeric(x)),
                    nrow = nrow(x), ncol = ncol(x))
      if (Sys.info()["sysname"] == "Darwin"){
        out[out == "00NA"] <- "  NA"
      }
      return(out)
    }
  }
  suppressWarnings(pop_list<-lapply(pop_list, plMake))
  # deal with missing data
  if (gp == 3){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
    }
  } else if (gp == 2){
    for(i in 1:npops){
      pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
    }
  }
  ##############################################################################
  # if only the refpop raw data is given
  if(npops == 1 && !is.null(pop_stats)){
    refPop <- pop_list[[1]]
    # read subject population stats
    trypopDF <- try(read.table(popStats, header = TRUE), silent = TRUE)
    if(is(trypopDF, "try-error")){ 
      message <- paste("[ERROR]",
                       "",
                       "There is a problem with 'pop_stats' file format",
                       "",
                       "See the package user manual for more details",
                       "of file format requirements",
                       "",
                       sep = "\n")
      cat(message)
      stop()
    } else {
      popDF <- read.table(popStats, header = TRUE)
    }
    # calculate refpop standard stats
    refalr <- rowMeans(replicate(NBS, AR(refPop)))
    refhe <- Hex(refPop)
    refStd <- data.frame(alr = refalr, hexp = refhe)
    refStd <- data.frame(pops = paste(pop_names, "-(ref)", sep = ""),
                         n = nrow(refPop),
                         alr = mean(refStd$alr, na.rm = TRUE),
                         alrse = sd(refStd$alr, na.rm = TRUE) / 
                           sqrt(length(na.omit(refStd$alr))),
                         he = mean(refStd$hexp, na.rm = TRUE),
                         hese = sd(refStd$hexp, na.rm = TRUE) / 
                           sqrt(length(na.omit(refStd$hexp)))
    )
    if(is.element("validloci", names(popDF))){
      refStd$validloci <- paste(loci_names, sep = "\t", collapse = "\t")
    }
    # extract subject population sizes
    popDF <- rbind(refStd, popDF)
    popSizes <- as.numeric(popDF$n)
    # extract valid locus information
    if(is.element("validloci", names(popDF))){
      vlocs <- as.character(popDF$validloci)
      vlocs <- sapply(vlocs, function(x){
        return(strsplit(x, split = "\\s+"))
      })
      validLoc <- lapply(vlocs, function(x){
        return(sapply(x, function(y){
          as.numeric(which(loci_names == y))
        }))
      })
    } else {
      validLoc <- lapply(1:length(popSizes), function(...){
        locs <- 1:nloci
        names(locs) <- loci_names 
        return(locs)
      })
    }    
    
    # calculate refPopStats based on each subject pop sample size
    ###########################################################################
    if(parallel){
      library("doParallel")
      cores <- detectCores()
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      clusterExport(cl, c("arHex", "NBS", "refPop", "popSizes", 
                          "gp", "nloci", "validLoc"), envir = environment())
      refPopStats <- parLapply(cl, seq_along(popSizes), function(i){
        inner <- list(ref = refPop[,validLoc[[i]]], size = popSizes[i],
                      gp = gp)
        outer <- replicate(NBS, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
      stopCluster(cl)
    } else {
      refPopStats <- lapply(1:length(popSizes), function(i){
        inner <- list(ref = refPop, size = popSizes[i],
                      gp = gp)
        outer <- replicate(NBS, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
    }
    ###########################################################################
    # calculate the means and s.e. for all refPopStats
    refs <- lapply(refPopStats, function(x){
      meanAlr <- mean(x$alr, na.rm = TRUE)
      seAlr <- sd(x$alr, na.rm = TRUE)/sqrt(length(na.omit(x$alr)))
      meanHexp <- mean(x$hexp, na.rm = TRUE)
      seHexp <- sd(x$hexp, na.rm = TRUE)/sqrt(length(na.omit(x$hexp)))
      list(alr = data.frame(mean = meanAlr, se = seAlr),
           hexp = data.frame(mean = meanHexp, se = seHexp))
    })
    
    # extract all subject pops info
    subs <- lapply(1:nrow(popDF), function(i){
      return(list(alr = data.frame(mean = popDF$alr[i], se = popDF$alrse[i]),
                  hexp = data.frame(mean = popDF$he[i], se = popDF$hese[i])))
    })
    ##############################################################################
    # calculate the ratio stats
    seRatCalc <- function(subs, refs){
      rat <- subs$mean/refs$mean
      seRat <- sqrt((rat^2) * (((subs$se/subs$mean)^2) + ((refs$se/refs$mean)^2)))
      return(seRat)
    }
    tsapply <- function(...) t(sapply(...))  
    divRatio <- tsapply(1:length(subs), function(i){
      # create a subset variable for refs
      alrRat <- subs[[i]]$alr$mean / refs[[i]]$alr$mean
      hexpRat <- subs[[i]]$hexp$mean / refs[[i]]$hexp$mean
      alrSErat <- seRatCalc(subs[[i]]$alr, refs[[i]]$alr) 
      hexpSErat <- seRatCalc(subs[[i]]$hexp, refs[[i]]$hexp)
      res <- c(pop = as.character(popDF$pops[i]),
               n = popSizes[i],
               alr = round(subs[[i]]$alr$mean, 4),
               alrSE = round(subs[[i]]$alr$se, 4),
               He = round(subs[[i]]$hexp$mean, 4),
               HeSE = round(subs[[i]]$hexp$se, 4),
               alrRatio = round(alrRat, 4), 
               alrSEratio = round(alrSErat, 4),
               heRatio = round(hexpRat, 4),
               heSEratio = round(hexpSErat, 4))
      return(res)
    })
    
    divRatio <- as.data.frame(divRatio)
  } else {
    ############################################################################
    # Subset pop_list into subject populations and reference population
    # reference population
    refPop <- pop_list[[refPos]]
    
    # Run AR and Hexpected function for each population other than the refpop
    # call them subject populations
    if(parallel){
      library("doParallel")
      cores <- detectCores()
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      clusterExport(cl, c("Hex", "AR", "NBS", "gp", "nloci", "pop_list"), 
                    envir = environment())
      subPopStats <- parLapply(cl, pop_list, function(x){
        # Calculate allelic richness
        # bootstrap first
        alrbs <- replicate(NBS, AR(x))
        # calculate the mean of the bootstraps per locus
        alr <- rowMeans(alrbs)
        # Calculate expected Het
        hex <- Hex(x)
        # create return obj
        return(data.frame(alr = alr, hexp = hex))
      })
    } else {
      subPopStats <- lapply(pop_list, function(x){
        # Calculate allelic richness
        # bootstrap first
        alrbs <- replicate(NBS, AR(x))
        # calculate the mean of the bootstraps per locus
        alr <- rowMeans(alrbs)
        # Calculate expected Het
        hex <- Hex(x)
        # create return obj
        return(data.frame(alr = alr, hexp = hex))
      })
    }
    
    # Check if any loci in each population is missing data
    validLocs <- lapply(subPopStats, function(x){
      which(!is.na(x[,1]))
    })
    # calculate the standardized alr and hex for the ref pop
    if(parallel){
      clusterExport(cl, c("arHex", "refPop", "NBS", "gp", "pop_sizes", 
                          "validLocs"), envir = environment())
      refPopStats <- parLapply(cl, seq_along(pop_list), function(i){
        inner <- list(ref = refPop[,validLocs[[i]]], size = pop_sizes[i],
                      gp = gp)
        outer <- replicate(NBS, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
      stopCluster(cl)
    } else {
      refPopStats <- lapply(seq_along(pop_list), function(i){
        inner <- list(ref = refPop, size = pop_sizes[i],
                      gp = gp)
        outer <- replicate(NBS, arHex(inner), simplify = FALSE)
        alrPre <- sapply(outer, function(x){
          return(x$alls)
        })
        alr <- rowMeans(alrPre)
        hexpPre <- sapply(outer, function(x){
          return(x$hexp)
        })
        hexp <- rowMeans(hexpPre)
        return(data.frame(alr = alr, hexp = hexp))
      })
    }
    # calculate the means and s.e. for all refPopStats
    refs <- lapply(refPopStats, function(x){
      meanAlr <- mean(x$alr, na.rm = TRUE)
      seAlr <- sd(x$alr, na.rm = TRUE)/sqrt(length(na.omit(x$alr)))
      meanHexp <- mean(x$hexp, na.rm = TRUE)
      seHexp <- sd(x$hexp, na.rm = TRUE)/sqrt(length(na.omit(x$hexp)))
      list(alr = data.frame(mean = meanAlr, se = seAlr),
           hexp = data.frame(mean = meanHexp, se = seHexp))
    })
    # calculate the means and s.e. for all subPopStats
    subs <- lapply(subPopStats, function(x){
      meanAlr <- mean(x$alr, na.rm = TRUE)
      seAlr <- sd(x$alr, na.rm = TRUE)/sqrt(length(na.omit(x$alr)))
      meanHexp <- mean(x$hexp, na.rm = TRUE)
      seHexp <- sd(x$hexp, na.rm = TRUE)/sqrt(length(na.omit(x$hexp)))
      list(alr = data.frame(mean = meanAlr, se = seAlr),
           hexp = data.frame(mean = meanHexp, se = seHexp))
    })
    ##############################################################################
    # calculate the ratio stats
    seRatCalc <- function(subs, refs){
      rat <- subs$mean/refs$mean
      seRat <- sqrt((rat^2) * (((subs$se/subs$mean)^2) + ((refs$se/refs$mean)^2)))
      return(seRat)
    }
    tsapply <- function(...) t(sapply(...))  
    divRatio <- tsapply(seq_along(pop_list), function(i){
      # create a subset variable for refs
      alrRat <- subs[[i]]$alr$mean / refs[[i]]$alr$mean
      hexpRat <- subs[[i]]$hexp$mean / refs[[i]]$hexp$mean
      alrSErat <- seRatCalc(subs[[i]]$alr, refs[[i]]$alr) 
      hexpSErat <- seRatCalc(subs[[i]]$hexp, refs[[i]]$hexp)
      res <- c(pop = pop_names[i],
               n = pop_sizes[i],
               alr = round(subs[[i]]$alr$mean, 4),
               alrSE = round(subs[[i]]$alr$se, 4),
               He = round(subs[[i]]$hexp$mean, 4),
               HeSE = round(subs[[i]]$hexp$se, 4),
               alrRatio = round(alrRat, 4), 
               alrSEratio = round(alrSErat, 4),
               heRatio = round(hexpRat, 4),
               heSEratio = round(hexpSErat, 4))
      return(res)
    })
    # add reference data to divRatio
    refPop <- divRatio[refPos,]
    divRatio <- divRatio[-refPos,]
    refPop[1] <- paste(refPop[1], "-(ref)", sep = "")
    divRatio <- as.data.frame(rbind(refPop, divRatio))
  }
  if(!is.null(outfile)){
    if(write_res){
      library("xlsx")
      # standard stats
      write.xlsx(divRatio, file = paste(of, "[divRatio].xlsx",sep = ""),
                 sheetName = "Diversity_ratios", col.names = TRUE,
                 row.names = FALSE, append = FALSE)
    } else {
      write.table(divRatio, "divRatio-out.txt", col.names = TRUE, 
                  row.names = FALSE, append = FALSE, sep = "\t", 
                  quote = FALSE)
    }
  }
  divRatio[,-1] <- apply(divRatio[,-1], 2, function(x){
    return(as.numeric(as.character(x)))
  })
  divRatio[,1] <- as.character(divRatio[,1])
  return(divRatio)
}
################################################################################
# end divRatio
################################################################################
#
#
#
#
#
#
################################################################################
# AR: calculates the number of allele per locus from pop_list (divRatio)
################################################################################
# This function accepts a list containing a single population sample from
# a standard pop_list object, usually the ref sample and an interger 
# representing the resample size used to bootstrap allelic richness
# Define Allelic richness function for a single population to
# be bootstrapped for a given sample size
AR <- function(x){
  if(length(x) == 2L){
    pl <- x$ref
    mSize <- x$size
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  } else {
    pl <- x
    mSize <- nrow(pl)
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  }
  nloci <- ncol(pop_list)
  gp = nchar(pop_list[1,1])/2
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss <- function(x){  # where x is object pop_list
      pl <- list()
      pl[[1]] <- matrix(substr(x, 1, 2), ncol = nloci)
      pl[[2]] <- matrix(substr(x, 3, 4), ncol = nloci)
      return(pl)
    }
  }
  pop_alleles <- pl_ss(pop_list)
  alln <- function(x){
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    })
  }
  allele_names <- alln(pop_alleles)
  Alls <- sapply(allele_names, function(x){
    length(x)
  })
  return(Alls)
}
################################################################################
# End AR
################################################################################
#
#
#
#
#
#
#
################################################################################
# Hex: calcuates expected heterozygosity from pop_list (divRatio)
################################################################################
# This function accepts a list containing a single population sample from
# a standard pop_list object, usually the ref sample and an interger 
# representing the resample size used to bootstrap expected heterozygosity
# Define a hexp function
Hex <- function(x){
  if(length(x) == 2L){
    pl <- x$ref
    mSize <- x$size
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ],ncol=ncol(x)))
    }
    pop_list <- bser(pl) # resample
  } else {
    pop_list <- x
  }
  nloci = ncol(pop_list)
  gp = nchar(pop_list[1,1])/2
  # split string genotypes
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3), ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6), ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,2), ncol=nloci)
      pl[[2]]<-matrix(substr(x,3,4), ncol=nloci)
      return(pl)
    }
  }
  pop_alleles <- pl_ss(pop_list)
  alln <- function(x){
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][, i], x[[2]][, i])), decreasing = FALSE))
    })
  }
  allele_names <- alln(pop_alleles)
  # Calculate expected He
  #if(npops == 1){
  #  loci_combi <- allele_names[,1]
  #} else {
  #  loci_combi <- apply(allele_names, 1, FUN = 'unlist')
  #}
  aaList <- function(x){
    return(sort(unique(x, decreasing = FALSE)))
  }
  # finter out unique alleles
  all_alleles <- lapply(allele_names, aaList)
  # Create allele frequency holders
  allele_freq <- lapply(1:ncol(pop_list), function(i){
    Nrow <- length(all_alleles[[i]])
    Ncol <- length(pop_list)
    mat <- matrix(rep(0,(Ncol * Nrow)), ncol = Ncol)
    rownames(mat) <- all_alleles[[i]]
    return(mat)
  })
  # rbind pop_alleles
  pa1 <- rbind(pop_alleles[[1]], pop_alleles[[2]])
  
  # Count alleles
  actabPre <- function(x){
    lapply(1:ncol(x), function(i){
      table(x[,i])
    })
  }
  actab <- actabPre(pa1) 
  # Count the number of individuals typed per locus per pop
  indtyppop1 <- function(x){
    apply(x, 2, function(y){
      length(na.omit(y))/2
    })
  }
  indtyppop <- indtyppop1(pa1)
  #calculate allele frequencies
  afCalcpop <- sapply(1:length(actab), function(i){
    actab[[i]]/(indtyppop[i] * 2)
  })
  # calculate heterozygosities
  Hexp <- sapply(afCalcpop, function(x){
    1 - (sum(x^2))
  })
  return(Hexp)
}
################################################################################
# End Hex
################################################################################
#
#
#
#
#
#
#
################################################################################
# arHex: calculates bootstrapped allelic richness and He (divRatio)
################################################################################
# This function will calculate bootstrapped allelic richness and expected
# heterozygosity for use in the Skrbinsek diversity standardization method
arHex <- function(x){
  gp = x$gp
  nloci = ncol(x$ref)
  if(length(x) == 4L){
    pl <- x$ref
    mSize <- x$size
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], 
                    ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  } else {
    pl <- x$ref
    mSize <- nrow(pl)
    bser<-function(x){
      return(matrix(x[sample(nrow(x), mSize, replace = TRUE), ], 
                    ncol = ncol(x)))
    }
    pop_list <- bser(pl) # resample
  }
  if (gp==3){
    pl_ss<-function(x){  # where x is object pop_list
      pl<-list()
      pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
      pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
      return(pl)
    }
  } else {
    pl_ss <- function(x){  # where x is object pop_list
      pl <- list()
      pl[[1]] <- matrix(substr(x, 1, 2), ncol = nloci)
      pl[[2]] <- matrix(substr(x, 3, 4), ncol = nloci)
      return(pl)
    }
  }
  pop_alleles <- pl_ss(pop_list)
  alln <- function(x){
    res <- sapply(1:ncol(x[[1]]), function(i){
      list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
    })
  }
  allele_names <- alln(pop_alleles)
  Alls <- sapply(allele_names, function(x){
    length(x)
  })
  
  #return(Alls)
  #############################################################################
  # Heterozygosity
  aaList <- function(x){
    return(sort(unique(x, decreasing = FALSE)))
  }
  # finter out unique alleles
  all_alleles <- lapply(allele_names, aaList)
  # Create allele frequency holders
  allele_freq <- lapply(1:ncol(pop_list), function(i){
    Nrow <- length(all_alleles[[i]])
    Ncol <- length(pop_list)
    mat <- matrix(rep(0,(Ncol * Nrow)), ncol = Ncol)
    rownames(mat) <- all_alleles[[i]]
    return(mat)
  })
  # rbind pop_alleles
  pa1 <- rbind(pop_alleles[[1]], pop_alleles[[2]])
  
  # Count alleles
  actabPre <- function(x){
    lapply(1:ncol(x), function(i){
      table(x[,i])
    })
  }
  actab <- actabPre(pa1) 
  # Count the number of individuals typed per locus per pop
  indtyppop1 <- function(x){
    apply(x, 2, function(y){
      length(na.omit(y))/2
    })
  }
  indtyppop <- indtyppop1(pa1)
  #calculate allele frequencies
  afCalcpop <- sapply(1:length(actab), function(i){
    actab[[i]]/(indtyppop[i] * 2)
  })
  # calculate heterozygosities
  Hexp <- sapply(afCalcpop, function(x){
    1 - (sum(x^2))
  })
  # return(Hexp)
  return(data.frame(alls = Alls, hexp = Hexp))
}
################################################################################
# End arHex
################################################################################
#
#
#
#
#
#
################################################################################
# bigDivPart - a wrapper function for the calculation of diff stats
################################################################################
#' @export
bigDivPart <- function(infile = NULL, outfile = NULL, WC_Fst = FALSE,
                       format = NULL){
  
  
  fstat = WC_Fst
  on = outfile
  if (!is.null(on) && format != "txt" && format != "xlsx") {
    stop("Please provide a valid output file format")
  }
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  dat <- fastScan(fname = infile)
  if (length(strsplit(dat[length(dat)], split = "\\s+")[[1]]) == 
        1) {
    dat <- dat[-(length(dat))]
  }
  rm(fastScan)
  z <- gc()
  rm(z)
  popLocation <- grep("^([[:space:]]*)POP([[:space:]]*)$", 
                      toupper(dat))
  pop_pos <- c(popLocation, (length(dat) + 1))
  loci_names <- as.vector(sapply(dat[2:(pop_pos[1] - 1)], function(x) {
    gsub(pattern = "\\s+", replacement = "", x)
  }))
  popSizes <- NULL
  for (i in 1:(length(pop_pos) - 1)) {
    popSizes[i] <- length((pop_pos[i] + 1):(pop_pos[(i + 
                                                       1)] - 1))
  }
  pops <- dat[-(c(1:(popLocation[1] - 1), popLocation))]
  popList <- lapply(seq_along(popSizes), function(i) {
    if (i == 1) {
      indx <- 1:popSizes[i]
    } else {
      indx <- (sum(popSizes[1:(i - 1)]) + 1):((sum(popSizes[1:(i - 1)])) +
                                                popSizes[i])
    }
    return(pops[indx])
  })
  npops <- length(popList)
  nloci <- length(loci_names)
  pop_sizes <- popSizes
  rm(dat, pops)
  z <- gc(reset = TRUE)
  rm(z)
  testStr <- strsplit(popList[[1]][1], split = "\\s+")[[1]]
  gpEst <- sapply(testStr, function(x) {
    if (is.character(x)) {
      nchar(x)/2
    } else {
      NA
    }
  })
  rm(testStr)
  gp <- as.numeric(names(sort(-table(gpEst)))[1])
  prePopList <- lapply(popList, function(x) {
    y <- array(data = NA, dim = c(length(x), (nloci + 1), 2))
    colnames(y) <- c("ind", loci_names)
    for (j in 1:length(x)) {
      data <- strsplit(x[j], split = "\\s+")[[1]]
      if (data[2] == ",") {
        data <- data[-2]
      }
      data[data == "NANA"] <- NA
      data[data == "0"] <- NA
      data[data == "000000"] <- NA
      data[data == "999999"] <- NA
      data[data == "-9-9"] <- NA
      data[data == "0000"] <- NA
      y[j, 2:(nloci + 1), 1] <- substr(data[2:(nloci + 1)], 1, gp)
      y[j, 2:(nloci + 1), 2] <- substr(data[2:(nloci + 1)], gp + 1, gp * 2)
      y[j, 1, 1] <- data[1]
      y[j, 1, 2] <- data[1]
    }
    return(y)
  })
  rm(popList)
  ind_names <- lapply(prePopList, function(x) {
    return(x[, 1, 1])
  })
  pop_names <- sapply(ind_names, function(x) {
    return(x[1])
  })
  nb <- bigPreDiv(prePopList, FALSE, nloci, npops, popSizes, 
                      fstat)
  stdOut <- data.frame(loci = c(loci_names, "Global"), 
                       H_st = c(nb$hst, NA), 
                       D_st = c(nb$dst, NA), 
                       G_st = c(nb$gst, nb$gst_all), 
                       G_hed_st = c(nb$gst_hedrick, nb$gst_all_hedrick), 
                       D_Jost = c(nb$djost, nb$djost_all))
  if (fstat) {
    estOut <- data.frame(loci = c(loci_names, "Global"), 
                         Harmonic_N = c(nb$locus_harmonic_N, NA), 
                         H_st_est = c(nb$hst_est, NA), 
                         D_st_est = c(nb$dst_est, NA), 
                         G_st_est = c(nb$gst_est, nb$gst_est_all), 
                         G_hed_st = c(nb$gst_est_hedrick, 
                                      nb$gst_est_all_hedrick), 
                         D_Jost = c(nb$djost_est, nb$djost_est_all), 
                         Fst_WC = nb$fstats[, 1], Fit_WC = nb$fstats[, 2])
  } else {
    estOut <- data.frame(loci = c(loci_names, "Global"), 
                         Harmonic_N = c(nb$locus_harmonic_N, NA), 
                         H_st_est = c(nb$hst_est, NA), 
                         D_st_est = c(nb$dst_est, NA), 
                         G_st_est = c(nb$gst_est, nb$gst_est_all), 
                         G_hed_st = c(nb$gst_est_hedrick, 
                                      nb$gst_est_all_hedrick), 
                         D_Jost = c(nb$djost_est, nb$djost_est_all))
  }
  if (!is.null(on)) {
    suppressWarnings(dir.create(path = paste(getwd(), "/", 
                                             on, "-[diveRsity]", "/", 
                                             sep = "")))
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
  }
  write_res <- is.element("xlsx", installed.packages()[, 1])
  if (!is.null(on)) {
    if (write_res && format == "xlsx") {
      require("xlsx")
      write.xlsx(stdOut, file = paste(of, "[bigDivPart].xlsx", 
                                      sep = ""), 
                 sheetName = "Standard_stats", col.names = TRUE, 
                 row.names = FALSE, append = FALSE)
      
      write.xlsx(estOut, file = paste(of, "[bigDivPart].xlsx", 
                                      sep = ""), 
                 sheetName = "Estimated_stats", col.names = TRUE, 
                 row.names = FALSE, append = TRUE)
    } else {
      std <- file(paste(of, "Standard-stats[bigDivPart].txt", 
                        sep = ""), "w")
      cat(paste(colnames(stdOut), sep = ""), "\n", sep = "\t", 
          file = std)
      stdOut <- as.matrix(stdOut)
      for (i in 1:nrow(stdOut)) {
        cat(stdOut[i, ], "\n", file = std, sep = "\t")
      }
      close(std)
      est <- file(paste(of, "Estimated-stats[bigDivPart].txt", 
                        sep = ""), "w")
      cat(paste(colnames(estOut), sep = ""), "\n", sep = "\t", 
          file = est)
      estOut <- as.matrix(estOut)
      for (i in 1:nrow(estOut)) {
        cat(estOut[i, ], "\n", file = est, sep = "\t")
      }
      close(est)
    }
  }
  list(standard = stdOut, estimates = estOut)
}
################################################################################
#  END - bigDivPart
################################################################################
#
#
#
#
#
#
#
################################################################################
# bigPreDiv - a function for the calculation of diff stats from big files
################################################################################
bigPreDiv <- function(prePopList, bs = FALSE, nloci, npops, 
                      popSizes, fstat){
  ps <- popSizes
  if (bs) {
    popList <- lapply(prePopList, function(x) {
      boot <- sample(1:length(x[, 1, 1]), replace = TRUE)
      return(x[boot, (2:(nloci + 1)), ])
    })
  } else {
    popList <- lapply(prePopList, function(x) {
      return(x[, (2:(nloci + 1)), ])
    })
  }
  indtyp <- lapply(popList, function(x) {
    apply(x, 2, function(y) {
      length(na.omit(y[, 1]))
    })
  })
  alls <- lapply(seq_along(popList), function(i) {
    apply(popList[[i]], 2, function(x) {
      return(unique(c(x[, 1], x[, 2])))
    })
  })
  all_alleles <- lapply(1:nloci, function(i) {
    alleles <- lapply(alls, function(x) {
      return(x[[i]])
    })
    return(sort(unique(unlist(alleles))))
  })
  #   obsAlls <- lapply(popList, function(x) {
  #     apply(x, 2, function(y) {
  #       als <- unique(c(na.omit(y[, 1]), na.omit(y[, 2])))
  #       counts <- sapply(als, function(z) {
  #         res <- length(which(y == z))
  #         return(res)
  #       })
  #     })
  #   })
  obsAlls <- lapply(popList, function(x) {
    lapply(1:ncol(x), function(i){
      alls <- c(x[,i,1], x[,i,2])
      return(table(alls))
    })
  })
  allele_freq <- lapply(1:nloci, function(i) {
    loc <- matrix(nrow = length(all_alleles[[i]]), ncol = npops)
    rownames(loc) <- all_alleles[[i]]
    for (j in 1:npops) {
      o <- obsAlls[[j]][[i]]
      n <- indtyp[[j]][i]
      loc[names(o), j] <- o/(2 * n)
    }
    loc[is.na(loc)] <- 0
    return(loc)
  })
  preLoc <- lapply(indtyp, function(x) {
    return(1/x)
  })
  loci_harm_N <- sapply(1:nloci, function(i) {
    loc <- sapply(1:npops, function(j) {
      return(preLoc[[j]][i])
    })
    return(npops/sum(loc))
  })
  loci_harm_N <- round(loci_harm_N, 2)
  indtypLoc <- lapply(1:nloci, function(i) {
    res <- sapply(1:npops, function(j) {
      return(indtyp[[j]][i])
    })
  })
  rm(indtyp)
  if (fstat) {
    badData <- sapply(indtypLoc, function(y) {
      is.element(0, y)
    })
    if (sum(badData) > 0) {
      nl <- nloci - (sum(badData))
    } else {
      nl <- nloci
    }
    gdData <- which(!badData)
    badData <- which(badData)
    all_genot <- matrix(data = NA, nrow = sum(ps),
                        ncol =  length(gdData))
    for (i in 1:npops) {
      if (i == 1) {
        res <- apply(popList[[i]], 2, function(y) {
          return(paste0(y[, 1], y[, 2]))
        })
        all_genot[1:ps[i],] <- res[, gdData]
        rm(res)
      } else {
        res <- apply(popList[[i]], 2, function(y) {
          return(paste0(y[, 1], y[, 2]))
        })
        all_genot[(sum(ps[1:(i - 1)]) + 1):sum(ps[1:i]), ] <- res[, gdData]
        rm(res)
      }
    }
    all_genot[all_genot == "NANA"] <- NA
    genoCount <- lapply(1:ncol(all_genot), function(i){
      table(all_genot[,i])
    })
    nameFormat <- function(x) {
      nms <- names(x)
      lgth <- nchar(nms[1])
      newNms <- sapply(nms, function(y) {
        paste(substr(y, 1, lgth/2), "/", substr(y, (lgth/2) + 
                                                  1, lgth), sep = "")
      })
      names(x) <- newNms
      return(x)
    }
    genoCount <- lapply(genoCount, nameFormat)
    h_sum <- list()
    for (i in 1:length(gdData)) {
      h_sum[[i]] <- vector()
      cnSplit <- strsplit(names(genoCount[[i]]), "/")
      for (j in 1:length(all_alleles[[gdData[i]]])) {
        het_id1 <- lapply(cnSplit, is.element, all_alleles[[gdData[i]]][j])
        het_id2 <- lapply(het_id1, sum)
        het_id1 <- which(het_id2 == 1)
        h_sum[[i]][j] <- sum(genoCount[[i]][het_id1])
      }
    }
    indtyp_tot <- lapply(indtypLoc, sum)
    kk_hsum <- lapply(1:ncol(all_genot), function(i) {
      list(h_sum[[i]], indtyp_tot[[gdData[i]]])
    })
    kk_hbar <- lapply(kk_hsum, function(x) {
      return(x[[1]]/x[[2]])
    })
    pdat <- lapply(1:length(all_genot[1, ]), function(i) {
      list(allele_freq[[gdData[i]]], indtypLoc[[gdData[i]]])
    })
    kk_p <- lapply(pdat, function(x) {
      if (is.null(x[[1]]) == FALSE) {
        apply(x[[1]], 1, function(y) {
          y * (2 * x[[2]])
        })
      }
    })
    res <- matrix(0, (nloci + 1), 2)
    colnames(res) <- c("Fst_WC", "Fit_WC")
    A <- vector()
    a <- vector()
    b <- vector()
    c <- vector()
    for (i in 1:length(gdData)) {
      kknbar <- indtyp_tot[[gdData[i]]]/npops
      kknC <- (indtyp_tot[[gdData[i]]] - sum(indtypLoc[[gdData[i]]]^2)/indtyp_tot[[gdData[i]]])/(npops - 1)
      kkptild <- kk_p[[i]]/(2 * indtypLoc[[gdData[i]]])
      kkptild[kkptild == "NaN"] <- NA
      kkpbar <- colSums(kk_p[[i]])/(2 * indtyp_tot[[gdData[i]]])
      kks2 <- colSums(indtypLoc[[gdData[i]]] * (kkptild - rep(kkpbar, each = npops))^2)/((npops - 1) * kknbar)
      kkA <- kkpbar * (1 - kkpbar) - (npops - 1) * kks2/npops
      kka <- kknbar * (kks2 - (kkA - (kk_hbar[[i]]/4))/(kknbar - 1))/kknC
      kkb <- (kknbar/(kknbar - 1))*(kkA-((2*kknbar-1)/(4*kknbar))*kk_hbar[[i]])
      #kkb <- kknbar * (kkA - (2 * (kknbar - 1)) * kk_hbar[[i]]/(4 * kknbar))/(kknbar - 1)
      kkc <- kk_hbar[[i]]/2
      A[i] <- sum(kkA, na.rm = TRUE)
      a[i] <- sum(kka, na.rm = TRUE)
      b[i] <- sum(kkb, na.rm = TRUE)
      c[i] <- sum(kkc, na.rm = TRUE)
      res[gdData[i], "Fst_WC"] <- round(sum(kka)/sum(kka + kkb + kkc), 4)
      res[gdData[i], "Fit_WC"] <- round(1 - sum(kkc)/sum(kka + kkb + kkc), 4)
    }
    res[res == "NaN"] <- NA
    res[res == 0] <- NA
    sumA <- sum(A, na.rm = TRUE)
    suma <- sum(a, na.rm = TRUE)
    sumb <- sum(b, na.rm = TRUE)
    sumc <- sum(c, na.rm = TRUE)
    res[(nloci + 1), "Fst_WC"] <- round(suma/(suma + sumb + sumc), 4)
    res[(nloci + 1), "Fit_WC"] <- round(1 - sumc/(suma + sumb + sumc), 4)
    z <- gc(reset = TRUE)
    rm(z)
    fst <- res
    rm(res)
  }
  ho <- lapply(popList, function(x) {
    apply(x, 2, function(y) {
      1 - (sum(na.omit(y[, 1] == y[, 2]))/length(na.omit(y[, 1])))
    })
  })
  he <- t(sapply(allele_freq, function(x) {
    apply(x, 2, function(y) {
      return(1 - sum(y^2))
    })
  }))
  mf <- lapply(allele_freq, function(x) {
    rowSums(x)/ncol(x)
  })
  ht <- sapply(mf, function(x) {
    1 - sum(x^2)
  })
  hs <- rowSums(he)/npops
  hs_est <- hs * ((2 * loci_harm_N)/((2 * loci_harm_N) - 1))
  ht_est <- ht + (hs_est/(2 * loci_harm_N * npops))
  ht_est[is.nan(ht_est)] <- NA
  hst <- round((ht - hs)/(1 - hs), 4)
  dst <- round(ht - hs, 4)
  gst <- round(dst/ht, 4)
  gst[is.nan(gst)] <- NA
  djost <- round((dst/(1 - hs)) * (npops/(npops - 1)), 4)
  djost[djost == 0] <- NA
  hst_est <- round((ht_est - hs_est)/(1 - hs_est), 4)
  dst_est <- round(ht_est - hs_est, 4)
  gst_est <- round(dst_est/ht_est, 4)
  gst_est[is.nan(gst_est)] <- NA
  gst_max <- ((npops - 1) * (1 - hs))/(npops - 1 + hs)
  gst_est_max <- (((npops - 1) * (1 - hs_est))/(npops - 1 + 
                                                  hs_est))
  gst_hedrick <- round(gst/gst_max, 4)
  gst_est_hedrick <- round(gst_est/gst_est_max, 4)
  gst_est_hedrick[gst_est_hedrick > 1] <- 1
  djost_est <- round((npops/(npops - 1)) * ((ht_est - hs_est)/(1 - hs_est)), 4)
  djost_est[djost_est == 0] <- NA
  ht_mean <- round(mean(ht, na.rm = TRUE), 4)
  hs_mean <- round(mean(hs), 4)
  gst_all <- round((ht_mean - hs_mean)/ht_mean, 4)
  gst_all_max <- round(((npops - 1) * (1 - hs_mean))/(npops - 
                                                        1 + hs_mean), 4)
  gst_all_hedrick <- round(gst_all/gst_all_max, 4)
  djost_all <- round(((ht_mean - hs_mean)/(1 - hs_mean)) * 
                       (npops/(npops - 1)), 4)
  hs_est_mean <- mean(hs_est, na.rm = TRUE)
  ht_est_mean <- mean(ht_est, na.rm = TRUE)
  gst_est_all <- round((ht_est_mean - hs_est_mean)/ht_est_mean, 
                       4)
  gst_est_all_max <- round((((npops - 1) * (1 - hs_est_mean))/(npops - 
                                                                 1 + hs_est_mean)), 4)
  gst_est_all_hedrick <- round(gst_est_all/gst_est_all_max, 
                               4)
  gst_est_all_hedrick[gst_est_all_hedrick > 1] <- 1
  if (nloci == 1) {
    djost_est_all <- round(djost_est, 4)
  } else {
    djost_est_all <- round(1/(1/mean(djost_est, na.rm = TRUE) + 
                                (var(djost_est, na.rm = TRUE) * (1/mean(djost_est, 
                                                                        na.rm = TRUE))^3)), 4)
  }
  djost_est[djost_est == 0] <- NaN
  djost[djost == 0] <- NaN
  if (fstat) {
    list(hst = hst, dst = dst, gst = gst, gst_hedrick = gst_hedrick, 
         djost = djost, locus_harmonic_N = loci_harm_N, hst_est = hst_est, 
         dst_est = dst_est, gst_est = gst_est, 
         gst_est_hedrick = gst_est_hedrick, 
         djost_est = djost_est, gst_all = gst_all, 
         gst_all_hedrick = gst_all_hedrick, 
         djost_all = djost_all, gst_est_all = gst_est_all, 
         gst_est_all_hedrick = gst_est_all_hedrick, 
         djost_est_all = djost_est_all, 
         fstats = fst)
  } else {
    list(hst = hst, dst = dst, gst = gst, gst_hedrick = gst_hedrick, 
         djost = djost, locus_harmonic_N = loci_harm_N, hst_est = hst_est, 
         dst_est = dst_est, gst_est = gst_est, 
         gst_est_hedrick = gst_est_hedrick, 
         djost_est = djost_est, gst_all = gst_all, 
         gst_all_hedrick = gst_all_hedrick, 
         djost_all = djost_all, gst_est_all = gst_est_all, 
         gst_est_all_hedrick = gst_est_all_hedrick, 
         djost_est_all = djost_est_all)
  }
}
################################################################################
# END - bigPreDiv
################################################################################
#
#
#
#
#
################################################################################
# arp2gen: arlequin file conversion to genepop
################################################################################
#' @export
arp2gen <- function(infile){
  # test if the file exists
  flForm <- strsplit(infile, split = "\\.")[[1]]
  if(length(flForm) > 2){
    stop("There were multiple '.' characters in your file name!")
  }
  tstfile <- paste(flForm[1], ".gen", sep = "")
  # define a fastscan function
  if(!file.exists(tstfile)){
    fastScan <- function(fname){
      s <- file.info(fname)$size 
      buf <-  readChar(fname, s, useBytes = TRUE)
      return(strsplit(buf, "\n", fixed = TRUE, 
                      useBytes = TRUE)[[1]])
    }
    
    # scan infile
    dat <- fastScan(infile)
    
    # strip needless whitespace
    dat <- gsub("^\\s+|\\s+$", "", dat)
    
    # some safeguards
    dataType <- grep("*datatype=*", tolower(dat))
    if(strsplit(dat[dataType], "=")[[1]][2] != "MICROSAT"){
      stop("Data are not in 'MICROSAT' format!")
    }
    
    # extract the relavant information (nloci, npops etc.)
    
    # missing data character
    missDataLine <- grep("*missingdata=*", tolower(dat))
    missData <- noquote(substr(dat[missDataLine],
                               nchar(dat[missDataLine]) - 1,
                               nchar(dat[missDataLine]) - 1))
    
    # samples sizes
    sampSizeLine <- grep("*samplesize=*", tolower(dat))
    if(length(sampSizeLine) > 1){
      sampNpos <- sapply(sampSizeLine, function(i){
        return(regexpr("=", dat[i])[1])
      })
    }
    popSizes <- as.numeric(substr(dat[sampSizeLine],
                                  start = sampNpos+1,
                                  stop = nchar(dat[sampSizeLine])))
    
    # number of population samples
    npops <- length(popSizes)
    
    # number of loci
    sampStrt <- grep("*sampledata=*", tolower(dat))
    
    # adjust sample starts for possible white space
    strts <- sapply(sampStrt, function(x){
      if(dat[(x+1)] == ""){
        return(x + 2)
      } else {
        return(x + 1)
      }
    })
    
    # define pop ends
    ends <- strts + ((popSizes * 2) - 1)
    
    nloci <- length(strsplit(dat[strts[1]], split = "\\s+")[[1]]) - 2
    
    # extract genotypes
    popGeno <- lapply(seq_along(strts), function(i){
      return(dat[strts[i]:ends[i]])
    })
    
    # check that popsizes are consistent
    popSzcheck <- sapply(popGeno, function(x) length(x)/2)
    if(!all(identical(popSzcheck, popSizes))){
      stop("Failed! Please make sure that your file is formatted correctly.")
    }
    
    # create a vector of odd indexes for each pop
    popIdx <- lapply(popGeno, function(x){
      return(seq(1, length(x), 2))
    })
    
    # paste alleles together
    popList <- lapply(seq_along(popGeno), function(i){
      al1 <- matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], 
                                    split = "\\s+")), nrow = popSizes[i],
                    byrow = TRUE)[,-(1:2)]
      al2 <- matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] + 1)], 
                                    split = "\\s+")), nrow = popSizes[i],
                    byrow = TRUE)
      tst <- matrix(paste(al1, al2, sep = ""), nrow = popSizes[i])
      tst <- cbind(paste(rep("pop", nrow(tst)), i, " ,", sep = ""), tst)
      # tidy up
      rm(al1, al2)
      z <- gc()
      rm(z)
      # replace missing data with genepop format
      if(nchar(tst[1,2]) == 4){
        tst[tst == paste(missData, missData, sep = "")] <- "0000"
      } else {
        tst[tst == paste(missData, missData, sep = "")] <- "000000"
      }
      out <- apply(tst, 1, function(x){
        return(paste(x, collapse = "\t"))
      })
      out <- c("POP", out)
      #     out <- rep(NA, nrow(tst))
      #     for(j in 1:nrow(tst)){
      #       out[j] <- paste(tst[j,], collapse = "\t")
      #     }
      # tidy up
      rm(tst)
      z <- gc()
      rm(z)
      return(out)
    })
    
    # A genepop file can not be written easily
    
    # Generate the outfile name
    outfile <- strsplit(infile, "\\.")[[1]]
    if(length(outfile) >= 2){
      outfile <- paste(outfile[-length(outfile)], collapse = ".")
    } else {
      outfile <- outfile[1]
    }
    
    # construct the file
    loci <- paste("locus", 1:nloci, sep = "")
    loci <- c(paste(outfile, "_gen_converted", sep = ""), loci)
    
    # outfile object
    of <- c(loci, unlist(popList))
    
    # define a file connection
    out <- file(paste(outfile, ".gen", sep = ""), "w")
    for(i in 1:length(of)){
      cat(of[i], "\n", file = out, sep = "")
    }
    close(out)
    return(TRUE)
  } else {
    return(NULL)
  }
}
################################################################################
# END - arp2gen
################################################################################
#
#
#
#
#
################################################################################
# divMigrate: an experimental function for detecting directional differentiation 
################################################################################
# A function to calculate pairwise directional differentiation
# a presented in the paper 'Directional genetic differentiation and
# asymmetric migration Lisa Sundqvist, Martin Zackrisson & David Kleinhans,
# 2013, arXiv pre-print (http://arxiv.org/abs/1304.0118)'
#' #@export
# divMigrate <- function(infile = NULL, stat = "d_jost"){
#   # check file format
#   cat("Caution! The method used in this function is still under development. \n")
# #   flForm <- strsplit(infile, split = "\\.")[[1]]
# #   ext <- flForm[[length(flForm)]]
# #   if(ext == "arp"){
# #     arp2gen(infile)
# #     cat("Arlequin file converted to genepop format!")
# #     infile <- paste(flForm[1], ".gen", sep = "")
# #   }
#   dat <- fileReader(infile)
#   rownames(dat) <- NULL
#   dat <- as.matrix(dat)
#   # determine genepop format
#   p1 <- which(toupper(dat[,1]) == "POP")[1] + 1
#   gp <- as.numeric(names(sort(-table(sapply(dat[p1, - 1], nchar)/2)))[1])
#   dat <- as.data.frame(dat)
#   rawData <- readGenepop(dat, gp = gp)
#   npops <- rawData$npops
#   nloci <- rawData$nloci
#   # generate pairwise hypothetical matrices (use allele_freq)
#   pw <- combn(npops, 2)
#   # calculate ht and hs
#   hths <- lapply(rawData$allele_freq, pwDivCalc, pw = pw, npops = npops)
#   # seperate ht and hs matrices
#   ht <- lapply(hths, "[[", 1)
#   hs <- lapply(hths, "[[", 2)
#   # tidy
#   rm(hths)
#   z <- gc()
#   rm(z)
#   # find the mean (use the Reduce function)
#   ht_mean <- Reduce(`+`, ht)/nloci
#   hs_mean <- Reduce(`+`, hs)/nloci
#   # calculate Dst
#   dst <- ht_mean - hs_mean
#   # calculate gst
#   gst <- dst/ht_mean
#   gst[gst < 0.0 | is.na(gst)] <- 0
#   # calculate D(Jost)
#   d_jost <- (2*(dst))/(1-hs_mean)
#   # calculate relative migration from d_jost
#   d_mig <- (1 - d_jost)/d_jost
#   # replace missing and negative values with 0
#   d_mig[d_mig < 0 | is.na(d_mig)] <- 0
#   dimnames(d_mig) <- list(rawData$pop_names,
#                           rawData$pop_names)
#   # standardize
#   d_mig <- d_mig/max(d_mig, na.rm = TRUE)
#   # test gst migration rate
#   gst_mig <- 0.5 * ((1/gst) - 1)
#   # fix inf
#   gst_mig[is.infinite(gst_mig)] <- 0
#   # standardise
#   gst_mig <- gst_mig/max(gst_mig, na.rm = TRUE)
#   # replace missing and negative values with 0
#   gst_mig[gst_mig < 0 | is.na(gst_mig)] <- 0
#   dimnames(gst_mig) <- list(rawData$pop_names,
#                             rawData$pop_names)
#   # test plot
#   #library("qgraph")
#   if(length(stat) == 2){
#     par(mfrow = c(2, 1 ))
#     qgraph(gst_mig, posCol = "black",
#            nodeNames = rawData$pop_names)
#     title(expression("G"["st"]))
#     qgraph(d_mig, posCol = "black",
#            nodeNames = rawData$pop_names)
#     title(expression("D"["Jost"]))
#     par(mfrow = c(1,1))
#   } else if(stat == "gst"){
#     qgraph(gst_mig, posCol = "black",
#            nodeNames = rawData$pop_names)
#     title(expression("G"["st"]))
#   } else if(stat == "d_jost"){
#     qgraph(d_mig, posCol = "black",
#            nodeNames = rawData$pop_names)
#     title(expression("D"["Jost"]))
#   }
#   list(D_mig =d_mig,
#        Gst_mig = gst_mig)
# }
################################################################################
# END - divMigrate
################################################################################
#
#
#
#
#
################################################################################
# pwDivCalc: a small function for calculating pairwise ht and hs 
################################################################################
pwDivCalc <- function(x, pw, npops){
  ht <- matrix(ncol = npops, nrow = npops)
  hs <- matrix(ncol = npops, nrow = npops)
  for(i in 1:ncol(pw)){
    gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
    f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
    ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
    ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
    hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
    hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
  }
  ht[is.nan(ht)] <- 0
  hs[is.nan(hs)] <- 0
  list(ht = ht, 
       hs = hs)
}
################################################################################
# END - pwDivCalc
################################################################################
# divPart development version
# includes improved performance for pairwise calculations

# Kevin Keenan 2013

# divPart, a wrapper function for the calculation of differentiation stats.
#' @export
fastDivPart <- function(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
                        WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE, 
                        bootstraps = 0, plot = FALSE, parallel = FALSE){
  
  ############################ Argument definitions ############################
  # define arguments for testing
  #   D <- "MyData_GP2D.gen"
  #   on <- NULL
  #   fst <- T
  #   bstrps <- 10
  #   bsls <- T
  #   bspw <- T
  #   plt <- F
  #   para <- T
  #   pWise <- T
  #   gp = 2
  # define
  D <- infile
  on <- outfile
  fst <- WC_Fst
  bstrps <- bootstraps
  bsls <- bs_locus
  bspw <- bs_pairwise
  plt <- plot
  para <- parallel
  pWise <- pairwise
  gp = gp
  ##############################################################################
  if(bsls==T && bstrps<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else if (bspw==T && bstrps<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else {
    #Use pre.div to calculate the standard global and locus stats
    accDat <- pre.divLowMemory(x <- list(infile = D,
                                         gp = gp,
                                         bootstrap = FALSE,
                                         locs = TRUE,
                                         fst = fst,
                                         min = FALSE))
    # create a directory for output
    if(!is.null(on)){
      suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                             "-[diveRsity]","/",sep="")))
    }
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
    wd <- getwd()
    write_res <- is.element("xlsx", installed.packages()[, 1])
    plot_res <- is.element("sendplot", installed.packages()[, 1])
    
    para_pack_inst<-is.element(c("parallel","doParallel","foreach","iterators"),
                               installed.packages()[,1])
    
    if(plt == TRUE && is.null(on)){
      writeWarn <- paste("", "[NOTE]",
                         "Your results can't be plotted as you have not",
                         "provided an argument for 'outfile'.",
                         "Analysis completed", sep="\n")
      cat(noquote(writeWarn))
    }
    para_pack <- all(para_pack_inst)
    if(write_res == FALSE){
      Warning1<-{paste(" "," ",
                       "[NOTE]",
                       "___________________________________________________________",
                       "Please install the package 'xlsx' if you would like your", 
                       "results written to an Excel workbook.",
                       "Alternatively, your result will automatically be written",
                       "to .txt files.",
                       "___________________________________________________________",
                       "To install 'xlsx' use:",
                       "> install.packages('xlsx', dependencies=TRUE)",
                       "See:",
                       "> ?install.packages - for usage details.",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning1))
    } 
    if(plot_res==F && plt==T){
      Warning2<-{paste(" "," "," ",
                       "[NOTE]  ",
                       "___________________________________________________________",
                       "Please install the package 'sendplot' to plot your results.",
                       "Use:",
                       "> install.packages('sendplot', dependencies = TRUE)",
                       "See:",
                       "> ?install.packages - for usage details",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning2))
    }
    if(fst == TRUE){
      namer<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
               "D_Jost_est","Fst_WC","Fit_WC")
    } else {
      namer<-c("Gst","G_hed_st","D_Jost","Gst_est","G_hed_st_est",
               "D_Jost_est")
    }
    
    ############################################################################
    # output file multilocus stats vector 
    # pre output table for global locus stats
    
    #standard
    pre_ot1 <- cbind(accDat$locus_names, round(as.numeric(accDat$hst), 4),
                     round(as.numeric(accDat$dst), 4),
                     round(as.numeric(accDat$gst), 4),
                     round(as.numeric(accDat$gst_hedrick), 4),
                     round(as.numeric(accDat$djost), 4))
    # Add global multi locus stats to output table
    ot1 <- rbind(pre_ot1, c("Global", "", "", accDat$gst_all, 
                            accDat$gst_all_hedrick, 
                            accDat$djost_all))
    colnames(ot1) <- c("loci", "H_st", "D_st", "G_st", "G_hed_st", "D_jost")
    #Estimated
    pre_ot2 <- cbind(accDat$locus_names,
                     round(as.numeric(accDat$locus_harmonic_N),4),
                     round(as.numeric(accDat$hst_est),4),
                     round(as.numeric(accDat$dst_est),4),
                     round(as.numeric(accDat$gst_est),4),
                     round(as.numeric(accDat$gst_est_hedrick),4),
                     round(as.numeric(accDat$djost_est),4))
    
    ot2 <- rbind(pre_ot2, c("Global", "", "", "", accDat$gst_est_all, 
                            accDat$gst_est_all_hedrick, 
                            accDat$djost_est_all))
    colnames(ot2) <- c("loci", "Harmonic_N", "H_st_est", "D_st_est",
                       "G_st_est", "G_hed_st_est", "D_Jost_est")
    if(fst == TRUE){
      ot2 <- cbind(ot2, accDat$fstats[, 2:3])
    }
    if(fst == TRUE){
      plot_data321 <- c("Overall","","","",accDat$gst_est_all,
                        accDat$gst_est_all_hedrick,
                        accDat$djost_est_all,
                        as.numeric(accDat$fstats["All",2]))
      
    } else {
      plot_data321<-c("Overall","","","",accDat$gst_est_all,
                      accDat$gst_est_all_hedrick,
                      accDat$djost_est_all)
    }
    if (!is.null(on)){
      if(write_res==TRUE){
        # write data to excel
        # Load dependencies
        require("xlsx")
        # standard stats
        write.xlsx(ot1,file=paste(of,"[fastDivPart].xlsx",sep=""),
                   sheetName="Standard_stats",col.names=T,
                   row.names=F,append=F)
        # Estimated stats
        write.xlsx(ot2,file=paste(of,"[fastDivPart].xlsx",sep=""),
                   sheetName="Estimated_stats",col.names=T,
                   row.names=F,append=T)
      } else {
        # text file alternatives
        std<-file(paste(of,"Standard-stats[fastDivPart].txt",sep=""), "w")
        cat(paste(colnames(ot1),sep=""),"\n",sep="\t",file=std)
        for(i in 1:nrow(ot1)){
          cat(ot1[i,],"\n",file=std,sep="\t")
        }
        close(std)
        est<-file(paste(of,"Estimated-stats[fastDivPart].txt",sep=""),"w")
        cat(paste(colnames(ot2),sep=""),"\n",sep="\t",file=est)
        for(i in 1:nrow(ot2)){
          cat(ot2[i,],"\n",file=est,sep="\t")
        }
        close(est)
      }
    }
    ot1out<-ot1[,-1]
    ot2out<-ot2[,-1]
    
    ot1out<-matrix(as.numeric(ot1[,2:6]),ncol=5)
    rownames(ot1out)<-ot1[,1]
    colnames(ot1out)<-colnames(ot1)[-1]
    
    ot2out<-matrix(as.numeric(ot2[,-1]),ncol=(ncol(ot2)-1))
    rownames(ot2out)<-ot2[,1]
    colnames(ot2out)<-colnames(ot2)[-1]
    if (para && !para_pack){
      Warning3<-{paste(" "," ",
                       "[NOTE]",
                       "___________________________________________________________",
                       "Please make sure the packages 'parallel', 'doParallel',",
                       "'foreach' and 'iterators' are installed. These are required",
                       " to run your analysis in parallel.",
                       "Your analysis will be run sequentially!",
                       "___________________________________________________________",
                       "To install these use:",
                       "> install.packages()",
                       "See:",
                       "> ?install.packages - for usage details.",
                       "___________________________________________________________",
                       sep="\n")
      }
      cat(noquote(Warning3))
    }
    
    ############################################################################
    ############################ Bootstrapper ##################################
    ############################################################################
    # Used only if bootstraps is greater than zero
    if(bsls == TRUE){
      
      if (para && para_pack) {
        
        if (para && para_pack) {
          #count cores
          library("doParallel")
          cores <- detectCores()
          cl<-makeCluster(cores)
          registerDoParallel(cl)
        }
        
        #vectorize prallele#
        gp_inls <- list(infile = D, gp = gp,
                        bootstrap = TRUE, 
                        locs = TRUE, fst = fst)
        # silence for memory efficiency
        #gp_in <- list()
        #for(i in 1:bstrps){
        #  gp_in[[i]] <- gp_inls
        #}
        
        # calculate stats from readGenepopX objects
        # export objects for parallel
        clusterExport(cl, c("gp_inls", "pre.divLowMemory"), 
                      envir = environment())
        # run parallel code
        bs_loc <- parLapply(cl, 1:bstrps, function(...){
          pre.divLowMemory(gp_inls)
        })
        # close the cluster connection
        stopCluster(cl)
        
        
        #vectorize data extraction#
        if(fst==TRUE){
          bs_glb <- do.call("rbind", lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all, 4),
              round(bs_loc[[x]]$gst_all_hedrick, 4),
              round(bs_loc[[x]]$djost_all, 4),
              round(bs_loc[[x]]$gst_est_all, 4),
              round(bs_loc[[x]]$gst_est_all_hedrick, 4),
              round(bs_loc[[x]]$djost_est_all, 4),
              as.numeric(bs_loc[[x]]$fstats["All", 2:3]))
          }))
        } else {
          bs_glb <- do.call("rbind", lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all, 4),
              round(bs_loc[[x]]$gst_all_hedrick, 4),
              round(bs_loc[[x]]$djost_all, 4),
              round(bs_loc[[x]]$gst_est_all, 4),
              round(bs_loc[[x]]$gst_est_all_hedrick, 4),
              round(bs_loc[[x]]$djost_est_all, 4))
          }))
        }
        bs_std <- lapply(1:accDat$nloci, function(x){
          do.call("rbind", lapply(1:length(bs_loc), function(y){
            c(round(bs_loc[[y]]$gst[x], 4),
              round(bs_loc[[y]]$gst_hedrick[x], 4),
              round(bs_loc[[y]]$djost[x], 4))
          }))
        })
        if(fst==TRUE){
          bs_est <- lapply(1:accDat$nloci, function(x){
            do.call("rbind", lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x], 4),
                round(bs_loc[[y]]$gst_est_hedrick[x], 4),
                round(bs_loc[[y]]$djost_est[x], 4),
                as.numeric(bs_loc[[y]]$fstats[x, 2:3]))
            }))
          })
        } else {
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4))
            }))
          })
        }
        rm(bs_loc)                  ###
        z<-gc(reset=T)                ### tidy up
        rm(z)                       ###
        
      } else {
        #vectorize non-parallel#
        
        gp_inls <- list(infile = D,
                        gp = gp,
                        bootstrap = TRUE, 
                        locs = TRUE, 
                        fst = fst)
        #gp_in<-list()
        #for(i in 1:bstrps){
        # gp_in[[i]]<-gp_inls
        #}
        # calculate stats from readGenepopX objects
        bs_loc <- lapply(1:bstrps, function(...){
          pre.divLowMemory(gp_inls)
        })
        
        
        if(fst==TRUE){
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4),
              as.numeric(bs_loc[[x]]$fstats[(accDat$nloci+1),2:3]))
          }))
        }else{
          bs_glb<-do.call("rbind",lapply(1:bstrps, function(x){
            c(round(bs_loc[[x]]$gst_all,4),
              round(bs_loc[[x]]$gst_all_hedrick,4),
              round(bs_loc[[x]]$djost_all,4),
              round(bs_loc[[x]]$gst_est_all,4),
              round(bs_loc[[x]]$gst_est_all_hedrick,4),
              round(bs_loc[[x]]$djost_est_all,4))
          }))
        }
        bs_std<-lapply(1:accDat$nloci, function(x){
          do.call("rbind",lapply(1:length(bs_loc), function(y){
            c(round(bs_loc[[y]]$gst[x],4),
              round(bs_loc[[y]]$gst_hedrick[x],4),
              round(bs_loc[[y]]$djost[x],4))}))
        })
        if(fst==TRUE){
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4),
                as.numeric(bs_loc[[y]]$fstats[x,2:3]))
            }))
          })
        } else {
          bs_est<-lapply(1:accDat$nloci, function(x){
            do.call("rbind",lapply(1:length(bs_loc), function(y){
              c(round(bs_loc[[y]]$gst_est[x],4),
                round(bs_loc[[y]]$gst_est_hedrick[x],4),
                round(bs_loc[[y]]$djost_est[x],4))
            }))
          })
        }
        rm(bs_loc)
        z<-gc(reset=T)
        rm(z)
        
      }
      
      
      #vectorize#
      if(fst == TRUE){
        bs_res <- lapply(1:8, function(x){
          matrix(ncol = 3, nrow = (accDat$nloci+1))
        })
      } else {
        bs_res<-lapply(1:6,function(x){matrix(ncol=3, nrow=(accDat$nloci+1))})
      }
      bs_join<-cbind(bs_std, bs_est)
      bs_cis <- apply(bs_join, 1, function(x){
        res <- lapply(x, function(y){
          apply(y, 2, function(z){
            ci <- as.vector(quantile(z, probs = c(0.025, 0.975), na.rm = TRUE))
            means <- mean(z, na.rm = TRUE)
            
            return(c(means, ci))
          })
        })
        ciM <- c(res$bs_std[1,], res$bs_est[1,])
        lci <- c(res$bs_std[2,], res$bs_est[2,])
        uci <- c(res$bs_std[3,], res$bs_est[3,])
        list(mu = ciM,
             lci = lci,
             uci = uci)
      })
      mu <- t(sapply(1:length(bs_cis), function(i){
        return(bs_cis[[i]]$mu)
      }))
      lci <- t(sapply(1:length(bs_cis), function(i){
        return(bs_cis[[i]]$lci)
      }))
      uci <- t(sapply(1:length(bs_cis), function(i){
        return(bs_cis[[i]]$uci)
      }))
      # calculate ci for global
      glb_mu <- apply(bs_glb, 2, function(x){
        return(mean(x, na.rm = TRUE))
      })
      glb_lci <- apply(bs_glb, 2, function(x){
        return(quantile(x, probs = 0.025, na.rm = TRUE))
      })
      glb_uci <- apply(bs_glb, 2, function(x){
        return(quantile(x, probs = 0.975, na.rm = TRUE))
      })
      # add glb ci to mu,  uci and lci
      mu <- rbind(mu, glb_mu)
      lci <- rbind(lci, glb_lci)
      uci <- rbind(uci, glb_uci)
      #ciCalc <- function(x){
      #  res <- lapply(x, function(y){
      #    apply(y, 2, function(z){
      #      return(quantile(z, probs = c(0.025, 0.975)))
      #    })
      #  })
      #  return(res)
      #}
      #ci <- function(x){
      #  (sd(na.omit(x))/sqrt(length(na.omit(x)))) * 1.96
      #}
      #bs_cis <- t(apply(bs_join, 1, ciCalc))
      #bs_cis<-rbind(bs_cis, apply(bs_glb, 2, ci))
      if(fst==TRUE){
        for(i in 1:8){
          bs_res[[i]][,1] <- round(mu[,i], 4)
          bs_res[[i]][,2] <- round(lci[,i], 4)
          bs_res[[i]][,3] <- round(uci[,i], 4)
          bs_res[[i]][is.na(bs_res[[i]])] <- 0
        }
      } else {
        for(i in 1:6){
          bs_res[[i]][,1] <- round(mu[,i], 4)
          bs_res[[i]][,2] <- round(lci[,i], 4)
          bs_res[[i]][,3] <- round(uci[,i], 4)
          bs_res[[i]][is.na(bs_res[[i]])] <- 0
        }
      }
      
      names(bs_res) <- namer
      
      bs_res1 <- bs_res
      if(fst){
        for(i in 1:8){
          dimnames(bs_res1[[i]])<-list(c(accDat$locus_names, "global"),
                                       c("Mean","Lower_CI", "Upper_CI"))
        }
      } else {
        for(i in 1:6){
          dimnames(bs_res1[[i]])<-list(c(accDat$locus_names,"global"),
                                       c("Mean","Lower_CI","Upper_CI"))
        }
      }
      # bs results output object header
      hdr <- matrix(c("locus", "Mean", "Lower_95%CI", "Upper_95%CI"), 
                    ncol=4)
      bs_out <- matrix(rbind(hdr, c(names(bs_res)[1], "", "", ""),
                             cbind(c(accDat$locus_names, "Overall"),
                                   bs_res[[1]])), ncol = 4)
      
      if(fst){
        for(i in 2:8){
          bs_out <- matrix(rbind(bs_out, c(names(bs_res)[i], "", "", ""),
                                 cbind(c(accDat$locus_names, "global"),
                                       bs_res[[i]])), ncol = 4)
        }
      } else {
        for(i in 2:6){
          bs_out<-matrix(rbind(bs_out,c(names(bs_res)[i],"","",""),
                               cbind(c(accDat$locus_names,"Global"),
                                     bs_res[[i]])),ncol=4)
        }
      }
      if(!is.null(on)){
        if(write_res==TRUE){
          write.xlsx(bs_out,file=paste(of,"[fastDivPart].xlsx",sep=""),
                     sheetName="Locus_bootstrap",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          bts<-file(paste(of,"Locus-bootstrap[fastDivPart].txt",sep=""), "w")
          cat(paste(colnames(bs_out),sep=""),"\n",sep="\t",file=bts)
          for(i in 1:nrow(bs_out)){
            cat(bs_out[i,],"\n",file=bts,sep="\t")
          }
          close(bts)
        }
      }
    }
    zzz<-gc()
    rm(zzz)
    if(plot_res==TRUE && plt==TRUE && bsls==TRUE){
      
      #vectorize#
      sorter<-function(x){
        z<-order(x[1:accDat$nloci,1],decreasing=F)
        #if(length(z) >= 200){
        #  z<-z[(length(z)-150):length(z)]
        #}
        return(z)
      }
      lso123<-lapply(bs_res, sorter)
      
      #
      names(lso123)<-namer
      plot.call_loci<-list()
      plot.extras_loci<-list()
      xy.labels_loci<-list()
      y.pos_loci<-list()
      x.pos_loci=1:accDat$nloci
      direct=of
      fn_pre_loci<-list()
      #Plot Gst_Nei
      plot.call_loci[[1]]=c("plot(bs_res[[4]][lso123[[4]],1],
                            ylim=c(0,(max(bs_res[[4]][,3])+
                            min(bs_res[[4]][,3]))),xaxt='n',
                            ylab=names(bs_res)[4],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")
      
      plot.extras_loci[[1]]=c("points(bs_res[[4]][lso123[[4]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[4]][lso123[[4]],2],
                              1:accDat$nloci,bs_res[[4]][lso123[[4]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[4]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
      
      xy.labels_loci[[1]]=data.frame(Locus_name=accDat$locus_names[lso123[[4]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[4]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[4]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[4]],1],4))
      
      y.pos_loci[[1]]=bs_res[[4]][lso123[[4]],1]
      fn_pre_loci[[1]]<-names(bs_res)[4]
      
      
      
      # Plot Gst_Hedrick
      plot.call_loci[[2]]=c("plot(bs_res[[5]][lso123[[5]],1],
                            ylim=c(0,1),xaxt='n',ylab=names(bs_res)[5],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")
      
      plot.extras_loci[[2]]=c("points(bs_res[[5]][lso123[[5]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[5]][lso123[[5]],2],
                              1:accDat$nloci,bs_res[[5]][lso123[[5]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[5]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
      
      xy.labels_loci[[2]]=data.frame(Locus_name=accDat$locus_names[lso123[[5]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[5]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[5]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[5]],1],4))
      
      y.pos_loci[[2]]=bs_res[[5]][lso123[[5]],1]
      fn_pre_loci[[2]]<-names(bs_res)[5]
      
      
      # Plot D_jost
      plot.call_loci[[3]]=c("plot(bs_res[[6]][lso123[[6]],1],
                            ylim=c(0,1),xaxt='n',ylab=names(bs_res)[6],type='n',
                            xlab='Loci \n (Hover over a point to see locus data)',
                            cex.lab=1.5,cex.axis=1.3,las=1)")
      
      plot.extras_loci[[3]]=c("points(bs_res[[6]][lso123[[6]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:accDat$nloci,bs_res[[6]][lso123[[6]],2],
                              1:accDat$nloci,bs_res[[6]][lso123[[6]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=c(0,bs_res[[6]][(accDat$nloci+1),2]),
                              lwd=1,lty=c(1,2),col=c('black','red'))")
      
      xy.labels_loci[[3]]=data.frame(Locus_name=accDat$locus_names[lso123[[6]]],
                                     Gst_Nei=round(bs_res[[4]][lso123[[6]],1],4),
                                     Gst_Hedrick=round(bs_res[[5]][lso123[[6]],1],4),
                                     D_jost=round(bs_res[[6]][lso123[[6]],1],4))
      
      y.pos_loci[[3]]=bs_res[[6]][lso123[[6]],1]
      fn_pre_loci[[3]]<-names(bs_res)[6]
      
      #plot(Fst)
      if(fst==TRUE){
        plot.call_loci[[4]]=c("plot(bs_res[[8]][lso123[[8]],1],
                              ylim=c(0,(max(bs_res[[8]][,3])+
                              min(bs_res[[8]][,3]))),xaxt='n',
                              ylab=names(bs_res)[8],type='n',
                              xlab='Loci \n (Hover over a point to see locus data)',
                              cex.lab=1.5,cex.axis=1.3,las=1)")
        
        plot.extras_loci[[4]]=c("points(bs_res[[8]][lso123[[8]],1],
                                pch=15,col='black',cex=1);
                                arrows(1:accDat$nloci,bs_res[[8]][lso123[[8]],2],
                                1:accDat$nloci,bs_res[[8]][lso123[[8]],3],code=3,
                                angle=90,length=0.05,lwd=0.1);
                                abline(h=c(0,bs_res[[8]][(accDat$nloci+1),2]),
                                lwd=1,lty=c(1,2),col=c('black','red'))")
        
        xy.labels_loci[[4]]=data.frame(Locus_name=accDat$locus_names[lso123[[8]]],
                                       Gst_Nei=round(bs_res[[4]][lso123[[8]],1],4),
                                       Gst_Hedrick=round(bs_res[[5]][lso123[[8]],1],4),
                                       D_jost=round(bs_res[[6]][lso123[[8]],1],4),
                                       Fst_WC=round(bs_res[[8]][lso123[[8]],1],4))
        
        y.pos_loci[[4]]=bs_res[[8]][lso123[[8]],1]
        fn_pre_loci[[4]]<-names(bs_res)[8]
      }
    }
    ############################################################################
    ################################## Pairwise ################################
    ############################################################################
    # population pair combinations
    
    # define new functions
    ############################################################################
    ############################################################################
    # pwCalc
    ############################################################################
    # New optimised function for the calculation of pairwise statistics
    # Returns a 3D array where each 'slot' represents the pairwise matrix
    # for Gst_est, G_st_est_hed and D_jost_est respectively
    
    # Kevin Keenan
    # 2013
    
    pwCalc <- function(infile, fst,  bs = FALSE){
      
      
      #   # uncomment for testing
      #   infile <- "pw_test.txt"
      #   source("readGenepopX.R")
      #   # read pwBasicCalc function
      #   source("pwBasicCalc.R")
      # define baseline info
      dat <- readGenepopX(list(infile = infile,
                               bootstrap = bs))
      if(fst){
        # calculate all fst
        fstat <- pwFstWC(dat)
        # extract locus theta and variance components
        locTheta <- lapply(fstat, "[[", 1)
        # sum res
        aLoc <- Reduce(`+`, lapply(fstat, "[[", 2))
        bLoc <- Reduce(`+`, lapply(fstat, "[[", 3))
        cLoc <- Reduce(`+`, lapply(fstat, "[[", 4))
        # calculate pw Fst across loci
        pwTheta <- aLoc/(aLoc+bLoc+cLoc)
        # clean up
        rm(aLoc, bLoc, cLoc, fstat)
        z <- gc()
        rm(z)
      }
      # extract allele frequencies
      af <- dat$allele_freq
      # extract harmonic mean sample sizes
      
      # make space in RAM
      dat$allele_freq <- NULL
      z <- gc()
      rm(z)
      # extract npops and nloci
      npops <- dat$npops
      nloci <- dat$nloci
      # define pairwise relationships
      pw <- combn(dat$npops, 2)
      # generate pairwise locus harmonic mean sample sizes
      indtyp <- dat$indtyp
      pwHarm <- lapply(indtyp, pwHarmonic, pw = pw)
      
      
      # calculate pairwise ht and hs
      hths <- mapply(pwBasicCalc, af, pwHarm,
                     MoreArgs = list(pw = pw, npops = dat$npops),
                     SIMPLIFY = FALSE)
      # seperate ht and hs
      #   htLoc <- lapply(hths, "[[", 1)
      #   hsLoc <- lapply(hths, "[[", 2)
      # seperate ht_est and hs_est
      hsEstLoc <- lapply(hths, "[[", 1)
      htEstLoc <- lapply(hths, "[[", 2)
      
      # clean up
      rm(hths)
      z <- gc()
      rm(z)
      
      # Calculate locus stats
      # Standard locus stats
      # locus Gst
      #   gstLoc <- mapply(FUN = gstCalc, ht = htLoc, hs = hsLoc, 
      #                    SIMPLIFY = FALSE)
      #   # locus G'st
      #   gstHedLoc <- mapply(FUN = gstHedCalc, ht = htLoc, hs = hsLoc,
      #                       SIMPLIFY = FALSE)
      #   # locus D_jost
      #   dLoc <- mapply(FUN = djostCalc, ht = htLoc, hs = hsLoc,
      #                  SIMPLIFY = FALSE)
      
      # Estimated locus stats
      # locus Gst_est
      gstLocEst <- mapply(FUN = gstCalc, ht = htEstLoc, 
                          hs = hsEstLoc, 
                          SIMPLIFY = FALSE)
      # locus G'st_est
      gstHedLocEst <- mapply(FUN = gstHedCalc, ht = htEstLoc, 
                             hs = hsEstLoc,
                             SIMPLIFY = FALSE)
      # locus D_jost_est
      dLocEst <- mapply(FUN = djostCalc, ht = htEstLoc, 
                        hs = hsEstLoc,
                        SIMPLIFY = FALSE)
      
      #   # calculate mean ht and hs
      #   htMean <- Reduce(`+`, htLoc)/nloci
      #   hsMean <- Reduce(`+`, hsLoc)/nloci
      # calculate mean ht_est and hs_est
      htEstMean <- Reduce(`+`, htEstLoc)/nloci
      hsEstMean <- Reduce(`+`, hsEstLoc)/nloci
      
      # calculate standard stats (uncomment for loc stats)
      
      #   # overall dst
      #   dstAll <- htMean - hsMean
      #   # overall gst (Nei 1973)
      #   gstAll <- (dstAll)/htMean
      #   # overall max gst (Hedricks 2005)
      #   gstAllMax <- ((2 - 1)*(1 - hsMean)) / ((2 - 1) + hsMean)
      #   # overall Hedricks' Gst
      #   gstAllHedrick <- gstAll/gstAllMax
      #   # Overall D_jost (Jost 2008)
      #   djostAll <- (dstAll/(1-hsMean))*(2/(2-1))
      
      # Calculate estimated stats
      
      # Overall estimated dst
      dstEstAll <- htEstMean - hsEstMean
      # Overall estimated Gst (Nei & Chesser, 1983)
      gstEstAll <- dstEstAll/htEstMean
      # Overall estimated max Gst (Hedricks 2005)
      gstEstAllMax <- ((2-1)*(1-hsEstMean))/(2-1+hsEstMean)
      # Overall estimated Hedricks' Gst
      gstEstAllHed <- gstEstAll/gstEstAllMax
      # Overall estimated D_Jost (Chao et al., 2008)
      if(nloci == 1){
        djostEstAll <- (2/(2-1))*((dstEstAll)/(1 - hsEstMean))
      } else {
        dLocEstMn <- Reduce(`+`, dLocEst)/nloci
        # calculate variance (convert dLocEst to an array)
        dLocEst1 <- array(unlist(dLocEst), 
                          dim = c(nrow(dLocEst[[1]]), 
                                  ncol(dLocEst[[1]]), 
                                  length(dLocEst)))
        dLocEstVar <- apply(dLocEst1, c(1,2), var)
        djostEstAll <- 1/((1/dLocEstMn)+((dLocEstVar*((1/dLocEstMn)^3))))
        # tidy up
        rm(dLocEstMn, dLocEstVar)
        z <- gc()
        rm(z)
      }
      
      # define a function to arrange locus stats into arrays
      #   arrDef <- function(x){
      #     return(array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
      #   }
      if(fst){
        resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll, pwTheta),
                        dim = c(nrow(gstEstAll),
                                ncol(gstEstAll),
                                4))
        lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 4))
        lstats[,,,1] <- unlist(gstLocEst)
        lstats[,,,2] <- unlist(gstHedLocEst)
        lstats[,,,3] <- unlist(dLocEst)
        lstats[,,,4] <- unlist(locTheta)
        #lstats <- array(list(gstLocEst, gstHedLocEst, dLocEst, locTheta)
        #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst, locTheta,
        #                 SIMPLIFY=FALSE)
      } else {
        resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll),
                        dim = c(nrow(gstEstAll),
                                ncol(gstEstAll), 3))
        lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 3))
        lstats[,,,1] <- unlist(gstLocEst)
        lstats[,,,2] <- unlist(gstHedLocEst)
        lstats[,,,3] <- unlist(dLocEst)
        #lstats <- list(gstLocEst, gstHedLocEst, dLocEst)
        #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst,
        #                 SIMPLIFY = FALSE)
      }
      # arrange loci into arrays
      #locOut <- lapply(lstats, arrDef)
      
      
      list(resArr = resArr,
           locOut = lstats)
    }
    ############################################################################
    # END - pwDivCalc
    ############################################################################
    # Calculate Weir & Cockerham's F-statistics (optimised)
    ############################################################################
    # pwFstWC: a function co calculate weir and cockerhams fis, fit, and fst
    ############################################################################
    pwFstWC <- function(rdat){
      #   rdat <- diveRsity::readGenepop("KK_test1v2.gen")
      pw <- combn(rdat$npops, 2)
      #   # account for loci with missing info for pops
      #   pwBadData <- function(indtyp, pw){
      #     out <- sapply(1:ncol(pw), function(i){
      #       is.element(0, indtyp[pw[,i]])
      #     })
      #   }
      #   badDat <- sapply(rdat$indtyp, pwBadData, pw = pw)
      #   if(any(badDat)){
      #     bd <- TRUE
      #   }
      #   # determine the number of loci per pw comparison
      #   nlocPw <- apply(badDat, 1, function(x){
      #     if(sum(x) > 0){
      #       nl <- rdat$nloci - sum(x)
      #     } else {
      #       nl <- rdat$nloci
      #     }
      #   })
      #   # define all good data
      #   gdDat <- lapply(1:nrow(badDat), function(i){
      #     which(!badDat[i,])
      #   })
      #   badDat <- lapply(1:nrow(badDat), function(i){
      #     which(badDat[i,])
      #   })
      # get all genotypes for each pw comparison
      allGenot <- apply(pw, 2, function(x){
        list(rdat$pop_list[[x[1]]], 
             rdat$pop_list[[x[2]]])
      })
      #   # filter bad data
      #   if(any(nlocPw != rdat$nloci)){
      #     idx <- which(nlocPw != rdat$nloci)
      #     for(i in idx){
      #       allGenot[[i]][[1]] <- allGenot[[i]][[1]][, gdDat[[i]]]
      #       allGenot[[i]][[2]] <- allGenot[[i]][[2]][, gdDat[[i]]]
      #     }
      #   }
      # unlist pw genotype data
      allGenot <- lapply(allGenot, function(x){
        return(do.call("rbind", x))
      })
      # identify unique genotypes
      #   genot <- lapply(allGenot, function(x){
      #     return(apply(x, 2, function(y){
      #       unique(na.omit(y))
      #     }))
      #   })
      # count number of genotypes per pw per loc
      
      genoCount <- lapply(allGenot, function(x){
        if(NCOL(x) == 1){
          return(list(table(x)))
        } else {
          lapply(1:ncol(x), function(i) table(x[,i]))
        }
      })
      
      
      #   genoCount <- lapply(allGenot, function(x){
      #     lapply(split(x,seq(NCOL(x))),table) # accounts for single loci
      #     #apply(x, 2, table)
      #   })
      
      
      # function to count heterozygotes
      htCount <- function(x){
        nms <- names(x)
        ncharGeno <- nchar(nms[1])
        alls <- cbind(substr(nms, 1, (ncharGeno/2)),
                      substr(nms, ((ncharGeno/2) + 1), ncharGeno))
        unqAlls <- unique(as.vector(alls))
        hetCounts <- sapply(unqAlls, function(a){
          idx <- which(rowSums(alls == a) == 1)
          return(sum(x[idx]))
        })
        # needed for no het count
        if(length(hetCounts) == 0L){
          return(NA)
        } else {
          return(hetCounts) 
        }
      }
      # hSum is the total observed hets per allele
      hSum <- lapply(genoCount, function(x){
        out <- lapply(x, htCount)
      })
      
      #   if(bd){
      #     # insert na for missing loci
      #     hSum <- lapply(seq_along(badDat), function(i){
      #       naPos <- badDat[[i]]
      #       idx <- c(seq_along(hSum[[i]]), (naPos - 0.5))
      #       return(c(hSum[[i]], rep(NA, length(naPos)))[order(idx)])
      #     }) 
      #   }
      # convert to locus orientated hSum
      hSum <- lapply(seq_along(hSum[[1]]), function(i){
        lapply(hSum, "[[", i)
      })
      
      # total ind typed per loc per pw
      indTypTot <- lapply(rdat$indtyp, function(x){
        return(apply(pw, 2, function(y){
          sum(x[y], na.rm = TRUE)
        }))
      })
      # nBar is the mean number of inds per pop
      nBar <- lapply(indTypTot, `/`, 2)
      
      # hbar per pw per loc
      hBar <- lapply(seq_along(hSum), function(i){
        divd <- indTypTot[[i]]
        return(mapply(`/`, hSum[[i]], divd, SIMPLIFY = FALSE))
      })
      
      # p per loc per pw
      pCalc <- function(x, y, pw){
        out <- lapply(seq_along(pw[1,]), function(i){
          return(cbind((x[,pw[1,i]]*(2*y[pw[1,i]])),
                       (x[,pw[2,i]]*(2*y[pw[2,i]]))))
        })
        return(out)
      }
      p <- mapply(FUN = pCalc, x = rdat$allele_freq, 
                  y = rdat$indtyp, 
                  MoreArgs = list(pw = pw), 
                  SIMPLIFY = FALSE)
      
      #   # convert p elements into array structure
      #   pArr <- lapply(p, function(x){
      #     d3 <- length(x)
      #     d2 <- 2
      #     d1 <- nrow(x[[1]])
      #     return(array(unlist(x), dim = c(d1, d2, d3)))
      #   })
      
      fstatCal <- function(indT, indtyp, hBar, nBar, p, pw, npops){
        #                 indT=indTypTot[[17]]
        #                 indtyp=rdat$indtyp[[17]]
        #                 hBar <- hBar[[17]]
        #                 nBar <- nBar[[17]]
        #                 p <- p[[17]]
        #                 pw <- pw
        #                 npops <- rdat$npops
        indLocPwSqSum <- sapply(seq_along(pw[1,]), function(i){
          return(sum(indtyp[pw[,i]]^2))
        })
        indtypPw <- lapply(1:ncol(pw), function(idx){
          return(indtyp[pw[,idx]])
        })
        nC <- indT - (indLocPwSqSum/indT)
        ptildCalc <- function(x,y){ 
          return(cbind((x[,1]/(2*y[1])),
                       (x[,2]/(2*y[2]))))
        }
        pTild <- mapply(FUN = ptildCalc, x = p, y = indtypPw,
                        SIMPLIFY = FALSE)
        pBar <- lapply(seq_along(p), function(i){
          return(rowSums((p[[i]])/(2*indT[i])))
        })
        s2 <- lapply(seq_along(pBar), function(i){
          pp <- (pTild[[i]]-pBar[[i]])^2
          pp <- cbind((pp[,1]*indtypPw[[i]][1]),
                      (pp[,2]*indtypPw[[i]][2]))
          pp <- rowSums(pp)
          return((pp/(1*nBar[i])))
        })
        A <- lapply(seq_along(pBar), function(i){
          return(pBar[[i]]*(1-pBar[[i]])-(1)*s2[[i]]/2)
        })
        # fix hBar for unequal lengths
        idx <- lapply(seq_along(A), function(i){
          out <- match(names(A[[i]]), names(hBar[[i]]))
          return(which(!is.na(out)))
        })
        A <- lapply(seq_along(A), function(i){
          return(A[[i]][idx[[i]]])
        })
        s2 <- lapply(seq_along(s2), function(i){
          return(s2[[i]][idx[[i]]])
        })
        a <- lapply(seq_along(s2), function(i){
          return(nBar[[i]]*(s2[[i]]-(A[[i]]-(hBar[[i]]/4))/(nBar[[i]]-1))/nC[[i]])
        })
        #     a <- lapply(seq_along(s2), function(i){
        #       return(a[[i]][idx[[i]]])
        #     })
        b <- lapply(seq_along(A), function(i){
          return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-((2*nBar[[i]]-1)/(4*nBar[[i]]))*hBar[[i]]))
          #return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-(2*(nBar[[i]]-1))*hBar[[i]]/(4*nBar[[i]])))
        })
        #     b <- lapply(seq_along(A), function(i){
        #       return(b[[i]][idx[[i]]])
        #     })
        cdat <- lapply(seq_along(A), function(i){
          return(hBar[[i]]/2)
        })
        #     cdat <- lapply(seq_along(A), function(i){
        #       return(cdat[[i]][idx[[i]]])
        #     })
        A <- sapply(A, sum)
        a <- sapply(a, sum)
        b <- sapply(b, sum)
        cdat <- sapply(cdat, sum)
        cdat[is.na(cdat)] <- 0
        theta <- a/(a+b+cdat)
        theta[is.nan(theta)] <- NA
        pwMat <- matrix(ncol = npops, nrow = npops)
        aMat <- matrix(ncol = npops, nrow = npops)
        bMat <- matrix(ncol = npops, nrow = npops)
        cMat <- matrix(ncol = npops, nrow = npops)
        for(i in 1:ncol(pw)){
          pwMat[pw[2,i], pw[1,i]] <- theta[i]
          aMat[pw[2,i], pw[1,i]] <- a[i]
          bMat[pw[2,i], pw[1,i]] <- b[i]
          cMat[pw[2,i], pw[1,i]] <- cdat[i]
        }
        #pwMat[is.nan(pwMat)] <- 0
        aMat[is.nan(aMat)] <- 0
        cMat[is.nan(bMat)] <- 0
        bMat[is.nan(bMat)] <- 0
        
        list(pwMat, aMat, bMat, cMat)
      }
      
      # run fstatCal for each locus
      pwLoc <- mapply(FUN = fstatCal, indT = indTypTot,
                      indtyp = rdat$indtyp, hBar = hBar,
                      nBar = nBar, p = p, 
                      MoreArgs = list(pw = pw, npops = rdat$npops),
                      SIMPLIFY = FALSE)
      return(pwLoc)
    }
    ############################################################################
    # END - pwDivCalc
    ############################################################################
    # pwBasicCalc: a small function for calculating pairwise ht and hs 
    ############################################################################
    pwBasicCalc <- function(af, sHarm, pw, npops){
      ht <- matrix(ncol = npops, nrow = npops)
      hs <- matrix(ncol = npops, nrow = npops)
      htEst <- matrix(ncol = npops, nrow = npops)
      hsEst <- matrix(ncol = npops, nrow = npops)
      for(i in 1:ncol(pw)){
        id1 <- pw[1,i]
        id2 <- pw[2,i]
        # locus ht
        ht[id2, id1] <- 1 - sum(((af[,id1] + af[,id2])/2)^2)
        # locus hs
        hs[id2, id1] <- 1 - sum((af[,id1]^2 + af[,id2]^2)/2)
        # locus hs_est
        hsEst[id2, id1] <- hs[id2, id1]*((2*sHarm[id2,id1])/(2*sHarm[id2,id1]-1))
        # locus ht_est
        htEst[id2, id1] <- ht[id2, id1] + (hsEst[id2, id1]/(4*sHarm[id2, id1]))
      }
      #   ht[is.nan(ht)] <- 0
      #   hs[is.nan(hs)] <- 0
      htEst[is.nan(htEst)] <- 0
      hsEst[is.nan(hsEst)] <- 0
      list(hsEst = hsEst,
           htEst = htEst)
    }
    ############################################################################
    # END - pwBasicCalc
    ############################################################################
    
    # define locus stat calculators
    gstCalc <- function(ht, hs){
      return((ht - hs)/ht)
    }
    
    gstHedCalc <- function(ht, hs){
      gstMax <- ((2-1)*(1-hs))/(2-1+hs)
      return(((ht-hs)/ht)/gstMax)
    }
    
    djostCalc <- function(ht, hs){
      return((2/1)*((ht-hs)/(1-hs)))
    }
    
    # calculate pairwise locus harmonic mean
    pwHarmonic <- function(lss, pw){
      np <- length(lss)
      lhrm <- matrix(ncol = np, nrow = np)
      pwSS <- cbind(lss[pw[1,]], lss[pw[2,]])
      lhrmEle <- (0.5 * ((pwSS[,1]^-1) + (pwSS[,2]^-1)))^-1
      for(i in 1:ncol(pw)){
        idx1 <- pw[1,i]
        idx2 <- pw[2,i]
        lhrm[idx2, idx1] <- lhrmEle[i]
      }
      return(lhrm)
    }
    ############################################################################
    # pwDivCalc: a small function for calculating pairwise ht and hs 
    ############################################################################
    pwDivCalc <- function(x, pw, npops){
      ht <- matrix(ncol = npops, nrow = npops)
      hs <- matrix(ncol = npops, nrow = npops)
      for(i in 1:ncol(pw)){
        gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
        f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
        ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
        ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
        hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
        hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
      }
      ht[is.nan(ht)] <- 0
      hs[is.nan(hs)] <- 0
      list(ht = ht, 
           hs = hs)
    }
    ############################################################################
    # END - pwDivCalc
    ############################################################################
    
    ############################################################################
    ############################################################################
    # working well 24/10/13
    if(pWise || bspw){
      # get pw names
      pw <- combn(accDat$npops, 2)
      popNms <- accDat$pop_names
      # for pw bootstrap table
      pw_nms <- paste(popNms[pw[1,]], popNms[pw[2,]], sep = " vs. ")
      
      pwStats <- pwCalc(D, fst, bs = FALSE)
      # extract stats
      gstPW <- pwStats$resArr[,,1]
      gstHPW <- pwStats$resArr[,,2]
      dPW <- pwStats$resArr[,,3]
      if(fst){
        thetaPW <- pwStats$resArr[,,4]
      }
      # clean up
      locstats <- pwStats$locOut
      rm(pwStats)
      z <- gc()
      rm(z)
      spc1 <- rep("", ncol(gstPW))
      if(fst){
        statNms <- c("Gst_est", "G'st_est", "Djost_est", "Fst_WC")
        outobj <- rbind(c(statNms[1], spc1), 
                        c("", popNms),
                        cbind(popNms, round(gstPW, 4)),
                        c(statNms[2], spc1),
                        c("", popNms),
                        cbind(popNms, round(gstHPW, 4)), 
                        c(statNms[3], spc1),
                        c("", popNms),
                        cbind(popNms, round(dPW, 4)), 
                        c(statNms[4], spc1),
                        c("", popNms),
                        cbind(popNms, round(thetaPW, 4)))
        outobj[is.na(outobj)] <- ""
        pwMatListOut <- list(gstPW, gstHPW, dPW, thetaPW)
        # add names to pwMatListOut
        names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
        # tidy up
        rm(gstPW, gstHPW, dPW, thetaPW)
        z <- gc()
        rm(z)
      } else {
        statNms <- c("Gst_est", "G'st_est", "Djost_est")
        outobj <- rbind(c(statNms[1], spc1), 
                        c("", popNms),
                        cbind(popNms, round(gstPW, 4)),
                        c(statNms[2], spc1),
                        c("", popNms),
                        cbind(popNms, round(gstHPW, 4)), 
                        c(statNms[3], spc1),
                        c("", popNms),
                        cbind(popNms, round(dPW, 4)))
        outobj[is.na(outobj)] <- ""
        pwMatListOut <- list(gstPW, gstHPW, dPW)
        # add names to pwMatListOut
        names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst")
        # tidy up
        rm(gstPW, gstHPW, dPW)
        z <- gc()
        rm(z)
      }
      if(!is.null(on)){
        if(write_res == TRUE){
          # write data to excel
          # Load dependencies
          # pw stats
          write.xlsx(outobj, file = paste(of, "[fastDivPart].xlsx", sep=""),
                     sheetName = "Pairwise-stats", col.names = FALSE,
                     row.names = FALSE, append = TRUE)
        } else {
          # text file alternatives
          pw_outer <- file(paste(of, "Pairwise-stats[fastDivPart].txt", sep=""), 
                           "w")
          for(i in 1:nrow(outobj)){
            cat(outobj[i,], "\n", file = pw_outer, sep = "\t")
          }
          close(std)
        }
      }
      for(i in 1:length(pwMatListOut)){
        dimnames(pwMatListOut[[i]]) <- list(popNms, popNms)
      }
      # convert locstats in to list format
      locstats <- lapply(apply(locstats, 4, list), function(x){
        lapply(apply(x[[1]], 3, list), function(y){
          out <- y[[1]]
          dimnames(out) <- list(popNms, popNms)
          return(out)
        })
      })
      # add names etc
      if(fst){
        # prepare locstats for output
        names(locstats) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
      } else {
        names(locstats) <- c("gstEst", "gstEstHed", "djostEst")
      }
      for(i in 1:length(locstats)){
        names(locstats[[i]]) <- accDat$locus_names
      }
      
    }
    
    #Bootstrap
    if(bspw == TRUE){
      if (para && para_pack) {
        library(parallel)
        cl <- makeCluster(detectCores())
        clusterExport(cl, c("pwCalc", "fst", "D", "readGenepopX",
                            "fileReader", "pwFstWC", "pwHarmonic",
                            "pwBasicCalc", "djostCalc", "gstCalc",
                            "gstHedCalc"), 
                      envir = environment())
        pwBsStat <- parLapply(cl, 1:bstrps, function(...){
          return(pwCalc(infile = D, fst, bs = TRUE))
        })
        stopCluster(cl)
      } else {
        pwBsStat <- lapply(1:bstrps, function(...){
          return(pwCalc(D, fst, bs = TRUE))
        })
      }
      
      
      # seperate each stat
      
      gstEst <- lapply(pwBsStat, function(x){
        x$resArr[,,1]
      })
      
      gstEstHed <- lapply(pwBsStat, function(x){
        x$resArr[,,2]
      })
      
      dEst <- lapply(pwBsStat, function(x){
        x$resArr[,,3]
      })
      
      if(fst){
        theta <- lapply(pwBsStat, function(x){
          x$resArr[,,4]
        })
      }
      pwBsLoc <- lapply(pwBsStat, "[[", 2)
      # tidy up
      rm(pwBsStat)
      z <- gc()
      rm(z)
      
      # convert bs lists to arrays for calculations
      if(fst){
        stats <- list(gstEst = array(unlist(gstEst),
                                     dim = c(nrow(gstEst[[1]]),
                                             nrow(gstEst[[1]]),
                                             bstrps)),
                      gstEstHed = array(unlist(gstEstHed),
                                        dim = c(nrow(gstEstHed[[1]]),
                                                nrow(gstEstHed[[1]]),
                                                bstrps)),
                      dEst = array(unlist(dEst),
                                   dim = c(nrow(dEst[[1]]),
                                           nrow(dEst[[1]]),
                                           bstrps)),
                      theta = array(unlist(theta),
                                    dim = c(nrow(theta[[1]]),
                                            nrow(theta[[1]]),
                                            bstrps)))
        
      } else {
        stats <- list(gstEst = array(unlist(gstEst),
                                     dim = c(nrow(gstEst[[1]]),
                                             nrow(gstEst[[1]]),
                                             bstrps)),
                      gstEstHed = array(unlist(gstEstHed),
                                        dim = c(nrow(gstEstHed[[1]]),
                                                nrow(gstEstHed[[1]]),
                                                bstrps)),
                      dEst = array(unlist(dEst),
                                   dim = c(nrow(dEst[[1]]),
                                           nrow(dEst[[1]]),
                                           bstrps)))
      }
      # tidy up
      if(fst){
        rm(dEst, gstEst, gstEstHed, theta)
        z <- gc()
        rm(z) 
      } else {
        # tidy up
        z <- gc()
        rm(z) 
      }
      # convert locus stats into arrays for CI calculations
      npops <- accDat$npops
      nloci <- accDat$nloci
      if(fst){
        locStats <- list(
          # nei's Gst
          gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,1])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Hedrick's Gst
          gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,2])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Jost's D
          dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,3])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Weir & Cockerham's Fst
          thetaLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,4])
          })), dim = c(npops, npops, nloci, bstrps))
        )
      } else {
        locStats <- list(
          # nei's Gst
          gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,1])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Hedrick's Gst
          gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,2])
          })), dim = c(npops, npops, nloci, bstrps)),
          # Jost's D
          dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
            return(x[,,,3])
          })), dim = c(npops, npops, nloci, bstrps))
        )
      }
      # convert locus bs stats into seperate loci
      locStats <- lapply(locStats, function(x){
        lapply(apply(x, 3, list), function(y){
          return(y[[1]])
        })
      })
      
      # calculate bias corrected CI
      
      biasCor <- function(param, bs_param){
        mnBS <- apply(bs_param, c(1,2), mean, na.rm = TRUE)
        mnBS[is.nan(mnBS)] <- NA
        mnBS <- mnBS - param
        bs_param <- sweep(bs_param, c(1:2), mnBS, "-")
        return(bs_param)
      }
      biasFix <- lapply(1:length(locstats), function(i){
        mapply(FUN = biasCor, param = locstats[[i]], bs_param = locStats[[i]],
               SIMPLIFY = FALSE)
      })
      # works well
      if (para && para_pack) {
        library(parallel)
        cl <- makeCluster(detectCores())
        bcLocLCI <- parLapply(cl, biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
          })
          return(loc)
        })
        bcLocUCI <- parLapply(cl, biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
          })
          return(loc)
        })
        stopCluster(cl)
      } else {
        bcLocLCI <- lapply(biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
          })
          return(loc)
        })
        bcLocUCI <- lapply(biasFix, function(x){
          loc <- lapply(x, function(y){
            apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
          })
          return(loc)
        })
      }
      # organise data into output format
      # define function
      dfSort <- function(act, Low, High, pw_nms){
        df <- data.frame(actual = as.vector(act[lower.tri(act)]),
                         lower = as.vector(Low[lower.tri(Low)]),
                         upper = as.vector(High[lower.tri(High)]))
        rownames(df) <- pw_nms
        df[is.nan(as.matrix(df))] <- NA
        return(df)
      }
      pwLocOutput <- lapply(1:length(locstats), function(i){
        mapply(FUN = dfSort, act = locstats[[i]], Low = bcLocLCI[[i]],
               High = bcLocUCI[[i]], MoreArgs = list(pw_nms = pw_nms),
               SIMPLIFY = FALSE)
      })
      if(fst){
        names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
      } else {
        names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst")
      }
      for(i in 1:length(pwLocOutput)){
        names(pwLocOutput[[i]]) <- accDat$locus_names
      }
      
      
      
      # uncomment for standard CIs
      #       # calculate the CIs per locus
      #       if (para && para_pack) {
      #         library(parallel)
      #         cl <- makeCluster(detectCores())
      #         locLCI <- parLapply(cl, locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #         # upper CIs
      #         locUCI <- parLapply(cl, locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #       } else {
      #         locLCI <- lapply(locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #         # upper CIs
      #         locUCI <- lapply(locStats, function(x){
      #           loc <- lapply(x, function(y){
      #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
      #           })
      #           return(loc)
      #         })
      #       }
      
      #         
      #       # organise locus stats into dataframe
      #       # pw sorter
      #       pwSorter <- function(x, pw){
      #         return(data.frame(actual = sapply(1:ncol(pw), function(i){
      #           x[[1]][pw[2,i], pw[1,i]]
      #           }), 
      #           Lower = sapply(1:ncol(pw), function(i){
      #             x[[2]][pw[2,i], pw[1,i]]
      #           }),
      #           Upper = sapply(1:ncol(pw), function(i){
      #             x[[3]][pw[2,i], pw[1,i]]
      #           })
      #         ))
      #       }
      #       locOutCol <- lapply(locOut, function(x){
      #         lapply(x, pwSorter, pw = pw)
      #       })
      
      
      
      
      
      
      # organise data
      # calculate bias for cis
      biasCalc <- function(param, bs_param, pw){
        #bias <- param
        for(i in 1:ncol(pw)){
          dat <- bs_param[pw[2,i], pw[1,i], ]
          t0 <- param[pw[2,i], pw[1,i]]
          mnBS <- mean(dat , na.rm = TRUE) - t0
          bs_param[pw[2,i], pw[1,i], ] <- bs_param[pw[2,i], pw[1,i], ] - mnBS
        }
        return(bs_param)
      }
      
      # try adjusting bootstrapped estimate using bias
      
      bcStats <- mapply(biasCalc, param = pwMatListOut, bs_param = stats, 
                        MoreArgs = list(pw = pw), SIMPLIFY = FALSE)
      
      # calculate the upper and lower 95% ci
      lowCI <- lapply(stats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.025, na.rm = TRUE))
      })
      
      # bias corrected
      bcLowCI <- lapply(bcStats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.025, na.rm = TRUE))
      })
      
      
      upCI <- lapply(stats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.975, na.rm = TRUE))
      })
      
      # bias corrected
      bcHighCI <- lapply(bcStats, function(x){
        return(apply(x, c(1,2), quantile, probs = 0.975, na.rm = TRUE))
      })
      
      
      statMean <- lapply(stats, function(x){
        return(apply(x, c(1,2), mean, na.rm = TRUE))
      })
      
      # bias corrected
      bcStatMean <- lapply(bcStats, function(x){
        return(apply(x, c(1,2), mean, na.rm = TRUE))
      }) 
      
      # tidy up
      rm(stats)
      z <- gc()
      rm(z)
      
      # organize ci and mean into output structure
      pw <- combn(ncol(lowCI[[1]]), 2)
      outOrg <- function(t0 ,t1 , t2, l1, l2, u1, u2, pw, pwNms){
        out <- matrix(ncol = 7, nrow = ncol(pw))
        colnames(out) <- c("actual", "mean", "BC_mean", "Lower_95%CI", 
                           "Upper_95%CI", "BC_Lower_95%CI", "BC_Upper_95%CI")
        rownames(out) <- pwNms
        for(i in 1:ncol(pw)){
          idx <- as.vector(rev(pw[,i]))
          out[i,] <- c(t0[idx[1], idx[2]], t1[idx[1], idx[2]],
                       t2[idx[1], idx[2]], l1[idx[1], idx[2]],
                       l2[idx[1], idx[2]], u1[idx[1], idx[2]],
                       u2[idx[1], idx[2]])
        }
        
        return(out)
      }
      outputStat <- mapply(FUN = outOrg, pwMatListOut, statMean,
                           bcStatMean, lowCI, upCI, bcLowCI, bcHighCI,  
                           MoreArgs = list(pw = pw, pwNms = pw_nms),
                           SIMPLIFY = FALSE)
      
      pw_res <- outputStat
      if(fst){
        names(pw_res) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
      } else {
        names(pw_res) <- c("gstEst", "gstEstHed", "djostEst")
      }
      
      # define pwWrite for output
      sprt <- lapply(names(pw_res), FUN = `c`, c("", "", "", "", "", "", ""))
      pwWrite <- lapply(pw_res, function(x){
        comparison <- rownames(x)
        cols <- colnames(x)
        rownames(x) <- NULL
        out <- cbind(comparison, round(x, 4))
        out <- rbind(colnames(out), out)
        colnames(out) <- NULL
        return(out)
      })
      pwWrite <- mapply(FUN = "rbind", sprt, pwWrite, SIMPLIFY = FALSE)
      pwWrite <- do.call("rbind", pwWrite)
      # write results
      if(!is.null(on)){
        if(write_res==TRUE){
          write.xlsx(pwWrite, file = paste(of, "[fastDivPart].xlsx", sep = ""),
                     sheetName = "Pairwise_bootstrap", col.names = FALSE,
                     row.names = FALSE, append = TRUE)
        } else {
          # text file alternatives
          pw_bts <- file(paste(of, "Pairwise-bootstrap[fastDivPart].txt", sep = ""),
                         "w")
          #cat(paste(colnames(pw_bs_out),sep=""),"\n",sep="\t",file=pw_bts)
          for(i in 1:nrow(pwWrite)){
            cat(pwWrite[i,], "\n", file = pw_bts, sep = "\t")
          }
          close(pw_bts)
        }
      } 
    }
    zzz<-gc()
    rm(zzz)
    ############################################################################
    #pw plotter
    if(plot_res==TRUE && plt==TRUE && bspw==TRUE){
      pwso <- list()
      for(i in 1:length(pw_res)){
        pwso[[i]] <- order(pw_res[[i]][, 1], decreasing = FALSE)
        #if(length(pwso[[i]]) >= 100){
        #  pwso[[i]]<-pwso[[i]][(length(pwso[[i]])-99):length(pwso[[i]])]
        #}
      }
      if(fst){
        names(pwso) <- namer[-c(1:3, length(namer))]
      } else {
        names(pwso) <- namer[-(1:3)]
      }
      
      # define plot parameters 
      plot.call_pw<-list()
      plot.extras_pw<-list()
      xy.labels_pw<-list()
      y.pos_pw<-list()
      x.pos_pw=1:length(pwso[[i]])
      fn_pre_pw<-list()
      direct=of
      #Plot Gst_Nei
      plot.call_pw[[1]]=c("plot(pw_res[[1]][pwso[[1]],1],
                          ylim=c(0,(max(pw_res[[1]][,3])+
                          min(pw_res[[1]][,3]))),xaxt='n',
                          ylab=names(pw_res)[1],type='n',
                          xlab='Pairwise comparisons 
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[1]]=c("points(pw_res[[1]][pwso[[1]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[1]]),pw_res[[1]][pwso[[1]],6],
                            1:length(pwso[[1]]),pw_res[[1]][pwso[[1]],6],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[5]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[1]] = data.frame(pairwise_name = pw_nms[pwso[[1]]],
                                     Gst_Nei = round(pw_res[[1]][pwso[[1]], 1],4),
                                     Gst_Hedrick = round(pw_res[[2]][pwso[[1]], 1],4),
                                     D_jost = round(pw_res[[3]][pwso[[1]], 1],4))
      
      y.pos_pw[[1]] = pw_res[[1]][pwso[[1]], 1]
      fn_pre_pw[[1]] <- names(pw_res)[1]
      
      
      
      # Plot Gst_Hedrick
      plot.call_pw[[2]]=c("plot(pw_res[[2]][pwso[[2]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(pw_res)[2],type='n',
                          xlab='Pairwise comparisons
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[2]]=c("points(pw_res[[2]][pwso[[2]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[2]]),pw_res[[2]][pwso[[2]],6],
                            1:length(pwso[[2]]),pw_res[[2]][pwso[[2]],7],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[6]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[2]] = data.frame(pairwise_name = pw_nms[pwso[[2]]],
                                     Gst_Nei = round(pw_res[[1]][pwso[[2]],1],4),
                                     Gst_Hedrick = round(pw_res[[2]][pwso[[2]],1],4),
                                     D_jost = round(pw_res[[3]][pwso[[2]],1],4))
      
      y.pos_pw[[2]] = pw_res[[2]][pwso[[2]],1]
      fn_pre_pw[[2]] <- names(pw_res)[2]
      
      
      # Plot D_jost
      plot.call_pw[[3]]=c("plot(pw_res[[3]][pwso[[3]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(pw_res)[3],type='n',
                          xlab='Pairwise comparisons 
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[3]]=c("points(pw_res[[3]][pwso[[3]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[3]]),pw_res[[3]][pwso[[3]],6],
                            1:length(pwso[[3]]),pw_res[[3]][pwso[[3]],7],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[7]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[3]]=data.frame(pairwise_name=pw_nms[pwso[[3]]],
                                   Gst_Nei=round(pw_res[[1]][pwso[[3]],1],4),
                                   Gst_Hedrick=round(pw_res[[2]][pwso[[3]],1],4),
                                   D_jost=round(pw_res[[3]][pwso[[3]],1],4))
      
      y.pos_pw[[3]]=pw_res[[3]][pwso[[3]],1]
      fn_pre_pw[[3]]<-names(pw_res)[3]
      #plot(Fst_WC)
      if(fst==TRUE){
        plot.call_pw[[4]]=c("plot(pw_res[[4]][pwso[[4]],1],
                            ylim=c(0,(max(pw_res[[4]][,3])+
                            min(pw_res[[4]][,3]))),xaxt='n',ylab=names(pw_res)[4],type='n',
                            xlab='Pairwise comparisons 
                            \n (Hover over a point to see pairwise info.)',
                            cex.lab=1.2,cex.axis=1.3,las=1)")
        
        plot.extras_pw[[4]]=c("points(pw_res[[4]][pwso[[4]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],6],
                              1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],7],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=as.numeric(plot_data321[7]),
                              lwd=1,lty=2,col='red')")
        
        xy.labels_pw[[4]]=data.frame(pairwise_name=pw_nms[pwso[[4]]],
                                     Gst_Nei=round(pw_res[[1]][pwso[[4]],1],4),
                                     Gst_Hedrick=round(pw_res[[2]][pwso[[4]],1],4),
                                     D_jost=round(pw_res[[3]][pwso[[4]],1],4),
                                     Fst_WC=round(pw_res[[4]][pwso[[4]],1],4))
        
        y.pos_pw[[4]]=pw_res[[4]][pwso[[4]],1]
        fn_pre_pw[[4]]<-names(pw_res)[4]
      }
    }
    ############################### Bootstrap end ##############################
    
    
    ################################# Plot resuts ###############################
    #make necessary data available
    if(plt==TRUE && plot_res==TRUE && bsls==TRUE && bspw==TRUE){
      pl<-list(bs_res=bs_res,
               pw_res=pw_res,
               accDat=accDat,
               lso123=lso123,
               pwso=pwso,
               plot.call_loci=plot.call_loci,
               plot.extras_loci=plot.extras_loci,
               xy.labels_loci=xy.labels_loci,
               x.pos_loci=x.pos_loci,
               y.pos_loci=y.pos_loci,
               fn_pre_loci=fn_pre_loci,
               direct=direct,
               plot_loci="TRUE",
               plot_pw="TRUE",
               plot.call_pw=plot.call_pw,
               plot.extras_pw=plot.extras_pw,
               xy.labels_pw=xy.labels_pw,
               y.pos_pw=y.pos_pw,
               fn_pre_pw=fn_pre_pw,
               x.pos_pw=x.pos_pw,
               pw=pw,
               plot_data321=plot_data321,
               fst=fst)
    } else if (plt==TRUE && plot_res==TRUE && bsls==TRUE && bspw==FALSE){
      pl<-list(bs_res=bs_res,
               accDat=accDat,
               lso123=lso123,
               plot.call_loci=plot.call_loci,
               plot.extras_loci=plot.extras_loci,
               xy.labels_loci=xy.labels_loci,
               x.pos_loci=x.pos_loci,
               y.pos_loci=y.pos_loci,
               fn_pre_loci=fn_pre_loci,
               direct=direct,
               plot_loci="TRUE",
               plot_pw="FALSE",
               plot_data321=plot_data321,
               fst=fst)
    } else if (plt==TRUE && plot_res==TRUE && bsls==FALSE && bspw==TRUE){
      pl<-list(pw_res=pw_res,
               accDat=accDat,
               pwso=pwso,
               plot.call_pw=plot.call_pw,
               plot.extras_pw=plot.extras_pw,
               xy.labels_pw=xy.labels_pw,
               x.pos_pw=x.pos_pw,
               y.pos_pw=y.pos_pw,
               fn_pre_pw=fn_pre_pw,
               direct=direct,
               plot_loci="FALSE",
               plot_pw="TRUE",
               pw=pw,plot_data321=plot_data321,
               fst=fst)
    }
    if(!is.null(on)){
      if (plt==TRUE && plot_res==TRUE){
        suppressWarnings(plotter(x=pl,img="1000x600"))
      }
    }
    zzz<-gc()
    rm(zzz)
    
    if(pWise | bspw){
      # Create mean pairwise values (for Erin Landguth 12/12)
      meanPairwise <- lapply(pwMatListOut, function(x){
        mean(x, na.rm = TRUE)
      })
      names(meanPairwise) <- names(pwMatListOut)
    }
    
    
    #############################################################################
    #Data for output
    if(bspw == TRUE && bsls == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_locus = bs_res1,
           bs_pairwise = pw_res,
           bs_pairwise_loci = pwLocOutput)
    } else if(bspw == TRUE && bsls == FALSE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_pairwise = pw_res,
           bs_pairwise_loci = pwLocOutput)
    } else if(bspw == FALSE && bsls == TRUE && pWise == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_locus = bs_res1,
           pw_locus = locstats)
    } else if(bspw == FALSE && bsls == FALSE && pWise == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           pw_locus = locstats)
    } else if(bspw == FALSE && bsls == TRUE && pWise == FALSE){
      list(standard = ot1out,
           estimate = ot2out,
           bs_locus = bs_res1)
    } else if(bspw == FALSE && bsls == FALSE && pWise == FALSE){
      list(standard = ot1out,
           estimate = ot2out)
    }
  }
}
################################################################################
# fastsDivPart end                                                             #
################################################################################
################################################################################
# fastsDivPart end                                                             #
################################################################################
#
#
#
#
#
#
#
################################################################################
# haploDiv: calculate Weir & Cockerham's Fst from haploid genotypes            #
################################################################################
# haploDiv function for calculating various statistics from haploid data
# try diploidization first
#' @export
haploDiv <- function(infile = NULL, outfile = NULL, pairwise = FALSE, 
                     bootstraps = 0){
  if(bootstraps != 0){
    bs_pairwise <- TRUE
    para <- TRUE
  } else {
    bs_pairwise <- FALSE
    para <- FALSE
  }
  haploFileReader <- function(x){
    fileReader <- function(infile){
      if (typeof(infile) == "list") {
        return(infile)
      } else if (typeof(infile) == "character") {
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if (ext == "arp") {
          convRes <- arp2gen(infile)
          if (!is.null(convRes)) {
            cat("Arlequin file converted to genepop format! \n")
            infile <- paste(flForm[1], ".gen", sep = "")
          } else {
            infile <- paste(flForm[1], ".gen", sep = "")
          }
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
          locs <- strsplit(dat[2], split = "\\s+")[[1]]
          if(length(locs != 1)){
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
        
        
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if (popLoc[1] == 3) {
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:length(dat)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x) {
          x <- unlist(strsplit(x, split = "\\s+"))
          if (is.element("", x)) {
            x <- x[-(which(x == ""))]
          }
          if (is.element(",", x)) {
            x <- x[-(which(x == ","))]
          }
          if (length(x) != 1 && length(x) != no_col) {
            x <- paste(x, collapse = "")
          }
          if (length(x) < no_col) {
            tabs <- paste(rep(NA, (no_col - length(x))), 
                          sep = "\t", collapse = "\t")
            line <- paste(x, tabs, sep = "\t")
            line <- unlist(strsplit(line, split = "\t"))
            return(line)
          } else {
            return(x)
          }
        })
      }
      out <- as.data.frame(t(dat1))
      rownames(out) <- NULL
      return(out)
    }
    z <- as.matrix(fileReader(x))
    nloc <- ncol(z)-1
    popStrt <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(z[,1])) + 1
    popEnd <- c(popStrt[-1] - 2, nrow(z))
    gp <- nchar(as.character(z[popStrt[1], 2]))
    diploFun <- function(strt, end, dat, gp){
      tst <- t(sapply(strt:end, function(i){
        if(gp == 2){
          paste(sprintf("%02g", as.numeric(dat[i, -1])), 
                sprintf("%02g", as.numeric(dat[i, -1])), sep = "")
        } else if(gp == 3){
          paste(sprintf("%03g", as.numeric(dat[i, -1])), 
                sprintf("%03g", as.numeric(dat[i, -1])), sep = "")
        } else {
          cat("There is a problem with your input file!")
        }
      }))
      return(tst)
    }
    
    diploidGeno <- mapply(diploFun, strt = popStrt, end = popEnd, 
                          MoreArgs = list(dat = z, gp = gp), SIMPLIFY = FALSE)
    dat <- as.matrix(z)
    for(i in 1:length(popEnd)){
      dat[popStrt[i]:popEnd[i], -1] <- diploidGeno[[i]]
    }
    list(data = as.data.frame(dat),
         gp = gp)
  }
  # read the file and diploidize
  dat <- haploFileReader(infile)
  out <- diveRsity::fastDivPart(infile = dat$data, outfile = NULL, 
                                gp = dat$gp, pairwise = pairwise, 
                                WC_Fst = TRUE, bootstraps = bootstraps,
                                bs_pairwise = bs_pairwise, parallel = para)
  if(pairwise && bootstraps > 0){
    output <- list(locus = out$estimate[-nrow(out$estimate), "Fst_WC"],
                   overall = out$estimate[nrow(out$estimate), "Fst_WC"],
                   pairwise = out$pairwise$thetaWC,
                   bs_pairwise = out$bs_pairwise$thetaWC) 
  } else if(pairwise && bootstraps == 0L){
    output <- list(locus = out$estimate[-nrow(out$estimate), "Fst_WC"],
                   overall = out$estimate[nrow(out$estimate), "Fst_WC"],
                   pairwise = out$pairwise$thetaWC) 
  } else{
    output <- list(locus = out$estimate[-nrow(out$estimate), "Fst_WC"],
                   overall = out$estimate[nrow(out$estimate), "Fst_WC"])
  }
  
  # write reuslts to file
  if(!is.null(outfile) && pairwise && bs_pairwise){
    pwmat <- round(output$pairwise, 4)
    pwmat[is.na(pwmat)] <- ""
    idx <- 1:ncol(pwmat)
    for(i in 1:length(idx)){
      pwmat[idx[i], idx[i]] <- "--"
    }
    # pairwise matrix
    fl <- file(paste(outfile, "-[pwMatrix].txt", sep = ""), "w")
    cat("Pairwise Fst (Weir & Cockerham, 1984)", "\n", sep = "\t", file = fl)
    cat("", "\n", sep = "\t", file = fl)
    cat(c("", colnames(output$pairwise)), "\n", sep = "\t", file = fl)
    for(i in 1:nrow(output$pairwise)){
      cat(c(colnames(output$pairwise)[i], pwmat[i, ]), 
          "\n", sep = "\t", file = fl)
    }
    close(fl)
    # Pairwise cis
    fl <- file(paste(outfile, "-[pwBootstrap].txt", sep = ""), "w")
    cat("Bootstrapped 95% Confidence intervals for Weir & Cockerham's (1984) Fst",
        "\n", sep = "\t", file = fl)
    cat("", "\n", file = fl)
    cat(c("","actual", "mean", "BC_mean", "lower", "upper", 
          "BC_lower", "BC_upper"), "\n", sep = "\t", file = fl)
    for(i in 1:nrow(output$bs_pairwise)){
      cat(c(rownames(output$bs_pairwise)[i], round(output$bs_pairwise[i, ], 4)), 
          "\n", sep = "\t", file = fl)
    }
    close(fl)
  } else if(!is.null(outfile) && pairwise && !bs_pairwise){
    # pairwise matrix
    pwmat <- round(output$pairwise, 4)
    pwmat[is.na(pwmat)] <- ""
    idx <- 1:ncol(pwmat)
    for(i in 1:length(idx)){
      pwmat[idx[i], idx[i]] <- "--"
    }
    fl <- file(paste(outfile, "-[pwMatrix].txt", sep = ""), "w")
    cat("Pairwise Fst (Weir & Cockerham, 1984)", "\n", sep = "\t", file = fl)
    cat("", "\n", sep = "\t", file = fl)
    cat(c("", colnames(output$pairwise)), "\n", sep = "\t", file = fl)
    for(i in 1:nrow(output$pairwise)){
      cat(c(colnames(output$pairwise)[i], pwmat[i, ]), 
          "\n", sep = "\t", file = fl)
    }
    close(fl)
  }
  return(output)
}
################################################################################
# haploDiv end                                                                 #
################################################################################
#
#
#
#
#
#
#' writeBoot: A function to write unprocessed bootstrapped matrices of 
#' differentiation and diversity partition statistics to file.
#' 
#' The function is based on fastDivPart from diveRsity.
#' 
#' Kevin Keenan 2014
#' 
#' 
writeBoot <- function(infile = NULL, outfile = NULL, gp = 3, bootstraps = 0, 
                      parallel = FALSE){
  ############################ Argument definitions ############################
  # define arguments for testing
  D <- infile
  on <- outfile
  gp <- gp
  fst <- TRUE
  bstrps <- bootstraps
  bsls <- FALSE
  bspw <- TRUE
  plt <- FALSE
  para <- parallel
  pWise <- TRUE
  
  ##############################################################################
  #Use pre.div to calculate the standard global and locus stats
  accDat <- pre.divLowMemory(list(infile = D,
                                  gp = gp,
                                  bootstrap = FALSE,
                                  locs = TRUE,
                                  fst = fst,
                                  min = FALSE))
  # create a directory for output
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[writeBoot]","/",sep="")))
  }
  of = paste(getwd(), "/", on, "-[writeBoot]", "/", sep = "")
  
  
  if(para){
    para_pack_inst<-is.element(c("parallel","doParallel","foreach","iterators"),
                               installed.packages()[,1])
  }
  
  para_pack <- all(para_pack_inst)
  ############################################################################
  ################################## Pairwise ################################
  ############################################################################
  # population pair combinations
  
  # define new functions
  ############################################################################
  ############################################################################
  # pwCalc
  ############################################################################
  # New optimised function for the calculation of pairwise statistics
  # Returns a 3D array where each 'slot' represents the pairwise matrix
  # for Gst_est, G_st_est_hed and D_jost_est respectively
  
  # Kevin Keenan
  # 2013
  
  pwCalc <- function(infile, fst,  bs = FALSE){
    
    
    #   # uncomment for testing
    #   infile <- "pw_test.txt"
    #   source("readGenepopX.R")
    #   # read pwBasicCalc function
    #   source("pwBasicCalc.R")
    # define baseline info
    dat <- readGenepopX(list(infile = infile,
                             bootstrap = bs))
    if(fst){
      # calculate all fst
      fstat <- pwFstWC(dat)
      # extract locus theta and variance components
      locTheta <- lapply(fstat, "[[", 1)
      # sum res
      aLoc <- Reduce(`+`, lapply(fstat, "[[", 2))
      bLoc <- Reduce(`+`, lapply(fstat, "[[", 3))
      cLoc <- Reduce(`+`, lapply(fstat, "[[", 4))
      # calculate pw Fst across loci
      pwTheta <- aLoc/(aLoc+bLoc+cLoc)
      # clean up
      rm(aLoc, bLoc, cLoc, fstat)
      z <- gc()
      rm(z)
    }
    # extract allele frequencies
    af <- dat$allele_freq
    # extract harmonic mean sample sizes
    
    # make space in RAM
    dat$allele_freq <- NULL
    z <- gc()
    rm(z)
    # extract npops and nloci
    npops <- dat$npops
    nloci <- dat$nloci
    # define pairwise relationships
    pw <- combn(dat$npops, 2)
    # generate pairwise locus harmonic mean sample sizes
    indtyp <- dat$indtyp
    pwHarm <- lapply(indtyp, pwHarmonic, pw = pw)
    
    
    # calculate pairwise ht and hs
    hths <- mapply(pwBasicCalc, af, pwHarm,
                   MoreArgs = list(pw = pw, npops = dat$npops),
                   SIMPLIFY = FALSE)
    # seperate ht and hs
    #   htLoc <- lapply(hths, "[[", 1)
    #   hsLoc <- lapply(hths, "[[", 2)
    # seperate ht_est and hs_est
    hsEstLoc <- lapply(hths, "[[", 1)
    htEstLoc <- lapply(hths, "[[", 2)
    
    # clean up
    rm(hths)
    z <- gc()
    rm(z)
    
    # Calculate locus stats
    # Standard locus stats
    # locus Gst
    #   gstLoc <- mapply(FUN = gstCalc, ht = htLoc, hs = hsLoc, 
    #                    SIMPLIFY = FALSE)
    #   # locus G'st
    #   gstHedLoc <- mapply(FUN = gstHedCalc, ht = htLoc, hs = hsLoc,
    #                       SIMPLIFY = FALSE)
    #   # locus D_jost
    #   dLoc <- mapply(FUN = djostCalc, ht = htLoc, hs = hsLoc,
    #                  SIMPLIFY = FALSE)
    
    # Estimated locus stats
    # locus Gst_est
    gstLocEst <- mapply(FUN = gstCalc, ht = htEstLoc, 
                        hs = hsEstLoc, 
                        SIMPLIFY = FALSE)
    # locus G'st_est
    gstHedLocEst <- mapply(FUN = gstHedCalc, ht = htEstLoc, 
                           hs = hsEstLoc,
                           SIMPLIFY = FALSE)
    # locus D_jost_est
    dLocEst <- mapply(FUN = djostCalc, ht = htEstLoc, 
                      hs = hsEstLoc,
                      SIMPLIFY = FALSE)
    
    #   # calculate mean ht and hs
    #   htMean <- Reduce(`+`, htLoc)/nloci
    #   hsMean <- Reduce(`+`, hsLoc)/nloci
    # calculate mean ht_est and hs_est
    htEstMean <- Reduce(`+`, htEstLoc)/nloci
    hsEstMean <- Reduce(`+`, hsEstLoc)/nloci
    
    # calculate standard stats (uncomment for loc stats)
    
    #   # overall dst
    #   dstAll <- htMean - hsMean
    #   # overall gst (Nei 1973)
    #   gstAll <- (dstAll)/htMean
    #   # overall max gst (Hedricks 2005)
    #   gstAllMax <- ((2 - 1)*(1 - hsMean)) / ((2 - 1) + hsMean)
    #   # overall Hedricks' Gst
    #   gstAllHedrick <- gstAll/gstAllMax
    #   # Overall D_jost (Jost 2008)
    #   djostAll <- (dstAll/(1-hsMean))*(2/(2-1))
    
    # Calculate estimated stats
    
    # Overall estimated dst
    dstEstAll <- htEstMean - hsEstMean
    # Overall estimated Gst (Nei & Chesser, 1983)
    gstEstAll <- dstEstAll/htEstMean
    # Overall estimated max Gst (Hedricks 2005)
    gstEstAllMax <- ((2-1)*(1-hsEstMean))/(2-1+hsEstMean)
    # Overall estimated Hedricks' Gst
    gstEstAllHed <- gstEstAll/gstEstAllMax
    # Overall estimated D_Jost (Chao et al., 2008)
    if(nloci == 1){
      djostEstAll <- (2/(2-1))*((dstEstAll)/(1 - hsEstMean))
    } else {
      dLocEstMn <- Reduce(`+`, dLocEst)/nloci
      # calculate variance (convert dLocEst to an array)
      dLocEst1 <- array(unlist(dLocEst), 
                        dim = c(nrow(dLocEst[[1]]), 
                                ncol(dLocEst[[1]]), 
                                length(dLocEst)))
      dLocEstVar <- apply(dLocEst1, c(1,2), var)
      djostEstAll <- 1/((1/dLocEstMn)+((dLocEstVar*((1/dLocEstMn)^3))))
      # tidy up
      rm(dLocEstMn, dLocEstVar)
      z <- gc()
      rm(z)
    }
    
    # define a function to arrange locus stats into arrays
    #   arrDef <- function(x){
    #     return(array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
    #   }
    if(fst){
      resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll, pwTheta),
                      dim = c(nrow(gstEstAll),
                              ncol(gstEstAll),
                              4))
      lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 4))
      lstats[,,,1] <- unlist(gstLocEst)
      lstats[,,,2] <- unlist(gstHedLocEst)
      lstats[,,,3] <- unlist(dLocEst)
      lstats[,,,4] <- unlist(locTheta)
      #lstats <- array(list(gstLocEst, gstHedLocEst, dLocEst, locTheta)
      #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst, locTheta,
      #                 SIMPLIFY=FALSE)
    } else {
      resArr <- array(c(gstEstAll, gstEstAllHed, djostEstAll),
                      dim = c(nrow(gstEstAll),
                              ncol(gstEstAll), 3))
      lstats <- array(NA, dim = c(dat$npops, dat$npops, dat$nloci, 3))
      lstats[,,,1] <- unlist(gstLocEst)
      lstats[,,,2] <- unlist(gstHedLocEst)
      lstats[,,,3] <- unlist(dLocEst)
      #lstats <- list(gstLocEst, gstHedLocEst, dLocEst)
      #lstats <- mapply(FUN = `list`, gstLocEst, gstHedLocEst, dLocEst,
      #                 SIMPLIFY = FALSE)
    }
    # arrange loci into arrays
    #locOut <- lapply(lstats, arrDef)
    
    
    list(resArr = resArr,
         locOut = lstats)
  }
  ############################################################################
  # END - pwDivCalc
  ############################################################################
  # Calculate Weir & Cockerham's F-statistics (optimised)
  ############################################################################
  # pwFstWC: a function co calculate weir and cockerhams fis, fit, and fst
  ############################################################################
  pwFstWC<-function(rdat){
    #   rdat <- diveRsity::readGenepop("KK_test1v2.gen")
    pw <- combn(rdat$npops, 2)
    #   # account for loci with missing info for pops
    #   pwBadData <- function(indtyp, pw){
    #     out <- sapply(1:ncol(pw), function(i){
    #       is.element(0, indtyp[pw[,i]])
    #     })
    #   }
    #   badDat <- sapply(rdat$indtyp, pwBadData, pw = pw)
    #   if(any(badDat)){
    #     bd <- TRUE
    #   }
    #   # determine the number of loci per pw comparison
    #   nlocPw <- apply(badDat, 1, function(x){
    #     if(sum(x) > 0){
    #       nl <- rdat$nloci - sum(x)
    #     } else {
    #       nl <- rdat$nloci
    #     }
    #   })
    #   # define all good data
    #   gdDat <- lapply(1:nrow(badDat), function(i){
    #     which(!badDat[i,])
    #   })
    #   badDat <- lapply(1:nrow(badDat), function(i){
    #     which(badDat[i,])
    #   })
    # get all genotypes for each pw comparison
    allGenot <- apply(pw, 2, function(x){
      list(rdat$pop_list[[x[1]]], 
           rdat$pop_list[[x[2]]])
    })
    #   # filter bad data
    #   if(any(nlocPw != rdat$nloci)){
    #     idx <- which(nlocPw != rdat$nloci)
    #     for(i in idx){
    #       allGenot[[i]][[1]] <- allGenot[[i]][[1]][, gdDat[[i]]]
    #       allGenot[[i]][[2]] <- allGenot[[i]][[2]][, gdDat[[i]]]
    #     }
    #   }
    # unlist pw genotype data
    allGenot <- lapply(allGenot, function(x){
      return(do.call("rbind", x))
    })
    # identify unique genotypes
    #   genot <- lapply(allGenot, function(x){
    #     return(apply(x, 2, function(y){
    #       unique(na.omit(y))
    #     }))
    #   })
    # count number of genotypes per pw per loc
    
    genoCount <- lapply(allGenot, function(x){
      if(NCOL(x) == 1){
        return(list(table(x)))
      } else {
        lapply(1:ncol(x), function(i) table(x[,i]))
      }
    })
    
    
    #   genoCount <- lapply(allGenot, function(x){
    #     lapply(split(x,seq(NCOL(x))),table) # accounts for single loci
    #     #apply(x, 2, table)
    #   })
    
    
    # function to count heterozygotes
    htCount <- function(x){
      nms <- names(x)
      ncharGeno <- nchar(nms[1])
      alls <- cbind(substr(nms, 1, (ncharGeno/2)),
                    substr(nms, ((ncharGeno/2) + 1), ncharGeno))
      unqAlls <- unique(as.vector(alls))
      hetCounts <- sapply(unqAlls, function(a){
        idx <- which(rowSums(alls == a) == 1)
        return(sum(x[idx]))
      })
      return(hetCounts)
    }
    # hSum is the total observed hets per allele
    hSum <- lapply(genoCount, function(x){
      out <- lapply(x, htCount)
    })
    
    #   if(bd){
    #     # insert na for missing loci
    #     hSum <- lapply(seq_along(badDat), function(i){
    #       naPos <- badDat[[i]]
    #       idx <- c(seq_along(hSum[[i]]), (naPos - 0.5))
    #       return(c(hSum[[i]], rep(NA, length(naPos)))[order(idx)])
    #     }) 
    #   }
    # convert to locus orientated hSum
    hSum <- lapply(seq_along(hSum[[1]]), function(i){
      lapply(hSum, "[[", i)
    })
    
    # total ind typed per loc per pw
    indTypTot <- lapply(rdat$indtyp, function(x){
      return(apply(pw, 2, function(y){
        sum(x[y])
      }))
    })
    # nBar is the mean number of inds per pop
    nBar <- lapply(indTypTot, `/`, 2)
    
    # hbar per pw per loc
    hBar <- lapply(seq_along(hSum), function(i){
      divd <- indTypTot[[i]]
      return(mapply(`/`, hSum[[i]], divd, SIMPLIFY = FALSE))
    })
    
    # p per loc per pw
    pCalc <- function(x, y, pw){
      out <- lapply(seq_along(pw[1,]), function(i){
        return(cbind((x[,pw[1,i]]*(2*y[pw[1,i]])),
                     (x[,pw[2,i]]*(2*y[pw[2,i]]))))
      })
      return(out)
    }
    p <- mapply(FUN = pCalc, x = rdat$allele_freq, 
                y = rdat$indtyp, 
                MoreArgs = list(pw = pw), 
                SIMPLIFY = FALSE)
    
    #   # convert p elements into array structure
    #   pArr <- lapply(p, function(x){
    #     d3 <- length(x)
    #     d2 <- 2
    #     d1 <- nrow(x[[1]])
    #     return(array(unlist(x), dim = c(d1, d2, d3)))
    #   })
    
    fstatCal <- function(indT, indtyp, hBar, nBar, p, pw, npops){
      #         indT=indTypTot[[1]]
      #         indtyp=rdat$indtyp[[1]]
      #         hBar <- hBar[[1]]
      #         nBar <- nBar[[1]]
      #         p <- p[[1]]
      #         pw <- pw
      #         npops <- rdat$npops
      indLocPwSqSum <- sapply(seq_along(pw[1,]), function(i){
        return(sum(indtyp[pw[,i]]^2))
      })
      indtypPw <- lapply(1:ncol(pw), function(idx){
        return(indtyp[pw[,idx]])
      })
      nC <- indT - (indLocPwSqSum/indT)
      ptildCalc <- function(x,y){ 
        return(cbind((x[,1]/(2*y[1])),
                     (x[,2]/(2*y[2]))))
      }
      pTild <- mapply(FUN = ptildCalc, x = p, y = indtypPw,
                      SIMPLIFY = FALSE)
      pBar <- lapply(seq_along(p), function(i){
        return(rowSums((p[[i]])/(2*indT[i])))
      })
      s2 <- lapply(seq_along(pBar), function(i){
        pp <- (pTild[[i]]-pBar[[i]])^2
        pp <- cbind((pp[,1]*indtypPw[[i]][1]),
                    (pp[,2]*indtypPw[[i]][2]))
        pp <- rowSums(pp)
        return((pp/(1*nBar[i])))
      })
      A <- lapply(seq_along(pBar), function(i){
        return(pBar[[i]]*(1-pBar[[i]])-(1)*s2[[i]]/2)
      })
      # fix hBar for unequal lengths
      idx <- lapply(seq_along(A), function(i){
        out <- match(names(A[[i]]), names(hBar[[i]]))
        return(which(!is.na(out)))
      })
      A <- lapply(seq_along(A), function(i){
        return(A[[i]][idx[[i]]])
      })
      s2 <- lapply(seq_along(s2), function(i){
        return(s2[[i]][idx[[i]]])
      })
      a <- lapply(seq_along(s2), function(i){
        return(nBar[[i]]*(s2[[i]]-(A[[i]]-(hBar[[i]]/4))/(nBar[[i]]-1))/nC[[i]])
      })
      #     a <- lapply(seq_along(s2), function(i){
      #       return(a[[i]][idx[[i]]])
      #     })
      b <- lapply(seq_along(A), function(i){
        return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-((2*nBar[[i]]-1)/(4*nBar[[i]]))*hBar[[i]]))
        #return((nBar[[i]]/(nBar[[i]]-1))*(A[[i]]-(2*(nBar[[i]]-1))*hBar[[i]]/(4*nBar[[i]])))
      })
      #     b <- lapply(seq_along(A), function(i){
      #       return(b[[i]][idx[[i]]])
      #     })
      cdat <- lapply(seq_along(A), function(i){
        return(hBar[[i]]/2)
      })
      #     cdat <- lapply(seq_along(A), function(i){
      #       return(cdat[[i]][idx[[i]]])
      #     })
      A <- sapply(A, sum)
      a <- sapply(a, sum)
      b <- sapply(b, sum)
      cdat <- sapply(cdat, sum)
      theta <- a/(a+b+cdat)
      pwMat <- matrix(ncol = npops, nrow = npops)
      aMat <- matrix(ncol = npops, nrow = npops)
      bMat <- matrix(ncol = npops, nrow = npops)
      cMat <- matrix(ncol = npops, nrow = npops)
      for(i in 1:ncol(pw)){
        pwMat[pw[2,i], pw[1,i]] <- theta[i]
        aMat[pw[2,i], pw[1,i]] <- a[i]
        bMat[pw[2,i], pw[1,i]] <- b[i]
        cMat[pw[2,i], pw[1,i]] <- cdat[i]
      }
      pwMat[is.nan(pwMat)] <- NA
      aMat[is.nan(aMat)] <- NA
      cMat[is.nan(bMat)] <- NA
      bMat[is.nan(bMat)] <- NA
      
      list(pwMat, aMat, bMat, cMat)
    }
    
    # run fstatCal for each locus
    pwLoc <- mapply(FUN = fstatCal, indT = indTypTot,
                    indtyp = rdat$indtyp, hBar = hBar,
                    nBar = nBar, p = p, 
                    MoreArgs = list(pw = pw, npops = rdat$npops),
                    SIMPLIFY = FALSE)
    return(pwLoc)
  }
  ############################################################################
  # END - pwDivCalc
  ############################################################################
  # pwBasicCalc: a small function for calculating pairwise ht and hs 
  ############################################################################
  pwBasicCalc <- function(af, sHarm, pw, npops){
    ht <- matrix(ncol = npops, nrow = npops)
    hs <- matrix(ncol = npops, nrow = npops)
    htEst <- matrix(ncol = npops, nrow = npops)
    hsEst <- matrix(ncol = npops, nrow = npops)
    for(i in 1:ncol(pw)){
      id1 <- pw[1,i]
      id2 <- pw[2,i]
      # locus ht
      ht[id2, id1] <- 1 - sum(((af[,id1] + af[,id2])/2)^2)
      # locus hs
      hs[id2, id1] <- 1 - sum((af[,id1]^2 + af[,id2]^2)/2)
      # locus hs_est
      hsEst[id2, id1] <- hs[id2, id1]*((2*sHarm[id2,id1])/(2*sHarm[id2,id1]-1))
      # locus ht_est
      htEst[id2, id1] <- ht[id2, id1] + (hsEst[id2, id1]/(4*sHarm[id2, id1]))
    }
    #   ht[is.nan(ht)] <- 0
    #   hs[is.nan(hs)] <- 0
    htEst[is.nan(htEst)] <- 0
    hsEst[is.nan(hsEst)] <- 0
    list(hsEst = hsEst,
         htEst = htEst)
  }
  ############################################################################
  # END - pwBasicCalc
  ############################################################################
  
  # define locus stat calculators
  gstCalc <- function(ht, hs){
    return((ht - hs)/ht)
  }
  
  gstHedCalc <- function(ht, hs){
    gstMax <- ((2-1)*(1-hs))/(2-1+hs)
    return(((ht-hs)/ht)/gstMax)
  }
  
  djostCalc <- function(ht, hs){
    return((2/1)*((ht-hs)/(1-hs)))
  }
  
  # calculate pairwise locus harmonic mean
  pwHarmonic <- function(lss, pw){
    np <- length(lss)
    lhrm <- matrix(ncol = np, nrow = np)
    pwSS <- cbind(lss[pw[1,]], lss[pw[2,]])
    lhrmEle <- (0.5 * ((pwSS[,1]^-1) + (pwSS[,2]^-1)))^-1
    for(i in 1:ncol(pw)){
      idx1 <- pw[1,i]
      idx2 <- pw[2,i]
      lhrm[idx2, idx1] <- lhrmEle[i]
    }
    return(lhrm)
  }
  ############################################################################
  # pwDivCalc: a small function for calculating pairwise ht and hs 
  ############################################################################
  pwDivCalc <- function(x, pw, npops){
    ht <- matrix(ncol = npops, nrow = npops)
    hs <- matrix(ncol = npops, nrow = npops)
    for(i in 1:ncol(pw)){
      gamma <- sum(sqrt(abs(x[,pw[1,i]] * x[,pw[2,i]])))^-1 
      f <- gamma * sqrt(x[,pw[1,i]] * x[,pw[2,i]])
      ht[pw[1,i],pw[2,i]] <- 1 - sum(((f + x[,pw[1,i]])/2)^2)
      ht[pw[2,i],pw[1,i]] <- 1 - sum(((f + x[,pw[2,i]])/2)^2)
      hs[pw[1,i],pw[2,i]] <- 1 - sum((f^2 + x[,pw[1,i]]^2)/2)
      hs[pw[2,i],pw[1,i]] <- 1 - sum((f^2 + x[,pw[2,i]]^2)/2)
    }
    ht[is.nan(ht)] <- 0
    hs[is.nan(hs)] <- 0
    list(ht = ht, 
         hs = hs)
  }
  ############################################################################
  # END - pwDivCalc
  ############################################################################
  
  ############################################################################
  ############################################################################
  # working well 24/10/13
  if(pWise || bspw){
    # get pw names
    pw <- combn(accDat$npops, 2)
    popNms <- accDat$pop_names
    # for pw bootstrap table
    pw_nms <- paste(popNms[pw[1,]], popNms[pw[2,]], sep = " vs. ")
    
    pwStats <- pwCalc(D, fst, bs = FALSE)
    # extract stats
    gstPW <- pwStats$resArr[,,1]
    gstHPW <- pwStats$resArr[,,2]
    dPW <- pwStats$resArr[,,3]
    if(fst){
      thetaPW <- pwStats$resArr[,,4]
    }
    # clean up
    locstats <- pwStats$locOut
    rm(pwStats)
    z <- gc()
    rm(z)
    #       spc1 <- rep("", ncol(gstPW))
    #       if(fst){
    #         statNms <- c("Gst_est", "G'st_est", "Djost_est", "Fst_WC")
    #         outobj <- rbind(c(statNms[1], spc1), 
    #                         c("", popNms),
    #                         cbind(popNms, round(gstPW, 4)),
    #                         c(statNms[2], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(gstHPW, 4)), 
    #                         c(statNms[3], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(dPW, 4)), 
    #                         c(statNms[4], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(thetaPW, 4)))
    #         outobj[is.na(outobj)] <- ""
    #         pwMatListOut <- list(gstPW, gstHPW, dPW, thetaPW)
    #         # add names to pwMatListOut
    #         names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    #         # tidy up
    #         rm(gstPW, gstHPW, dPW, thetaPW)
    #         z <- gc()
    #         rm(z)
    #       } else {
    #         statNms <- c("Gst_est", "G'st_est", "Djost_est")
    #         outobj <- rbind(c(statNms[1], spc1), 
    #                         c("", popNms),
    #                         cbind(popNms, round(gstPW, 4)),
    #                         c(statNms[2], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(gstHPW, 4)), 
    #                         c(statNms[3], spc1),
    #                         c("", popNms),
    #                         cbind(popNms, round(dPW, 4)))
    #         outobj[is.na(outobj)] <- ""
    #         pwMatListOut <- list(gstPW, gstHPW, dPW)
    #         # add names to pwMatListOut
    #         names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst")
    #         # tidy up
    #         rm(gstPW, gstHPW, dPW)
    #         z <- gc()
    #         rm(z)
    #       }
    # 
    #       for(i in 1:length(pwMatListOut)){
    #         dimnames(pwMatListOut[[i]]) <- list(popNms, popNms)
    #       }
    
    # convert locstats in to list format
    #       locstats <- lapply(apply(locstats, 4, list), function(x){
    #         lapply(apply(x[[1]], 3, list), function(y){
    #           out <- y[[1]]
    #           dimnames(out) <- list(popNms, popNms)
    #           return(out)
    #         })
    #       })
    # add names etc
    #       if(fst){
    #         # prepare locstats for output
    #         names(locstats) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    #       } else {
    #         names(locstats) <- c("gstEst", "gstEstHed", "djostEst")
    #       }
    #       for(i in 1:length(locstats)){
    #         names(locstats[[i]]) <- accDat$locus_names
    #       }
    
  }
  
  #Bootstrap
  if(bspw == TRUE){
    if (para && para_pack) {
      library(parallel)
      cl <- makeCluster(detectCores())
      clusterExport(cl, c("pwCalc", "fst", "D", "readGenepopX",
                          "fileReader", "pwFstWC", "pwHarmonic",
                          "pwBasicCalc", "djostCalc", "gstCalc",
                          "gstHedCalc"), 
                    envir = environment())
      pwBsStat <- parLapply(cl, 1:bstrps, function(...){
        return(pwCalc(infile = D, fst, bs = TRUE))
      })
      stopCluster(cl)
    } else {
      pwBsStat <- lapply(1:bstrps, function(...){
        return(pwCalc(D, fst, bs = TRUE))
      })
    }
    
    
    # seperate each stat
    
    gstEst <- lapply(pwBsStat, function(x){
      x$resArr[,,1]
    })
    
    gstEstHed <- lapply(pwBsStat, function(x){
      x$resArr[,,2]
    })
    
    dEst <- lapply(pwBsStat, function(x){
      x$resArr[,,3]
    })
    
    if(fst){
      theta <- lapply(pwBsStat, function(x){
        x$resArr[,,4]
      })
    }
    #pwBsLoc <- lapply(pwBsStat, "[[", 2)
    # tidy up
    rm(pwBsStat)
    z <- gc()
    rm(z)
    
    # convert bs lists to arrays for calculations
    if(fst){
      stats <- list(gstEst = array(unlist(gstEst),
                                   dim = c(nrow(gstEst[[1]]),
                                           nrow(gstEst[[1]]),
                                           bstrps)),
                    gstEstHed = array(unlist(gstEstHed),
                                      dim = c(nrow(gstEstHed[[1]]),
                                              nrow(gstEstHed[[1]]),
                                              bstrps)),
                    dEst = array(unlist(dEst),
                                 dim = c(nrow(dEst[[1]]),
                                         nrow(dEst[[1]]),
                                         bstrps)),
                    theta = array(unlist(theta),
                                  dim = c(nrow(theta[[1]]),
                                          nrow(theta[[1]]),
                                          bstrps)))
      
    } else {
      stats <- list(gstEst = array(unlist(gstEst),
                                   dim = c(nrow(gstEst[[1]]),
                                           nrow(gstEst[[1]]),
                                           bstrps)),
                    gstEstHed = array(unlist(gstEstHed),
                                      dim = c(nrow(gstEstHed[[1]]),
                                              nrow(gstEstHed[[1]]),
                                              bstrps)),
                    dEst = array(unlist(dEst),
                                 dim = c(nrow(dEst[[1]]),
                                         nrow(dEst[[1]]),
                                         bstrps)))
    }
    # tidy up
    if(fst){
      rm(dEst, gstEst, gstEstHed, theta)
      z <- gc()
      rm(z) 
    } else {
      # tidy up
      z <- gc()
      rm(z) 
    }
    #       # convert locus stats into arrays for CI calculations
    #       npops <- accDat$npops
    #       nloci <- accDat$nloci
    #       if(fst){
    #         locStats <- list(
    #           # nei's Gst
    #           gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,1])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Hedrick's Gst
    #           gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,2])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Jost's D
    #           dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,3])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Weir & Cockerham's Fst
    #           thetaLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,4])
    #           })), dim = c(npops, npops, nloci, bstrps))
    #         )
    #       } else {
    #         locStats <- list(
    #           # nei's Gst
    #           gstLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,1])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Hedrick's Gst
    #           gstHedLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,2])
    #           })), dim = c(npops, npops, nloci, bstrps)),
    #           # Jost's D
    #           dJostLocStat = array(unlist(lapply(pwBsLoc, function(x){
    #             return(x[,,,3])
    #           })), dim = c(npops, npops, nloci, bstrps))
    #         )
    #       }
    #       # convert locus bs stats into seperate loci
    #       locStats <- lapply(locStats, function(x){
    #         lapply(apply(x, 3, list), function(y){
    #           return(y[[1]])
    #         })
    #       })
    #       
    #       # calculate bias corrected CI
    #       
    #       biasCor <- function(param, bs_param){
    #         mnBS <- apply(bs_param, c(1,2), mean, na.rm = TRUE)
    #         mnBS[is.nan(mnBS)] <- NA
    #         mnBS <- mnBS - param
    #         bs_param <- sweep(bs_param, c(1:2), mnBS, "-")
    #         return(bs_param)
    #       }
    #       biasFix <- lapply(1:length(locstats), function(i){
    #         mapply(FUN = biasCor, param = locstats[[i]], bs_param = locStats[[i]],
    #                SIMPLIFY = FALSE)
    #       })
    #       # works well
    #       if (para && para_pack) {
    #         library(parallel)
    #         cl <- makeCluster(detectCores())
    #         bcLocLCI <- parLapply(cl, biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #         bcLocUCI <- parLapply(cl, biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #         stopCluster(cl)
    #       } else {
    #         bcLocLCI <- lapply(biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #         bcLocUCI <- lapply(biasFix, function(x){
    #           loc <- lapply(x, function(y){
    #             apply(y, c(1,2), quantile, probs = 0.975, na.rm = TRUE)
    #           })
    #           return(loc)
    #         })
    #       }
    #       # organise data into output format
    #       # define function
    #       dfSort <- function(act, Low, High, pw_nms){
    #         df <- data.frame(actual = as.vector(act[lower.tri(act)]),
    #                          lower = as.vector(Low[lower.tri(Low)]),
    #                          upper = as.vector(High[lower.tri(High)]))
    #         rownames(df) <- pw_nms
    #         df[is.nan(as.matrix(df))] <- NA
    #         return(df)
    #       }
    #       pwLocOutput <- lapply(1:length(locstats), function(i){
    #         mapply(FUN = dfSort, act = locstats[[i]], Low = bcLocLCI[[i]],
    #                High = bcLocUCI[[i]], MoreArgs = list(pw_nms = pw_nms),
    #                SIMPLIFY = FALSE)
    #       })
    #       if(fst){
    #         names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    #       } else {
    #         names(pwLocOutput) <- c("gstEst", "gstEstHed", "djostEst")
    #       }
    #       for(i in 1:length(pwLocOutput)){
    #         names(pwLocOutput[[i]]) <- accDat$locus_names
    #       }
    
    
    pwMatListOut <- list(gstPW, gstHPW, dPW, thetaPW)
    # add names to pwMatListOut
    names(pwMatListOut) <- c("gstEst", "gstEstHed", "djostEst", "thetaWC")
    
    
    
    
    # organise data
    # calculate bias for cis
    biasCalc <- function(param, bs_param, pw){
      #bias <- param
      for(i in 1:ncol(pw)){
        dat <- bs_param[pw[2,i], pw[1,i], ]
        t0 <- param[pw[2,i], pw[1,i]]
        mnBS <- mean(dat , na.rm = TRUE) - t0
        bs_param[pw[2,i], pw[1,i], ] <- bs_param[pw[2,i], pw[1,i], ] - mnBS
      }
      return(bs_param)
    }
    
    # try adjusting bootstrapped estimate using bias
    
    bcStats <- mapply(biasCalc, param = pwMatListOut, bs_param = stats, 
                      MoreArgs = list(pw = pw), SIMPLIFY = FALSE)
    
    ## Write results
    # we need three files for each of the four statistics calculated
    fileNms <- list()
    for(i in 1:4){
      fileNms[[i]] <- vector()
      fileNms[[i]][1] <- paste(of, names(stats)[i], 
                               "-actual.txt", sep = "")
      fileNms[[i]][2] <- paste(of, names(stats)[i], 
                               "-uncorrected.txt", sep = "")
      fileNms[[i]][3] <- paste(of, names(stats)[i], 
                               "-corrected.txt", sep = "")
    }
    # open file connection and write results
    for(i in 1:4){
      actual <- file(fileNms[[i]][1], "w")
      uncor <- file(fileNms[[i]][2], "w")
      corr <- file(fileNms[[i]][3], "w")
      # write standard statistics
      pwMatListOut[[i]][is.na(pwMatListOut[[i]])] <- ""
      for(j in 1:nrow(pwMatListOut[[i]])){
        cat(pwMatListOut[[i]][j,], "\n", file = actual, sep = "\t")
      }
      close(actual)
      # write standard bootstraps
      # generate a bootstrap list
      out1 <- lapply(apply(stats[[i]], 3, list), "[[", 1)
      out1 <- do.call("rbind", out1)
      out1[is.na(out1)] <- ""
      out1 <- out1[,-ncol(out1)]
      for(j in 1:nrow(out1)){
        cat(out1[j,], "\n", file = uncor, sep = "\t")
      }
      close(uncor)
      rm(out1)
      out2 <- lapply(apply(bcStats[[i]], 3, list), "[[", 1)
      out2 <- do.call("rbind", out2)
      out2[is.na(out2)] <- ""
      out2 <- out2[,-ncol(out2)]
      for(j in 1:nrow(out2)){
        cat(out2[j,], "\n", file = corr, sep = "\t")
      }
      close(corr)
    }
  }
}
################################################################################
# writeBoot end                                                                #
################################################################################
#
#
#
#
#
################################################################################
# gpSampler                                                                    #
################################################################################
#' @export
gpSampler <- function(infile = NULL, samp_size = 10, outfile = NULL){
  dat <- fileReader(infile)
  rownames(dat) <- NULL
  dat <- as.matrix(dat)
  # determine genepop format
  p1 <- which(toupper(dat[,1]) == "POP")[1] + 1
  gp <- as.numeric(names(sort(-table(sapply(dat[p1, - 1], nchar)/2)))[1])
  dat <- as.data.frame(dat)
  rawData <- readGenepop(dat, gp = gp)  
  # resample
  if (length(samp_size) == 1){
    samp_size <- rep(samp_size, rawData$npops)
  }
  idx <- lapply(seq_along(samp_size), function(i){
    sample(samp_size[i], size = samp_size[i], replace = FALSE)
  })
  pop_list <- lapply(seq_along(idx), function(i){
    samp <- rawData$pop_list[[i]][idx[[i]], ]
    blnk <- rep("\t", ncol(samp))
    return(rbind(blnk, samp))
  })
  
  ind_vectors <- lapply(seq_len(rawData$npops), function(i){
    return(c("POP", paste(rep("pop_", samp_size[i]), i, "_", 
                          1:samp_size[i], ",", sep = "")))
  })
  pre_data <- matrix(rep("\t", ((rawData$nloci + 1) * (rawData$nloci + 1))),
                     ncol = (rawData$nloci + 1))
  pre_data[1, ] <- c("Title", rep("\t", rawData$nloci))
  for(i in 2:(rawData$nloci+1)){
    pre_data[i,1] <- rawData$loci_names[(i-1)]
  }
  # add pop to ind_vectors and pop_list
  popOut <- do.call("rbind", pop_list)
  indVect <- do.call("c", ind_vectors)
  output <- rbind(pre_data, cbind(indVect, popOut))
  rownames(output) <- NULL
  if(gp == 3){
    output[is.na(output)] <- "000000"
  } else {
    output[is.na(output)] <- "0000"
  }
  #bs_data_file<-data.frame(bs_data_file)
  out <- file(paste(outfile,".gen",sep=""), "w")
  for(i in 1:nrow(output)){
    if(i == nrow(output)){
      cat(output[i,], file = out,  sep="\t")
    } else {
      cat(output[i,], "\n", file = out,  sep="\t") 
    }
  }
  close(out)
}
################################################################################
# gpSampler                                                                    #
################################################################################
#
#
#
#
#
#
################################################################################
# polyIn                                                                       #
################################################################################
#' A function for calculating informativeness for the inference of ancestry
#' for loci of any ploidy.
#' Kevin Keenan 2014
#' @export 
polyIn <- function(infile = NULL, pairwise = FALSE, parallel = FALSE){
  if(is.null(infile)){
    stop("Please provide and input file!")
  }
  polyReader <- function(infile, parallel){
    # read data
    fastScan <- function(fname) {
      s <- file.info(fname)$size
      buf <- readChar(fname, s, useBytes = TRUE)
      return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
    }
    dat <- fastScan(infile)
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
    # npops
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    npops <- length(popLoc)
    no_col <- popLoc[1] - 1
    # get genotypes
    strt <- popLoc + 1
    ends <- c(popLoc[-1] - 1, length(dat))
    genoRet <- function(strt, ends, x){
      out <- strsplit(x[strt:ends], split = "\\s+")
      return(do.call("rbind", out))
    }
    genos <- mapply(genoRet, strt = strt, ends = ends, 
                    MoreArgs = list(x = dat), SIMPLIFY = FALSE)
    indNames <- lapply(genos, function(x){
      return(x[,1])
    })
    indNames <- do.call("c", indNames)
    genos <- lapply(genos, function(x){
      x[x == "-9"] <- NA
      return(x[,-1])
    })
    # ploidy estimation
    ploidy <- round(mean(nchar(na.omit(genos[[1]][,1]))))
    # convert genotypes into allele arrays
    if(parallel){
      library(parallel)
      cl <- makeCluster(detectCores())
      clusterExport(cl, "ploidy", envir = environment())
      allArr <- parLapply(cl, genos, function(x){
        alls <- strsplit(x, split = "")
        # fix NAs
        alls <- lapply(alls, function(y){
          if(any(is.na(y))){
            y <- rep(NA, ploidy)
            return(y)
          } else {
            return(y)
          }
        })
        # alls array
        alls <- matrix(alls, ncol = ncol(x), nrow = nrow(x))
        alls <- lapply(1:ncol(alls), function(i){
          do.call("rbind", alls[,i])
        })
        af <- lapply(alls, function(al){
          ct <- table(al)
          return(as.vector(ct/sum(ct)))
        })
      })
      # create allele frequency matrices
      afMat <- parLapply(cl, 1:length(locs), function(i){
        dat <- lapply(allArr, "[[", i)
        return(do.call("cbind", dat))
      })
      stopCluster(cl)
    } else {
      allArr <- lapply(genos, function(x){
        alls <- strsplit(x, split = "")
        # fix NAs
        alls <- lapply(alls, function(y){
          if(any(is.na(y))){
            y <- rep(NA, ploidy)
            return(y)
          } else {
            return(y)
          }
        })
        # alls array
        alls <- matrix(alls, ncol = ncol(x), nrow = nrow(x))
        alls <- lapply(1:ncol(alls), function(i){
          do.call("rbind", alls[,i])
        })
        af <- lapply(alls, function(al){
          ct <- table(al)
          return(as.vector(ct/sum(ct)))
        })
      })
      # create allele frequency matrices
      afMat <- lapply(1:length(locs), function(i){
        dat <- lapply(allArr, "[[", i)
        return(do.call("cbind", dat))
      })
    }
    
    names(afMat) <- locs
    return(afMat)
  }
  
  inFunc <- function(af, pw = FALSE){
    if(pw){
      combs <- combn(ncol(af), 2)
      afComb <- lapply(1:ncol(combs), function(i){
        return(af[,combs[,i]])
      })
      Out <- lapply(afComb, function(x){
        sum(apply(x, 1, function(y){
          trm1 <- -mean(y) * log(mean(y))
          trm2 <- sum((y/length(y))*log(y))
          return(sum(sum(trm1 + trm2)))
        }))
      })
      inOut <- matrix(NA, ncol = ncol(af), nrow = ncol(af))
      for(i in 1:length(Out)){
        inOut[combs[2,i], combs[1,i]] <- Out[[i]]
      }
    } else {
      inOut <- apply(af, 1, function(x){
        trm1 <- -mean(x) * log(mean(x))
        trm2 <- sum((x/length(x))*log(x))
        return(sum(sum(trm1 + trm2)))
      })
      inOut <- sum(inOut)
    }
    return(inOut)
  }
  afs <- polyReader(infile, parallel)
  locNames <- names(afs)
  if(parallel){
    library(parallel)
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("inFunc", "pairwise"), envir = environment())
    ins <- parLapply(cl, afs, inFunc, pw = pairwise)
    stopCluster(cl)
  } else {
    ins <- lapply(afs, inFunc, pw = pairwise)
  }
  if(!pairwise){
    ins <- unlist(ins)
    names(ins) <- locNames
  } else {
    names(ins) <- locNames
  }
  return(ins)
}
################################################################################
# polyIn                                                                       #
################################################################################
#
#
#
#
#
#
#
#
################################################################################
# snp2gp: a function for converting SNP data to genepop format                 #
################################################################################
snp2gp <- function(infile, prefix_length = 2){
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  if(is.list(infile)){
    infile <- as.matrix(infile)
    dat <- sapply(1:nrow(infile), function(i){
      out <- paste(infile[i,], collapse = "\t")
      return(out)
    })
    dat <- c(paste(colnames(infile), collapse = "\t"), dat)
  } else {
    dat <- fastScan(infile)
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
  pre <- strsplit(infile, split = "\\.")[[1]][1]
  pop_list <- cbind(do.call("c", indNames),
                    do.call("rbind", pop_list))
  pop_list <- rbind(c(paste(pre, "-converted", sep = ""), 
                      rep(NA, nloci)),
                    c(c(paste(locs[1:(nloci-1)], ",", sep = ""), 
                        locs[nloci]), NA), pop_list)
  pop_list[is.na(pop_list)] <- "\t"
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
################################################################################
# end snp2gp                                                                   #
################################################################################
#
#
#
#
#
#
################################################################################
# NEW STYLE FUNCTIONS                                                          #
################################################################################
## Uncomment for non c++ version
#' ### diffCalcRcpp: an Rcpp version of diveRsity::diffCalc
#' 
#' __Kevin Keenan__ (2014)
#' 
#' 
#' ```{r echo = FALSE, message = TRUE, cache = FALSE}
#'  library(rCharts, warn = FALSE)
#'  knitr::opts_knit$set(self.contained=TRUE)
#'  knitr::opts_chunk$set(tidy = FALSE, message = TRUE, echo=TRUE)
#' ```
#' 
#' This code will be the basis for a faster more efficinet pairwise bootstrap
#' routine. The major difference will be that all of the data will be resampled
#' n time and the pairwise statistics will be calculated following this step.
#' The original method pairs the data then bootstraps each pair seperatly. This
#' is much less efficient.
#' 
#' __Kevin Keenan__ (2014)
#' 
diffCalc <- function(infile = NULL, outfile = NULL, fst = FALSE, 
                     pairwise = FALSE, bs_locus = FALSE, 
                     bs_pairwise = FALSE, boots = NULL, 
                     para = FALSE){
  #' # Calculate diversity statistics functions
  #' 
  #' __Kevin Keenan__ (2014)
  #'
  #' ### D_Jost
  dCalc <- function(ht, hs, n = NULL){
    if(!is.null(n)){
      return(((ht-hs)/(1-hs)) * (n/(n-1)))
    } else {
      return(2* ((ht-hs)/(1-hs)))
    }
  }
  #' ### Gst (Nei)
  gstCalc <- function(ht, hs){
    return((ht - hs)/ht)
  }
  #' ### G'st (Hedrick)
  GstCalc <- function(ht, hs, n = NULL){
    if(!is.null(n)){
      htmax <- ((n - 1) + hs)/n
    } else {
      htmax <- (1+hs)/2
    }
    return(((ht-hs)/ht)/((htmax-hs)/htmax))
  }
  #' ### Locus and overall WC stats
  thetaCalc <- function(a, b, cdat){
    return(a/(a+b+cdat))
  }
  #' ### bias correction
  #' 
  #' Function for correcting the bias associated with bootstrapped
  #' statistics
  #' 
  #' __Kevin Keenan__ (2014)
  #' 
  # Calculate CIs (bias corrected)
  diffCalcbCor <- function(act, bs){
    if(is.matrix(bs)){
      mn <- rowMeans(bs, na.rm = TRUE) - act
      mn[is.nan(mn)] <- NA
      bc <- mapply(`-`, split(bs, row(bs)), mn)
      return(bc)
    } else {
      mn <- mean(bs, na.rm = TRUE) - act
      mn[is.nan(mn)] <- NA
      return(bs - mn)
    }
  }
  #' # Calculate the harmonic sample size for pair of populations
  #' 
  #' __Kevin Keenan__ (2014)
  diffCalcHarm <- function(idt, pw){
    ps <- apply(pw, 2, function(x){
      return((0.5*((idt[x[1]]^-1) + idt[x[2]]^-1))^-1)
    })
    return(ps)    
  }
  #   #data(Test_data, package = "diveRsity")
  #   #Test_data[is.na(Test_data)] <- ""
  #   #Test_data[Test_data == "0"] <- "000000"
  #   source("R/diffCalcbCor.R")
  #   source("R/diffCalcHarm.R")
  #   source("R/misc.R")
  #   infile <- "speed-tests/MyData.gen"#Test_data
  #   outfile <- "test"
  #   fst = TRUE
  #   pairwise = TRUE
  #   bs_locus = TRUE
  #   bs_pairwise = TRUE
  #   boots = 10
  #   para = FALSE
  
  
  bs <- boots
  # read infile
  ip <- rgp(infile)
  # define some parameters
  
  # pop sizes
  ps = sapply(ip$genos, function(x){
    dim(x)[1]
  })
  
  # npops
  np = ncol(ip$af[[1]])
  
  # pairwise matrix index
  pw <- combn(np, 2)
  
  # set up parallel cluster
  if(para){
    library(parallel)
    ncor <- detectCores()
  }
  
  # create resample indexes
  if(!is.null(bs)){
    idx <- lapply(1:bs, function(i){
      lapply(1:np, function(j){
        sample(ps[j], size = ps[j], replace = TRUE)
      })
    }) 
  }
  
  #' ### Calculate point estimates (locus and global)
  #' Calculate gst, Gst, Fst and D for loci (across samples) and 
  #' global (across loci and samples).
  #' 
  
  #####################################
  # Point calculations
  #####################################
  preStats <- statCalc(rsDat = ip$genos, al = ip$af, fst = fst, bs = FALSE)
  # identify loci with unscored samples and fix for glbWCcpp
  #lapply
  # tabMerge
  tabMerge <- function(...){
    ip <- unlist(list(...))
    idx <- names(ip) != "NA"
    ip <- ip[idx]
    out <- sapply(split(ip, names(ip)), sum)
    if(length(out) == 0L){
      ot <- NA
      names(ot) <- "NA"
      return(ot)
    } else {
      return(out)
    }
  }
  if(fst){
    # locus variance components (working)
    hsum <- lapply(preStats$hsum, tabMerge)
    varC <- mapply("glbWCcpp", hsum = hsum, af = preStats$alOut,
                   indtyp = preStats$indtyp, SIMPLIFY = FALSE)
    rm(hsum)
    # Locus F-statistics
    locFst <- sapply(varC, function(x){
      return(x$a/(x$a+x$b+x$c))
    })
    locFit <- sapply(varC, function(x){
      return(1 - (x$c/(x$a + x$b + x$c)))
    })
    locFis <- sapply(varC, function(x){
      return(1 - (x$c/(x$b + x$c)))
    })
    
    # global variance components
    glba <- mean(sapply(varC, "[[", "a"), na.rm = TRUE)
    glbb <- mean(sapply(varC, "[[", "b"), na.rm = TRUE)
    glbc <- mean(sapply(varC, "[[", "c"), na.rm = TRUE)
    
    # F-statistics
    glbTheta <- glba/(glba + glbb + glbc)
    glbFit <- 1 - (glbc/(glba + glbb + glbc))
    glbFis <- 1 - (glbc/(glbb + glbc)) 
  }
  
  # calculate harmonic mean sample size per locus
  hrm <- function(x){
    return(length(x)/(sum(1/x)))
  }
  harmLoc <- sapply(preStats$indtyp, hrm)
  # locus global stats (replace with Rcpp function)
  # calculate heterozygosities
  hthsLoc <- mapply(varFunc, af = preStats$alOut, sHarm = harmLoc, 
                    SIMPLIFY = FALSE)
  # Jost's D (loci + global) #
  dLoc <- dCalc(ht = sapply(hthsLoc, "[[", "htEst"),
                hs = sapply(hthsLoc, "[[", "hsEst"),
                n = np)
  mnD <- mean(dLoc, na.rm = TRUE)
  varD <- var(dLoc, na.rm = TRUE)
  dGlb <- 1/((1/mnD)+varD*(1/mnD)^3)
  
  # gst (loci + global) #
  gLoc <- gstCalc(sapply(hthsLoc, "[[", "htEst"),
                  sapply(hthsLoc, "[[", "hsEst"))
  gGlb <- gstCalc(mean(sapply(hthsLoc, "[[", "htEst"), na.rm = TRUE),
                  mean(sapply(hthsLoc, "[[", "hsEst"), na.rm = TRUE))
  
  # Gst (loci + global) #
  GLoc <- GstCalc(sapply(hthsLoc, "[[", "htEst"),
                  sapply(hthsLoc, "[[", "hsEst"),
                  n = np)
  GGlb <- GstCalc(mean(sapply(hthsLoc, "[[", "htEst"), na.rm = TRUE),
                  mean(sapply(hthsLoc, "[[", "hsEst"), na.rm = TRUE),
                  n = np)
  
  #####################################
  # END
  #####################################
  
  
  #' ### Create Est structure
  #' 
  
  #####################################
  # Arange points estimates for output
  #####################################
  # non mem prob
  if(!fst){
    est <- data.frame(loci = c(ip$locs, "Global"),
                      gst = round(c(gLoc, gGlb), 4),
                      Gst = round(c(GLoc, GGlb), 4),
                      D = round(c(dLoc, dGlb), 4))
    rownames(est) <- NULL
  } else {
    est <- data.frame(loci = c(ip$locs, "Global"),
                      gst = round(c(gLoc, gGlb), 4),
                      Gst = round(c(GLoc, GGlb), 4),
                      D = round(c(dLoc, dGlb), 4),
                      Fis = round(c(locFis, glbFis), 4),
                      Fst = round(c(locFst, glbTheta), 4),
                      Fit = round(c(locFit, glbFit), 4))
    rownames(est) <- NULL
  }
  
  #####################################
  # END
  #####################################
  
  
  #' #### Output format (preview)
  #' ```{r, cache=FALSE, echo=FALSE, results='asis'}
  #'  dt <- Datatables$new()
  #'  dt$addTable(est)
  #'  dt$print("table1", include_assets = TRUE, cdn = TRUE)
  #'  
  #'```
  #'
  #' <br>
  #'
  #' ### Calculate locus and global bootstraps
  #' Calculate the bootstrapped \(95\%\) Confidence intervals for locus and 
  #' global gst, Gst, Fst and D.
  #' 
  #' 
  #####################################
  # Bootstrap data
  ####################################
  if(bs_locus || bs_pairwise){
    # pre-statistics
    # to run the bootstraps
    # replace al with 0
    al <- lapply(ip$af, function(x){
      nms <- rownames(x)
      x[x != 0] <- 0
      rownames(x) <- nms
      return(x)
    })
    if(para){
      cl <- makeCluster(ncor)
      clusterExport(cl, c("myTab", "al"), envir = environment())
      bsDat <- parLapply(cl, idx, statCalc, rsDat = ip$genos, al = al, 
                         fst = fst)
      stopCluster(cl)
    } else {
      bsDat <- lapply(idx, statCalc, rsDat = ip$genos, al = al, fst = fst)
    }
    indtyp <- lapply(bsDat, "[[", "indtyp")
    # clean up
    rm(idx)
  }
  #####################################
  # END
  ####################################
  
  #####################################
  # Locus and global bootstrapping
  #####################################
  
  # added Rcpp
  if(bs_locus){
    if(para){
      cl <- makeCluster(ncor)
    }
    
    #####################################
    # W&C locus stats
    #####################################
    # Low memory overhead (Takes 12.1 sec for 1000 bs or "pw_test.gen")
    if(fst){
      if(para){
        clusterExport(cl, c("glbWCcpp", "tabMerge"), 
                      envir = environment())
        bsVarC <- parLapply(cl, bsDat, function(x){
          hsum <- lapply(x$hsum, tabMerge)
          stat <- mapply("glbWCcpp", hsum = hsum, af = x$alOut, 
                         indtyp = x$indtyp, SIMPLIFY = FALSE)
          bsFst <- sapply(stat, function(x){
            return(x$a/(x$a+x$b+x$c))
          })
          bsFit <- sapply(stat, function(x){
            return(1 - (x$c/(x$a + x$b + x$c)))
          })
          bsFis <- sapply(stat, function(x){
            return(1 - (x$c/(x$b + x$c)))
          })
          glba <- mean(sapply(stat, "[[", "a"), na.rm = TRUE)
          glbb <- mean(sapply(stat, "[[", "b"), na.rm = TRUE)
          glbc <- mean(sapply(stat, "[[", "cdat"), na.rm = TRUE)
          glbst <- glba/(glba + glbb + glbc)
          glbit <- 1 - (glbc/(glba + glbb + glbc))
          glbis <- 1 - (glbc/(glbb + glbc))
          list(bsFstLoc = bsFst, bsFitLoc = bsFit, bsFisLoc = bsFis, 
               bsFstAll = glbst, bsFitAll = glbit, bsFisAll = glbis)
        })
      } else {
        bsVarC <- lapply(bsDat, function(x){
          hsum <- lapply(x$hsum, tabMerge)
          stat <- mapply(glbWCcpp, hsum = hsum, af = x$alOut, 
                         indtyp = x$indtyp, SIMPLIFY = FALSE)
          bsFst <- sapply(stat, function(y){
            return(y$a / (y$a + y$b + y$c))
          })
          bsFit <- sapply(stat, function(y){
            return(1 - (y$c / (y$a + y$b + y$c)))
          })
          bsFis <- sapply(stat, function(y){
            return(1 - (y$c / (y$b + y$c)))
          })
          glba <- mean(sapply(stat, "[[", "a"), na.rm = TRUE)
          glbb <- mean(sapply(stat, "[[", "b"), na.rm = TRUE)
          glbc <- mean(sapply(stat, "[[", "c"), na.rm = TRUE)
          glbst <- glba/(glba + glbb + glbc)
          glbit <- 1 - (glbc/(glba + glbb + glbc))
          glbis <- 1 - (glbc/(glbb + glbc))
          list(bsFstLoc = bsFst, bsFitLoc = bsFit, bsFisLoc = bsFis, 
               bsFstAll = glbst, bsFitAll = glbit, bsFisAll = glbis)
        })
      }
      # Extract bootstrap statistics
      bsFstL <- sapply(bsVarC, "[[", "bsFstLoc")
      bsFisL <- sapply(bsVarC, "[[", "bsFisLoc")
      bsFitL <- sapply(bsVarC, "[[", "bsFitLoc")
      bsFstA <- sapply(bsVarC, "[[", "bsFstAll")
      bsFisA <- sapply(bsVarC, "[[", "bsFisAll")
      bsFitA <- sapply(bsVarC, "[[", "bsFitAll")
      # Calculate bias corrected bs values
      bsFstL <- diffCalcbCor(locFst, bsFstL)
      bsFisL <- diffCalcbCor(locFis, bsFisL)
      bsFitL <- diffCalcbCor(locFit, bsFitL)
      bsFstA <- diffCalcbCor(glbTheta, bsFstA)
      bsFisA <- diffCalcbCor(glbFis, bsFisA)
      bsFitA <- diffCalcbCor(glbFit, bsFitA)
      # Calculate 95% CIs
      LCI <- data.frame(Fst = apply(bsFstL, 2, quantile, probs = 0.025, 
                                    na.rm = TRUE),
                        Fis = apply(bsFisL, 2, quantile, probs = 0.025, 
                                    na.rm = TRUE),
                        Fit = apply(bsFitL, 2, quantile, probs = 0.025, 
                                    na.rm = TRUE))
      UCI <- data.frame(Fst = apply(bsFstL, 2, quantile, probs = 0.975, 
                                    na.rm = TRUE),
                        Fis = apply(bsFisL, 2, quantile, probs = 0.975, 
                                    na.rm = TRUE),
                        Fit = apply(bsFitL, 2, quantile, probs = 0.975, 
                                    na.rm = TRUE))
      # clean up
      #rm(bsFstL, bsFisL, bsFitL, bsVarC)
      #z <- gc(reset = TRUE)
      #rm(z)
    }
    #####################################
    # END W&C locus stats
    #####################################
    
    
    #####################################
    # Het based locus stats
    #####################################
    # Calculate heterozygosity based stats
    # harmonic sample sizes
    
    allHarm <- lapply(bsDat, function(x){
      sapply(x$indtyp, function(y){
        return(((1/length(y))*(sum(y^-1)))^-1)
      })
    })
    
    # add allHarm to bsDat
    listAdd <- function(bsDat, sHarm){
      bsDat$nHarm <- sHarm
      return(bsDat)
    }
    
    bsDat <- mapply(listAdd, bsDat = bsDat, sHarm = allHarm, SIMPLIFY = FALSE)
    
    # Very low memory over head (takes 1.6 sec for 1000 boostraps using 
    # "pw_test.gen")
    # Calculate stats
    if(para){
      clusterExport(cl, c("varFunc", "dCalc", "gstCalc", "GstCalc"),
                    envir = environment())
      hetVar <- parLapply(cl, bsDat, function(x){
        stat <- mapply("varFunc", af = x$alOut, sHarm = x$nHarm, 
                       SIMPLIFY = FALSE)
        ht <- sapply(stat, "[[", "htEst")
        hs <- sapply(stat, "[[", "hsEst")
        n <- ncol(x$alOut[[1]])
        bsDLoc <- dCalc(ht = ht, hs = hs, n = n)
        mnD <- mean(bsDLoc, na.rm = TRUE)
        varD <- var(bsDLoc, na.rm = TRUE)
        bsDAll <- 1/((1/mnD) + varD * (1/mnD)^3)
        bsgLoc <- gstCalc(ht = ht, hs = hs)
        bsgAll <- gstCalc(ht = mean(ht, na.rm = TRUE), 
                          hs = mean(hs, na.rm = TRUE))
        bsGLoc <- GstCalc(ht = ht, hs = hs, n = n)
        bsGAll <- GstCalc(ht = mean(ht, na.rm = TRUE),
                          hs = mean(hs, na.rm = TRUE), n = n)
        list(bsDLoc = bsDLoc, bsgLoc = bsgLoc, bsGLoc = bsGLoc, 
             bsDAll = bsDAll, bsgAll = bsgAll, bsGAll = bsGAll)
        })
    } else {
      hetVar <- lapply(bsDat, function(x){
        stat <- mapply("varFunc", af = x$alOut, sHarm = x$nHarm, 
                       SIMPLIFY = FALSE)
        ht <- sapply(stat, "[[", "htEst")
        hs <- sapply(stat, "[[", "hsEst")
        n <- ncol(x$alOut[[1]])
        bsDLoc <- dCalc(ht = ht, hs = hs, n = n)
        mnD <- mean(bsDLoc, na.rm = TRUE)
        varD <- var(bsDLoc, na.rm = TRUE)
        bsDAll <- 1/((1/mnD) + varD * (1/mnD)^3)
        bsgLoc <- gstCalc(ht = ht, hs = hs)
        bsgAll <- gstCalc(ht = mean(ht, na.rm = TRUE), 
                          hs = mean(hs, na.rm = TRUE))
        bsGLoc <- GstCalc(ht = ht, hs = hs, n = n)
        bsGAll <- GstCalc(ht = mean(ht, na.rm = TRUE),
                          hs = mean(hs, na.rm = TRUE), n = n)
        list(bsDLoc = bsDLoc, bsgLoc = bsgLoc, bsGLoc = bsGLoc, 
             bsDAll = bsDAll, bsgAll = bsgAll, bsGAll = bsGAll)
      })
   }
    
    # Extract bootstrap statistics
    bsDL <- sapply(hetVar, "[[", "bsDLoc")
    bsgL <- sapply(hetVar, "[[", "bsgLoc")
    bsGL <- sapply(hetVar, "[[", "bsGLoc")
    bsDA <- sapply(hetVar, "[[", "bsDAll")
    bsgA <- sapply(hetVar, "[[", "bsgAll")
    bsGA <- sapply(hetVar, "[[", "bsGAll")
    # Calculate bias corrected bs values
    bsDL <- diffCalcbCor(dLoc, bsDL)
    bsgL <- diffCalcbCor(gLoc, bsgL)
    bsGL <- diffCalcbCor(GLoc, bsGL)
    bsdA <- diffCalcbCor(dGlb, bsDA)
    bsgA <- diffCalcbCor(gGlb, bsgA)
    bsGA <- diffCalcbCor(GGlb, bsGA)
    
    # Calculate 95% CIs
    if(fst){
      # lower  
      LCI$D <- apply(bsDL, 2, quantile, probs = 0.025, na.rm = TRUE)
      LCI$gst <- apply(bsgL, 2, quantile, probs = 0.025, na.rm = TRUE)
      LCI$Gst <- apply(bsGL, 2, quantile, probs = 0.025, na.rm = TRUE)
      # upper
      UCI$D <- apply(bsDL, 2, quantile, probs = 0.975, na.rm = TRUE)
      UCI$gst <- apply(bsgL, 2, quantile, probs = 0.975, na.rm = TRUE)
      UCI$Gst <- apply(bsGL, 2, quantile, probs = 0.975, na.rm = TRUE)
    } else {
      # lower  
      LCI <- data.frame(D = apply(bsDL, 2, quantile, probs = 0.025, na.rm = TRUE))
      LCI$gst <- apply(bsgL, 2, quantile, probs = 0.025, na.rm = TRUE)
      LCI$Gst <- apply(bsGL, 2, quantile, probs = 0.025, na.rm = TRUE)
      # upper
      UCI <- data.frame(D = apply(bsDL, 2, quantile, probs = 0.975, na.rm = TRUE))
      UCI$gst <- apply(bsgL, 2, quantile, probs = 0.975, na.rm = TRUE)
      UCI$Gst <- apply(bsGL, 2, quantile, probs = 0.975, na.rm = TRUE)
    }
    
    #####################################
    # END: Het based locus stats
    #####################################
    
    # global CIs
    if(fst){
      glbBS <- data.frame(gst = bsgA, Gst = bsGA,  d = bsDA, fst = bsFstA, 
                          fit = bsFitA, fis = bsFisA)
      rm(bsFstA, bsFisA, bsFitA, bsDA, bsgA, bsGA)
    } else {
      glbBS <- data.frame(gst = bsgA, Gst = bsGA, d = bsDA) 
      rm(bsDA, bsgA, bsGA)
    }
    glbLCI <- apply(glbBS, 2, quantile, probs = 0.025, na.rm = TRUE)
    glbUCI <- apply(glbBS, 2, quantile, probs = 0.975, na.rm = TRUE)
    if(fst){
      glbOut <- data.frame(stat = c("gst", "Gst", "D", "Fst", "Fit", "Fis"),
                           actual = c(gGlb, GGlb, dGlb, glbTheta, glbFit, glbFis),
                           lower = glbLCI, upper = glbUCI, row.names = NULL) 
    } else {
      glbOut <- data.frame(stat = c("gst", "Gst", "D"),
                           actual = c(gGlb, GGlb, dGlb),
                           lower = glbLCI, upper = glbUCI, row.names = NULL) 
    }
   if(para){
     stopCluster(cl)
   }
  }
  
  
  #####################################
  # END
  #####################################
  
  popnms <- sapply(ip$indnms, "[[", 1)
  pwpops <- paste(popnms[pw[1,]], " vs ", popnms[pw[2,]])
  
  ######################################
  # Pairwise bootstrapping
  #####################################
  #' ### Calculate PW bootstrap statistics
  
  if(pairwise || bs_pairwise){
    #' ### Define heterozygosity based stats (Gst, G'st, D_Jost)
    
    # Calculate non-bootstrap parameters
    
    # calculate harmonic sample sizes from preStats
    preNharm <- lapply(preStats$indtyp, diffCalcHarm, pw = pw)
    # add to preStats
    preStats$sHarm <- preNharm
    # calculate heterozygosities
    nbshths <- mapply(pwHCalc, af = preStats$alOut, sHarm = preStats$sHarm,
                      MoreArgs = list(pw = pw-1), SIMPLIFY = FALSE)
    # convert to locus focus
    nbshths <- lapply(c("hsEst", "htEst"), function(x){
      lapply(nbshths, "[[", x)
    })
    names(nbshths) <- c("hsEst", "htEst")
    # calculate mean ht and hs
    hsMn <- vapply(1:ncol(pw), function(i){
      mean(vapply(nbshths$hsEst, "[[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
    }, FUN.VALUE = numeric(1))
    
    htMn <- vapply(1:ncol(pw), function(i){
      mean(vapply(nbshths$htEst, "[[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
    }, FUN.VALUE = numeric(1))
    
    # calculate non-bootstrap pairwise stats
    
    # Jost's D
    pwDLoc <- mapply(dCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwDall <- apply(pwDLoc, 1, function(x){
      mnD <- mean(x, na.rm = TRUE)
      vrD <- mean(x, na.rm = TRUE)
      return(1/((1/mnD) + vrD * (1/mnD)*3))
    })
    
    # gst
    pwgLoc <- mapply(gstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwgAll <- gstCalc(ht = htMn, hs = hsMn)
    
    # Gst
    pwGLoc <- mapply(GstCalc, ht = nbshths$htEst, hs = nbshths$hsEst)
    pwGAll <- GstCalc(ht = htMn, hs = hsMn)
    
    if(fst){
      # theta
      # calculate standard statistics (non-bootstrap)
      hsum <- lapply(preStats$hsum, function(x){
        op <- lapply(1:ncol(pw), function(i){
          return(tabMerge(x[pw[,i]]))
        })
        return(op)
      })
      pwVar <- mapply("pwWCcpp", hsum1 = hsum, af1 = preStats$alOut,
                      indtyp1 = preStats$indtyp, MoreArgs = list(pw = pw-1),
                      SIMPLIFY = FALSE)
      # convert to locus focus
      pwVar <- lapply(c("a", "b", "c"), function(x){
        return(lapply(pwVar, "[[", x))
      })
      names(pwVar) <- c("a", "b", "cdat")
      
      # loci
      pwFstLoc <- mapply(thetaCalc, a = pwVar$a, b = pwVar$b, cdat = pwVar$cdat)
      
      # Global
      mnA <- vapply(1:ncol(pw), function(i){
        mean(vapply(pwVar$a, "[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
      mnB <- vapply(1:ncol(pw), function(i){
        mean(vapply(pwVar$b, "[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
      mnC <- vapply(1:ncol(pw), function(i){
        mean(vapply(pwVar$cdat, "[", i, FUN.VALUE = numeric(1)), na.rm = TRUE)
      }, FUN.VALUE = numeric(1))
      pwFstAll <- mnA/(mnA + mnB + mnC)
    }
  }
  
  if(bs_pairwise){
    # set up parallel cluster
    if(para){
      cl <- makeCluster(ncor)
    }
    # calculate harmonic sample sizes (list[bs]:list[loc]:vector[pw])
    # low memory overhead (takes 2.7 sec for 1000 bs using "pw_test.gen")
    if(para){
      clusterExport(cl, c("diffCalcHarm", "pw"), envir = environment())
      sHarm <- parLapply(cl, indtyp, function(x){
        lapply(x, diffCalcHarm, pw = pw)
      })
      #stopCluster(cl)
    } else {
      sHarm <- lapply(indtyp, function(x){
        lapply(x, diffCalcHarm, pw = pw)
      })
    }
    rm(indtyp)
    # combine sHarm with bsDat
    listAdd <- function(bsDat, sHarm){
      bsDat$sHarm <- sHarm
      bsDat$nHarm <- NULL
      return(bsDat)
    }
    
    bsDat <- mapply(listAdd, bsDat = bsDat, sHarm = sHarm, SIMPLIFY = FALSE)
    rm(sHarm)
    
    # calculate data
    if(para){
      clusterExport(cl, c("pw", "pwHCalc"), envir = environment())
      hths <- parLapply(cl, bsDat, function(x){
        lapply(1:length(x$alOut), function(i){
          pwHCalc(x$alOut[[i]], sHarm = x$sHarm[[i]], pw = pw-1)
        })
        # mapply(pwHetCalc, af = x$alOut, sHarm = x$sHarm, 
        #       MoreArgs = list(pw = pw), SIMPLIFY = FALSE)
        #gc()
      })
    } else {
      hths <- lapply(bsDat, function(x){
        lapply(1:length(x$alOut), function(i){
          pwHCalc(x$alOut[[i]], sHarm = x$sHarm[[i]], pw = pw-1)
        })
        # mapply(pwHetCalc, af = x$alOut, sHarm = x$sHarm, 
        #       MoreArgs = list(pw = pw), SIMPLIFY = FALSE)
      })
    }
    
    # convert hths to stat focus
    hths <- lapply(hths, function(x){
      hsEst <- lapply(x, "[[", 1)
      htEst <- lapply(x, "[[", 2)
      #hs <- lapply(x, "[[", 3)
      #ht <- lapply(x, "[[", 4)
      return(list(hsEst = hsEst, htEst = htEst))#, hs = hs, ht = ht))
    })
    
    
    # calculate locus estimator (parallel is slower!)
    pwDLocbs <- sapply(hths, function(x){
      return(mapply(dCalc, ht = x$htEst, hs = x$hsEst, SIMPLIFY = TRUE))
    }, simplify = "array")
    # overall d bootstraps
    if(para){
      pwDAllbs <- parApply(cl, pwDLocbs, c(1,3), function(x){
        mn <- mean(x, na.rm = TRUE)
        vr <- var(x, na.rm = TRUE)
        return(1/((1/mn)+vr*(1/mn)^3))  
      })
    } else {
      pwDAllbs <- apply(pwDLocbs, c(1,3), function(x){
        mn <- mean(x, na.rm = TRUE)
        vr <- var(x, na.rm = TRUE)
        return(1/((1/mn)+vr*(1/mn)^3))  
      })
    }
    rm(pwDLocbs)
    
    # calculate locus estimtor
    #     pwgLocbs <- sapply(hths, function(x){
    #       return(mapply(gstCalc, ht = x$htEst, hs = x$hsEst, SIMPLIFY = TRUE))
    #     }, simplify = "array")
    # overall gst (returns pw values in rows and bootstraps reps in cols)
    pwgAllbs <- sapply(hths, function(x){
      ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
      hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
      return(gstCalc(ht, hs))
    }, simplify = "array")
    
    # calculate locus estimator
    #     pwGLocbs <- sapply(hths, function(x){
    #       return(mapply(GstCalc, ht = x$htEst, hs = x$hsEst, SIMPLIFY = TRUE))
    #     }, simplify = "array")
    # overall Gst (returns pw values in rows and bootstraps reps in cols)
    pwGAllbs <- sapply(hths, function(x){
      ht <- rowMeans(matrix(unlist(x$htEst), ncol = length(x$htEst)))
      hs <- rowMeans(matrix(unlist(x$hsEst), ncol = length(x$hsEst)))
      return(GstCalc(ht, hs))
    }, simplify = "array")
    
    #' ### Calculate Weir & Cockerham variance components 
    if(fst){
      # Calculate bootstrap statistics
      if(para){
        #cl <- makeCluster(detectCores())
        clusterExport(cl, c("pw", "pwWCcpp", "tabMerge"), envir = environment())
        wcVar <- parLapply(cl, bsDat, function(x){
          hsum <- lapply(x$hsum, function(y){
            op <- lapply(1:ncol(pw), function(i){
              return(tabMerge(y[pw[,i]]))
            })
            return(op)
          })
          return(mapply("pwWCcpp", hsum1 = hsum, indtyp1 = x$indtyp, af1 = x$alOut, 
                        MoreArgs = list(pw = pw-1), SIMPLIFY = FALSE))
        })
        #stopCluster(cl)
      } else {
        wcVar <- lapply(bsDat, function(x){
          hsum <- lapply(x$hsum, function(y){
            op <- lapply(1:ncol(pw), function(i){
              return(tabMerge(y[pw[,i]]))
            })
            return(op)
          })
          return(mapply("pwWCcpp", hsum1 = hsum, indtyp1 = x$indtyp, af1 = x$alOut, 
                        MoreArgs = list(pw = pw-1), SIMPLIFY = FALSE))
        })
      }
      
      # Calculate bootstrap theta
      #       pwFstLocbs <- sapply(wcVar, function(x){
      #         a <- lapply(x, "[[", "a")
      #         b <- lapply(x, "[[", "b")
      #         cdat <- lapply(x, "[[", "cdat")
      #         return(mapply(thetaCalc, a = a, b = b, cdat = cdat, 
      #                SIMPLIFY = TRUE))
      #       }, simplify = "array")
      pwFstAllbs <- sapply(wcVar, function(x){
        a <- rowMeans(sapply(x, "[[", "a"))
        b <- rowMeans(sapply(x, "[[", "b"))
        cdat <- rowMeans(sapply(x, "[[", "c"))
        return(thetaCalc(a, b, cdat))
      }, simplify = "array")
      # clean up
      rm(wcVar, bsDat)
      z <- gc()
      
      
      
      #################################
      # calculate WC bias corrected CIs
      #################################
      
      # loci
      #       pwFstLocbs <- mapply(diffCalcbCor, 
      #                            act = split(pwFstLoc, row(pwFstLoc)),
      #                            bs = lapply(1:dim(pwFstLocbs)[1], function(i){
      #                              return(pwFstLocbs[i,,])}), SIMPLIFY = "array")
      #       # transpose array
      #       pwFstLocbs <- aperm(pwFstLocbs)
      #       
      #       pwFstLCI <- apply(pwFstLocbs, c(1,2), quantile, prob = 0.025, 
      #                         na.rm = TRUE)
      #       pwFstUCI <- apply(pwFstLocbs, c(1,2), quantile, prob = 0.975, 
      #                         na.rm = TRUE)
      
      # global
      pwFstAllbs <- t(mapply(diffCalcbCor, act = pwFstAll, 
                             bs = split(pwFstAllbs, row(pwFstAllbs))))
      pwFstAllLCI <- apply(pwFstAllbs, 1, quantile, prob = 0.025, 
                           na.rm = TRUE)
      pwFstAllUCI <- apply(pwFstAllbs, 1, quantile, prob = 0.975, 
                           na.rm = TRUE)
      
      #################################
      # END
      #################################
    }
    
    #' ### Heterozygosity based CIs
    
    #################################
    # Calculate het based bs stats
    #################################
    
    # Jost's D  
    # loci
    #     pwDLocbs <- mapply(diffCalcbCor, act = split(pwDLoc, row(pwDLoc)),
    #                        bs = lapply(1:dim(pwDLocbs)[1], function(i){
    #                          return(pwDLocbs[i,,])}), SIMPLIFY = "array")
    #     pwDLocbs <- aperm(pwDLocbs)
    #     # CIs
    #     pwDLocLCI <- apply(pwDLocbs, c(1,2), quantile, prob = 0.025, 
    #                        na.rm = TRUE)
    #     pwDLocUCI <- apply(pwDLocbs, c(1,2), quantile, prob = 0.975, 
    #                        na.rm = TRUE)
    
    # global
    pwDAllbs <- t(mapply(diffCalcbCor, act = pwDall, 
                         bs = split(pwDAllbs, row(pwDAllbs))))
    # CIs
    pwDAllLCI <- apply(pwDAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwDAllUCI <- apply(pwDAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    # gst
    # loci
    #     pwgLocbs <- mapply(diffCalcbCor, act = split(pwgLoc, row(pwgLoc)),
    #                        bs = lapply(1:dim(pwgLocbs)[1], function(i){
    #                          return(pwgLocbs[i,,])}), SIMPLIFY = "array")
    #     pwgLocbs <- aperm(pwgLocbs)
    #     # CIs
    #     pwgLocLCI <- apply(pwgLocbs, c(1,2), quantile, prob = 0.025, 
    #                        na.rm = TRUE)
    #     pwgLocUCI <- apply(pwgLocbs, c(1,2), quantile, prob = 0.975, 
    #                        na.rm = TRUE)
    
    # global
    pwgAllbs <- t(mapply(diffCalcbCor, act = pwgAll, 
                         bs = split(pwgAllbs, row(pwgAllbs))))
    # CIs
    pwgAllLCI <- apply(pwgAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwgAllUCI <- apply(pwgAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    # Gst
    # loci
    #     pwGLocbs <- mapply(diffCalcbCor, act = split(pwGLoc, row(pwGLoc)),
    #                        bs = lapply(1:dim(pwGLocbs)[1], function(i){
    #                          return(pwGLocbs[i,,])}), SIMPLIFY = "array")
    #     pwGLocbs <- aperm(pwGLocbs)
    #     # CIs
    #     pwGLocLCI <- apply(pwGLocbs, c(1,2), quantile, prob = 0.025, 
    #                        na.rm = TRUE)
    #     pwGLocUCI <- apply(pwGLocbs, c(1,2), quantile, prob = 0.975, 
    #                        na.rm = TRUE)
    
    # global
    pwGAllbs <- t(mapply(diffCalcbCor, act = pwGAll, 
                         bs = split(pwGAllbs, row(pwGAllbs))))
    # CIs
    pwGAllLCI <- apply(pwGAllbs, 1, quantile, prob = 0.025, na.rm = TRUE)
    pwGAllUCI <- apply(pwGAllbs, 1, quantile, prob = 0.975, na.rm = TRUE)
    
    #################################
    # END
    #################################
    # stop cluster
    if(para){
      stopCluster(cl) 
    }
  }
  
  
  
  
  #' ### organise outputs
  # est is already available from above
  
  #################################
  # Global bootstrap results
  #################################
  op <- list(std_stats = est)
  if(bs_locus){
    op$global_bs <- data.frame(stat = glbOut[,1], round(glbOut[,-1], 4))
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Locus Confidence intervals
  #################################
  if(bs_locus){
    if(fst){
      statnms <- c("gst", "Gst", "D", "Fst", "Fis", "Fit")
    } else {
      statnms <- c("gst", "Gst", "D")
    }
    locCI <- lapply(statnms, function(x){
      return(data.frame(locus = ip$locs,
                        actual = est[-nrow(est) ,x],
                        lower = round(LCI[, x], 4),
                        upper = round(UCI[, x], 4),
                        row.names = NULL))
    })
    names(locCI) <- statnms
    op$bs_locus <- locCI
    
    rm(locCI)
    z <- gc(reset = TRUE)
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Pairwise stats
  #################################
  if(pairwise){
    # locus
    if(fst){
      # gst
      pwgLoc <- data.frame(round(t(pwgLoc), 4))
      dimnames(pwgLoc) <- list(ip$locs, pwpops)
      op$pw_locus$gst <- pwgLoc
      rm(pwgLoc)
      # Gst
      pwGLoc <- data.frame(round(t(pwGLoc), 4))
      dimnames(pwGLoc) <- list(ip$locs, pwpops)
      op$pw_locus$Gst <- pwGLoc
      rm(pwGLoc)
      # D
      pwDLoc <- data.frame(round(t(pwDLoc), 4))
      dimnames(pwDLoc) <- list(ip$locs, pwpops)
      op$pw_locus$D <- pwDLoc
      rm(pwDLoc)
      # Fst
      pwFstLoc <- data.frame(round(t(pwFstLoc), 4))
      dimnames(pwFstLoc) <- list(ip$locs, pwpops)
      op$pw_locus$Fst <- pwFstLoc
      rm(pwFstLoc)
    } else {
      # gst
      pwgLoc <- data.frame(round(t(pwgLoc), 4))
      dimnames(pwgLoc) <- list(ip$locs, pwpops)
      op$pw_locus$gst <- pwgLoc
      rm(pwgLoc)
      # Gst
      pwGLoc <- data.frame(round(t(pwGLoc), 4))
      dimnames(pwGLoc) <- list(ip$locs, pwpops)
      op$pw_locus$Gst <- pwGLoc
      rm(pwGLoc)
      # D
      pwDLoc <- data.frame(round(t(pwDLoc), 4))
      dimnames(pwDLoc) <- list(ip$locs, pwpops)
      op$pw_locus$D <- pwDLoc
      rm(pwDLoc)
    }
    # global
    if(fst){
      # gst
      op$pairwise <- list(gst = matrix(NA, nrow = np, ncol = np))
      op$pairwise$gst[lower.tri(op$pairwise$gst)] <- round(pwgAll, 4)
      dimnames(op$pairwise$gst) <- list(popnms, popnms)
      #rm(pwgAll)
      # Gst
      op$pairwise$Gst <- op$pairwise$gst
      op$pairwise$Gst[lower.tri(op$pairwise$Gst)] <- round(pwGAll, 4)
      dimnames(op$pairwise$Gst) <- list(popnms, popnms)
      #rm(pwGAll)
      # D
      op$pairwise$D <- op$pairwise$gst
      op$pairwise$D[lower.tri(op$pairwise$D)] <- round(pwDall, 4)
      dimnames(op$pairwise$D) <- list(popnms, popnms)
      #rm(pwDall)
      # Fst
      op$pairwise$Fst <- op$pairwise$gst
      op$pairwise$Fst[lower.tri(op$pairwise$Fst)] <- round(pwFstAll, 4)
      dimnames(op$pairwise$Fst) <- list(popnms, popnms)
      #rm(pwFstAll)
    } else {
      # gst
      op$pairwise <- list(gst = matrix(NA, nrow = np, ncol = np))
      op$pairwise$gst[lower.tri(op$pairwise$gst)] <- round(pwgAll, 4)
      dimnames(op$pairwise$gst) <- list(popnms, popnms)
      #rm(pwgAll)
      # Gst
      op$pairwise$Gst <- op$pairwise$gst
      op$pairwise$Gst[lower.tri(op$pairwise$Gst)] <- round(pwGAll, 4)
      dimnames(op$pairwise$Gst) <- list(popnms, popnms)
      #rm(pwGAll)
      # D
      op$pairwise$D <- op$pairwise$gst
      op$pairwise$D[lower.tri(op$pairwise$D)] <- round(pwDall, 4)
      dimnames(op$pairwise$D) <- list(popnms, popnms)
      #rm(pwDall)
    }
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Pairwise CI
  #################################
  if(bs_pairwise){
    if(fst){
      # gst
      op$bs_pairwise <- list(gst = data.frame(populations = pwpops,
                                              actual = round(pwgAll, 4),
                                              lower = round(pwgAllLCI, 4),
                                              upper = round(pwgAllUCI, 4),
                                              row.names = NULL))
      rm(pwgAll, pwgAllLCI, pwgAllUCI)
      # Gst
      op$bs_pairwise$Gst <- data.frame(populations = pwpops,
                                       actual = round(pwGAll, 4),
                                       lower = round(pwGAllLCI, 4),
                                       upper = round(pwGAllUCI, 4),
                                       row.names = NULL)
      rm(pwGAll, pwGAllLCI, pwGAllUCI)
      # D
      op$bs_pairwise$D <- data.frame(populations = pwpops,
                                     actual = round(pwDall, 4),
                                     lower = round(pwDAllLCI, 4),
                                     upper = round(pwDAllUCI, 4),
                                     row.names = NULL)
      rm(pwDall, pwDAllLCI, pwDAllUCI)
      # Fst
      op$bs_pairwise$Fst <- data.frame(populations = pwpops,
                                       actual = round(pwFstAll, 4),
                                       lower = round(pwFstAllLCI, 4),
                                       upper = round(pwFstAllUCI, 4),
                                       row.names = NULL)
      rm(pwFstAll, pwFstAllLCI, pwFstAllUCI)
    } else {
      # gst
      op$bs_pairwise <- list(gst = data.frame(populations = pwpops,
                                              actual = round(pwgAll, 4),
                                              lower = round(pwgAllLCI, 4),
                                              upper = round(pwgAllUCI, 4),
                                              row.names = NULL))
      rm(pwgAll, pwgAllLCI, pwgAllUCI)
      # Gst
      op$bs_pairwise$Gst <- data.frame(populations = pwpops,
                                       actual = round(pwGAll, 4),
                                       lower = round(pwGAllLCI, 4),
                                       upper = round(pwGAllUCI, 4),
                                       row.names = NULL)
      rm(pwGAll, pwGAllLCI, pwGAllUCI)
      # D
      op$bs_pairwise$D <- data.frame(populations = pwpops,
                                     actual = round(pwDall, 4),
                                     lower = round(pwDAllLCI, 4),
                                     upper = round(pwDAllUCI, 4),
                                     row.names = NULL)
      rm(pwDall, pwDAllLCI, pwDAllUCI)
      z <- gc(reset = TRUE)
    }
  }
  #################################
  # END
  #################################
  
  
  #################################
  # Write results to file
  #################################
  if(!is.null(outfile)){
    # set up an output folder
    opf <- paste(getwd(), "/", outfile, "-[diffCalc]/", sep = "")
    dir.create(opf, showWarnings = FALSE)
    
    outnms <- names(op)
    # define a write function
    out <- sapply(outnms, function(x){
      if(x == "std_stats" || x == "global_bs"){
        ot <- paste(colnames(op[x][[1]]), collapse = "\t")
        preot <- apply(op[x][[1]], 1, paste, collapse = "\t")
        ot <- c(ot, preot)
        #fl <- file(paste(x, ".txt", sep = ""), "w")
        #cat(ot, sep = "\n", file = fl)
        #close(fl)
        writeLines(paste(ot, collapse = "\n"), 
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "pairwise"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          dat <- op[x][[1]][y][[1]]
          dat[is.na(dat)] <- ""
          dimnames(dat) <- list(NULL, NULL)
          opt <- apply(dat, 1, paste0, collapse = "\t", na.rm = "")
          popnmsOut <- paste(popnms, "\t", sep = "")
          opt <- mapply(paste, popnmsOut, opt, 
                        MoreArgs = list(collapse = "\t")) 
          opt <- c(y, "", paste("pops", paste(popnms, collapse = "\t"), 
                                sep = "\t"), opt, "")
          return(opt)
        })
        if(fst){
          ot <- c("Pairwise stats", "Fst = Weir & Cockerham's theta, (1984)",
                  "D = Jost, (2008)", "gst = Nei & Chesser, (1983)",
                  "Gst = Hedrick, (2005)", "", unlist(ot))
        } else {
          ot <- c("Pairwise stats", "D = Jost, (2008)", 
                  "gst = Nei & Chesser, (1983)",
                  "Gst = Hedrick, (2005)", "", 
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"), 
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "bs_locus"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          ot1 <- c("", y, "", paste(colnames(op[x][[1]][y][[1]]), 
                                    collapse = "\t"))
          ot2 <- apply(op[x][[1]][y][[1]], 1, paste, collapse = "\t")
          return(c(ot1, ot2))
        })
        if(fst){
          ot <- c("Locus 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "D = Jost, 2008",
                  "Fst = Weir & Cockerham, 1984",
                  "Fis = Weir & Cockerham, 1984",
                  "Fit = Weir & Cockerham, 1984", "",
                  unlist(ot))
        } else {
          ot <- c("Locus 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "D = Jost, 2008",
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"), 
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "pw_locus"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          ot1 <- c("", y, "", paste("Loci", paste(colnames(op[x][[1]][y][[1]]), 
                                                  collapse = "\t"), sep = "\t"))
          ot2 <- apply(op[x][[1]][y][[1]], 1, paste, collapse = "\t")
          ot2 <- mapply(`paste`, rownames(op[x][[1]][y][[1]]), ot2,
                        MoreArgs = list(sep = "\t"))
          return(c(ot1, ot2))
        })
        if(fst){
          ot <- c("Locus Pairwise estimates", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "D = Jost, 2008",
                  "Fst = Weir & Cockerham, 1984",
                  unlist(ot))
        } else {
          ot <- c("Locus Pairwise estimates", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "D = Jost, 2008",
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"),
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      } else if(x == "bs_pairwise"){
        statnms <- names(op[x][[1]])
        ot <- lapply(statnms, function(y){
          ot1 <- c("", y, "", paste(colnames(op[x][[1]][y][[1]]), 
                                    collapse = "\t"))
          ot2 <- apply(op[x][[1]][y][[1]], 1, paste, collapse = "\t")
          return(c(ot1, ot2))
        })
        if(fst){
          ot <- c("Pairwise 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "D = Jost, 2008",
                  "Fst = Weir & Cockerham, 1984",
                  unlist(ot))
        } else {
          ot <- c("Pairwise 95% CIs", "",
                  "gst = Nei & Chesser, 1983",
                  "Gst = Hedrick, 2005",
                  "D = Jost, 2008",
                  unlist(ot))
        }
        writeLines(paste(ot, collapse = "\n"), 
                   paste(opf, x, ".txt", sep = ""))
        ot <- NULL
      }
    })
    rm(out)
    #z <- gc(reset = TRUE)
  }
  #################################
  # END
  #################################
  
  return(op)
}
################################################################################
# end diffCalc function                                                        #
################################################################################
################################################################################
# New inCalc                                                                   #
################################################################################
#
#
#
#
#' New inCalc function for diveRsity package
#' 
#' Kevin Keenan 2014
#' 
# infile <- "Test.txt"
# bootstraps = 10
# pairwise = FALSE
# outfile <- "out"
# parallel = TRUE
# xlsx <- FALSE
inCalc <- function(infile = NULL, outfile = NULL, pairwise = FALSE, 
                   xlsx = FALSE, bootstraps = NULL, parallel = FALSE){
  #source("rgp.R")
  #source("inFunc.R")
  # calculate basic information
  alf <- rgp(infile)
  locs <- alf$locs
  pnms <- sapply(alf$indnms, function(x){return(x[1])})
  af <- lapply(alf$af, function(x){
    colnames(x) <- pnms
    return(x)
  })
  names(af) <- locs
  if(pairwise){
    inStatPW <- lapply(af, inFunc, pw = TRUE) 
  }
  inStatGLB <- lapply(af, inFunc, pw = FALSE)
  ####-- Global bootstrap --####
  if(!is.null(bootstraps)){
    # population sizes
    ps <- sapply(alf$genos, function(x){
      return(dim(x)[1])
    })
    # generate resample indexes
    idx <- lapply(1:bootstraps, function(x){
      lapply(ps, function(y){
        sample(y, y, replace = TRUE)
      })
    })
    # define a bootstrap function
    paraFunc <- function(genos, idx, af){
      # generate resamples
      rsFun <- function(x, y){
        return(x[y,,])
      }
      rsDat <- mapply(rsFun, x = genos, y = idx, SIMPLIFY = FALSE)
      
      # calculate allele frequecies
      alf <- lapply(rsDat, function(x){
        apply(x, 2, function(y){
          table(c(y[,1], y[,2]))/(length(na.omit(y[,1]))*2)
        })
      })
      # count loci
      nloci <- dim(rsDat[[1]])[2]
      # sort frequencies by loci
      locAl <- lapply(1:nloci, function(i){
        lapply(alf, "[[", i)
      })
      
      alSort <- function(x, y){
        idx <- lapply(x, function(z){
          match(names(z), rownames(y))
        })
        for(i in 1:length(idx)){
          y[idx[[i]], i] <- x[[i]]
        }
        return(y)
      }
      
      # generate allele frequency output
      alOut <- mapply(alSort, x = locAl, y = af, SIMPLIFY = FALSE)
      return(alOut)
    }
    
    if(parallel){
      library(parallel)
      cl <- makeCluster(detectCores())
      clusterExport(cl, c("inFunc", "paraFunc", "alf", "pairwise"), 
                    envir = environment())
      glbbs <- parLapply(cl, idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = FALSE))
      })
      if(!pairwise){
        stopCluster(cl)
      }
    } else {
      glbbs <- lapply(idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = FALSE))
      })
    }
    ####-- calculate global 95% CIs --####
    glbbs <- sapply(glbbs, function(x){
      unlist(x)
    })
    glbLowCI <- apply(glbbs, 1, quantile, probs = 0.025)
    glbUpCI <- apply(glbbs, 1, quantile, probs = 0.975) 
  }
  # organise data
  # global In
  glbDat <- data.frame(Locus = locs, Global_In = round(unlist(inStatGLB), 4))
  rownames(glbDat) <- NULL
  if(!is.null(bootstraps)){
    glbDat$lower_ci <- round(glbLowCI, 4)
    glbDat$upper_ci <- round(glbUpCI, 4)
  }
  ####-- Pairwise data --####
  if(pairwise){
    # pairwise In
    pw <- combn(ncol(af[[1]]), 2)
    pwmat <- round(do.call("rbind", inStatPW), 4)
    opHeader <- paste(colnames(af[[1]])[pw[1,]], " vs ", 
                      colnames(af[[1]])[pw[2,]], sep = "")
    colnames(pwmat) <- opHeader
    
    pwDat <- data.frame(Locus = locs)
    pwDat <- cbind(pwDat, pwmat)
    rownames(pwDat) <- NULL
    row1 <- paste(colnames(pwDat), collapse = "\t")
  }
  ### Bootstrap code ###
  if(!is.null(bootstraps) && pairwise){
    if(parallel){
      system.time({
        inbs <- parLapply(cl, idx, function(x){
          af <- paraFunc(alf$genos, x, alf$af)
          return(lapply(af, inFunc, pw = pairwise))
        })
        stopCluster(cl)
      })
    } else {
      inbs <- lapply(idx, function(x){
        af <- paraFunc(alf$genos, x, alf$af)
        return(lapply(af, inFunc, pw = pairwise))
      })
    }
    ####-- Calculate pairwise CIs --####
    inbs <- sapply(inbs, function(x){
      return(do.call("rbind", x))
    }, simplify = "array")
    # calculate CIs
    lowCI <- round(as.data.frame(apply(inbs, c(1,2), quantile, probs = 0.025)),
                   4)
    upCI <- round(as.data.frame(apply(inbs, c(1,2), quantile, probs = 0.975)),
                  4)
    # add names
    lowCI <- cbind(glbDat$Locus, lowCI)
    upCI <- cbind(glbDat$Locus, upCI)
    rownames(lowCI) <- NULL
    colnames(lowCI) <- c("Locus", opHeader)
    rownames(upCI) <- NULL
    colnames(upCI) <- c("Locus", opHeader)
  }
  
  
  ####-- Write the data to file --####
  
  # Set up directory
  outDir <- paste(getwd(), "/", outfile, "-[inCalc]/", sep = "")
  if(!file.exists(outDir)){
    dir.create(path = outDir, showWarnings = FALSE)
  }
  if(xlsx){
    if(!is.element("xlsx",installed.packages()[,1])){
      stop("The 'xlsx' package must be installed. Writing to text...")
    }
    require("xlsx")
    outF <- paste(outDir, "outfile-[in.Calc].xlsx", sep = "")
    write.xlsx(glbDat, file = outF, sheetName = "Global In", 
               col.names = TRUE, row.names = FALSE, append = FALSE)
    if(pairwise){
      write.xlsx(pwDat, file = outF, sheetName = "Pairwise In", 
                 col.names = TRUE, row.names = FALSE, append = TRUE)
    }
    if(!is.null(bootstraps)){
      write.xlsx(lowCI, file = outF, sheetName = "Lower CI (PW)", 
                 col.names = TRUE, row.names = FALSE, append = TRUE)
      write.xlsx(upCI, file = outF, sheetName = "Upper CI (PW)", 
                 col.names = TRUE, row.names = FALSE, append = TRUE)
    }
  } else {
    rn <- paste(colnames(glbDat), collapse = "\t")
    glbDatOut <- apply(glbDat, 1, paste, collapse = "\t")
    glbDatOut <- c(rn, glbDatOut)
    of1 <- paste(outDir, "Global-[in.Calc].txt", sep = "")
    fl1 <- file(of1, "w")
    for(i in 1:length(glbDatOut)){
      cat(glbDatOut[i], sep = "\n", file = fl1)
    }
    close(fl1)
    if(pairwise){
      pwDatOut <- apply(pwDat, 1, paste, collapse = "\t")
      pwDatOut <- c(row1, pwDatOut)
      of2 <- paste(outDir, "pairwise-[in.Calc].txt", sep = "")
      fl2 <- file(of2, "w")
      for(i in 1:length(pwDatOut)){
        cat(pwDatOut[i], sep = "\n", file = fl2)
      }
      close(fl2)
    }
    if(!is.null(bootstraps) && pairwise){
      of3 <- paste(outDir, "Lower_CI-[in.Calc].txt", sep = "")
      fl3 <- file(of3, "w")
      of4 <- paste(outDir, "Upper_CI-[in.Calc].txt", sep = "")
      fl4 <- file(of4, "w")
      lowCIout <- c(row1, apply(lowCI, 1, paste, collapse = "\t"))
      upCIout <- c(row1, apply(upCI, 1, paste, collapse = "\t"))
      for(i in 1:length(lowCIout)){
        cat(lowCIout[i], sep = "\n", file = fl3)
        cat(upCIout[i], sep = "\n", file = fl4)
      }
      close(fl3)
      close(fl4)
    }
  }
  ####-- Function outputs --####
  if(!pairwise && is.null(bootstraps)){
    list(global = glbDat)
  } else if(pairwise && is.null(bootstraps)){
    list(global = glbDat,
         pairwise = pwDat)
  } else if (pairwise && !is.null(bootstraps)){
    list(global = glbDat,
         pairwise = pwDat,
         lower_CI = lowCI,
         upper_CI = upCI)    
  } else {
    list(global = glbDat)
  }
}
#' actual inCalc sub-funciton
#' 
#' Kevin Keenan 2014
inFunc <- function(af, pw = FALSE){
  if(pw){
    combs <- combn(ncol(af), 2)
    afComb <- lapply(1:ncol(combs), function(i){
      return(af[,combs[,i]])
    })
    np <- length(unique(c(combs[1,], combs[2,])))
    out <- apply(combs, 2, function(x){
      y <- af[,x]
      y <- y[rowSums(y) != 0,]
      inOut <- apply(y, 1, function(z){
        p_j <- sum(z)/length(z)
        trm1 <- (-p_j*(log(p_j)))
        lgx <- log(z)
        lgx[is.infinite(lgx)] <- 0
        trm2 <- sum(z/(length(z))*lgx) 
        return(trm1 + trm2)
      })
      return(sum(inOut))
    })
    return(out)
  } else {
    inOut <- apply(af, 1, function(x){
      p_j <- sum(x)/length(x)
      trm1 <- (-p_j*(log(p_j)))
      lgx <- log(x)
      lgx[is.infinite(lgx)] <- 0
      trm2 <- sum(x/(length(x))*lgx)  
      return(trm1 + trm2)
    })
    inOut <- sum(inOut)
  }
  return(inOut)
}
################################################################################
# rpg: a new faster, memory efficient function for reading genepop files       #
################################################################################
#' New readgenepop format
#' 
#' Kevin Keenan 2014
rgp <- function(infile){
  fastScan <- function(fname) {
    s <- file.info(fname)$size
    buf <- readChar(fname, s, useBytes = TRUE)
    return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
  }
  if(is.list(infile)){
    infile <- as.matrix(infile)
    dat <- apply(infile, 1, function(x){
      x <- x[!is.na(x)]
      return(paste(x, collapse = "\t"))
    })
    #dat <- c(paste(colnames(infile), collapse = "\t"), dat)
  } else {
    dat <- fastScan(infile)
    # strip whitespace from the end of lines
    dat <- sapply(dat, function(x){
      return(sub("\\s+$", "", x))
    })
    names(dat) <- NULL
  }
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
  genoRet <- function(strt, ends, x){
    out <- strsplit(x[strt:ends], split = "\\s+")
    x <- do.call("rbind", c(out, deparse.level = 0))
    if(round(mean(nchar(x[,2]))) == 1L){
      x[,1] <- paste(x[,1], x[,2], sep = "")
      x <- x[,(-2)]
    }
    x[x == "-9"] <- NA
    x[x == "0000"] <- NA
    x[x == "000000"] <- NA
    # output
    list(ls = x[,(-1)],
         nms = as.vector(x[,1]))
  }
  genos <- mapply(genoRet, strt = strt, ends = ends, 
                  MoreArgs = list(x = dat), SIMPLIFY = FALSE)
  indNames <- lapply(genos, "[[", 2)
  #indNames <- do.call("c", indNames)
  genos <- lapply(genos, "[[", 1)
  # detect genepop format
  gp <- round(mean(nchar(na.omit(genos[[1]][,1:2]))/2))
  # convert genotypes to arrays
  genos <- lapply(genos, function(x){
    al1 <- substr(x, 1, gp)
    al2 <- substr(x, (gp+1), (gp*2))
    out <- array(NA, dim = c(nrow(x), ncol(x), 2))
    out[,,1] <- al1
    out[,,2] <- al2
    return(out)
  })
  
  # calculate allele frequencies, obs alleles, popSizes
  # define function
  statFun <- function(x, cl = NULL){
    # if(!is.null(cl)){
    # tab <- parLapply(cl, 1:dim(x)[2], function(i){return(table(x[,i,]))})
    # } else {
    #tab <- lapply(1:dim(x)[2], function(i){return(table(x[,i,]))})
    #}
    popSizes <- apply(x, 2, function(y){
      length(na.omit(y[,1])) * 2
    })
    af <- lapply(1:dim(x)[2], function(i){
      y <- as.vector(na.omit(x[,i,]))
      nms <- unique(y)[order(unique(y))]
      ot <- myTab(y)
      names(ot) <- nms
      return(ot)
    })
    popSizes <- popSizes/2
    list(af = af, ps = popSizes)
  }
  # rearrange data by loci
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
  # calculate stats
  obsAllSize <- lapply(genos, statFun)
  # get individual stats
  af <- lapply(obsAllSize, function(x){
    out <- x$af
    x$af <- NULL
    return(out)
  })
  #obs <- lapply(obsAllSize, function(x){
  #  return(x$obs)
  #})
  ps <- lapply(obsAllSize, function(x){
    out <- x$ps
    x$ps <- NULL
    return(out)
  })
  af <- lapply(1:(nloci), function(i){
    return(lapply(af, "[[", i))
  })
  #obs <- lapply(1:(nloci), function(i){
  #  return(lapply(obs, "[[", i))
  #})
  ps <- lapply(1:(nloci), function(i){
    return(sapply(ps, "[", i))
  })
  af <- lapply(af, check)
  # names(af) <- locs
  #obs <- lapply(obs, check)
  gc()
  list(af = af, genos = genos, ps = ps, gp = gp,
       indnms = indNames, locs = locs)
}
################################################################################
# end rpg                                                                      #
################################################################################
#' ### StatCalc: Function for the calculation of allele frequencies etc.
#' 
#' __Kevin Keenan__ (2014)
#' 
#' Subfunction for diffCalcRcpp
#' # Calculate pre-diversity/differentiation statistics
#' 
#' This function accepts the output from the `rpg` function and returns
#' some basic statistics which can then be used to calculate Fst etc.
#' 
#' If `statCalc` is passed `idx`, the index rows will be sampled from `plAr`
#' and the statistics describing this re-sample will be returned. Thus, the 
#' function can be used for bootstrapping data.
#' 
# myTab <- function(x){
#   x <- as.character(na.omit(x))
#   mtch <- unique(x)
#   return(sapply(mtch, function(y){
#     return(sum(x == y))
#   }))
# }

# library(Rcpp)
# sourceCpp("src/myTab.cpp")
#' __Kevin Keenan__ (2014)
####-- prestats function --####
statCalc <- function(rsDat, idx = NULL, al, fst, bs = TRUE){
  # generate resamples
  if(bs){
    rsFun <- function(x, y){
      return(x[y,,])
    }
    rsDat <- mapply(rsFun, x = rsDat, y = idx, SIMPLIFY = FALSE) 
  }
  
  # calculate allele frequecies
  
  alf <- lapply(rsDat, function(x){
    apply(x, 2, function(y){
      if(all(is.na(y))){
        return(NA)
      } else {
        y <- as.vector(na.omit(y))
        nms <- unique(y)[order(unique(y))]
        ot <- myTab(y)
        names(ot) <- nms
        return(ot)
      }
    })
  })
  nloci <- length(al)
  
  # organise allele frequencies
  alf <- lapply(1:length(al), function(i){
    lapply(alf, "[[", i)
  })
  
  alSort <- function(x, y){
    idx <- lapply(x, function(z){
      match(names(z), rownames(y))
    })
    for(i in 1:length(idx)){
      y[idx[[i]], i] <- x[[i]]
    }
    return(y)
  }
  
  # generate allele frequency output
  alOut <- mapply(alSort, x = alf, y = al, SIMPLIFY = FALSE)
  # calculate harmonic sizes (only for estimators)
  popSizes <- lapply(rsDat, function(x){
    lgths <- apply(x, 2, function(y){
      nrow(na.omit(y))
    })
    return(lgths)
  })
  ps <- do.call("cbind", popSizes)
  
  # inds per locus
  indtyp <- lapply(1:nloci, function(i){
    vapply(rsDat, function(x){
      op <- length(x[!is.na(x[,i,1]),i,1])
      if(op == 0L){
        return(NA)
      } else {
        return(op) 
      }
    }, FUN.VALUE = numeric(1))
  })
  
  # if wc fst = true calculate
  if(fst){
    hsums <- lapply(rsDat, function(x){
      hts <- lapply(1:dim(x)[2], function(i){
        out <- x[,i,1] == x[,i,2]
        return(out)
      })
      gts <- lapply(1:dim(x)[2], function(i){x[,i,]})
      alls <- lapply(gts, function(y){unique(y[!is.na(y)])})
      #     htcount <- function(gts, hts){
      #      return(table(gts[!hts,]))
      #     }
      htCount <- function(gts, hts, alls){
        if(all(is.na(gts))){
          out <- 0
          names(out) <- "NA"
          return(out)
        } else {
          gts <- gts[!hts, ]
          ht <- sapply(alls, function(al){
            return(sum((gts == al), na.rm = TRUE))
          })
          names(ht) <- alls
          return(ht)
        }
      }
      htcounts <- mapply(htCount, gts = gts, hts = hts, alls = alls, 
                         SIMPLIFY = FALSE)
      return(htcounts)
    })
    # convert hsums to locus focus
    hsums <- lapply(1:nloci, function(i){
      lapply(hsums, "[[", i)
    })
    list(alOut = alOut, ps = ps, hsums = hsums, indtyp = indtyp)
  } else {
    list(alOut = alOut, ps = ps, indtyp = indtyp)
  }
}
##########
# END
##########
#
#
#
#
# preamble ----

#' divMigrate: an experimental function for detecting directional differentiation 
#' A function to calculate pairwise directional differentiation
#' a presented in the paper 'Directional genetic differentiation and
#' asymmetric migration Lisa Sundqvist, Martin Zackrisson & David Kleinhans,
#' 2013, arXiv pre-print (http://arxiv.org/abs/1304.0118)'
#' @export

# function definition ----
divMigrate <- function(infile = NULL, nbs = 0, filter_threshold = 0, 
                       plot = FALSE, plot_col = "darkblue", para = FALSE){
  # preabmle ----
  #nbs <- 1000
  cat("Caution! The method used in this function is still under development. \n")
  #infile <- paste(getwd(), "/test_folder/test_1_1_gen_converted.gen", sep = "")
  
  # read data ----
  #data(Test_data, package = "diveRsity")
  #Test_data[is.na(Test_data)] <- ""
  #Test_data[Test_data == "0"] <- "000000"
  #infile <- Test_data#"test_folder/test_1_1_gen_converted.gen"
  dat <- rgp(infile)
  npops <- length(dat$genos)
  nloci <- length(dat$af)
  
  # generate pw combos ----
  pw <- combn(npops, 2)
  
  # calculate ht and hs ----
  #library(Rcpp) # comment out for package
  #sourceCpp("src/pwHt.cpp") # comment out for package
  hths <- lapply(dat$af, pwHt, pw = pw-1)
  # seperate ht and hs matrices
  ht <- lapply(hths, "[[", "ht")
  hs <- lapply(hths, "[[", "hs")
  
  # Calculate D ----
  # function for locus d
  d <- function(ht, hs){
    return(((ht-hs)/(1-hs))*2)
  }
  # locus d
  dloc <- mapply(`d`, ht = ht, hs = hs, SIMPLIFY = "array")
  # set any nan values to 1
  dloc[is.nan(dloc)] <- 1
  # calculate multilocus d
  hrmD <- apply(dloc, c(1,2), function(x){
    mn <- mean(x, na.rm = TRUE)
    vr <- var(x, na.rm = TRUE)
    return(1/((1/mn) + vr * (1/mn)^3))
  })
  # calculate migration
  dMig <- (1 - hrmD) / hrmD
  # fix infinities
  dMig[is.infinite(dMig)] <- NA
  # calculate relative migration
  dRel <- dMig/max(dMig, na.rm = TRUE)
  dRel[is.nan(dRel)] <- NA
  # plotting network
  dRelPlt <- dRel
  dRelPlt[dRelPlt < filter_threshold] <- 0
  if(nbs != 0L && plot){
    par(mfrow = c(2,1))
  }
  if(plot){
    qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                   legend = TRUE, posCol = plot_col, 
                   edge.labels = TRUE, mar = c(2,2,5,5))
    title(paste("\n Relative migration network (Filter threshold = ", 
                filter_threshold, ")", sep = ""))
    pdf("Relative_migration-[divMigrate].pdf", paper = "a4r")
    qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                   legend = TRUE, posCol = plot_col, 
                   edge.labels = TRUE, mar = c(2,2,5,5))
    title(paste("\n Relative migration network (Filter threshold = ", 
                filter_threshold, ")", sep = ""))
    if(nbs == 0){
      dev.off() 
    } 
  }
  # Bootstrapping ----
  
  if(nbs != 0L){
    # generate bootstrap indexes ----
    ps <- sapply(dat$indnms, length)
    idx <- lapply(1:nbs, function(i){
      lapply(ps, function(x){
        return(sample(x, size = x, replace = TRUE))
      })
    })
    
    # calculate bootstrap D ----
    # load bs function
    #source("R/bsFun.R")
    # run bootstrap function
    if(para){
      library(parallel)
      cl <- makeCluster(detectCores())
      clusterExport(cl, c("bsFun", "dat", "pw"), envir = environment())
      bsD <- parSapply(cl, idx, function(x){
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, pw = pw))
      }, simplify = "array")
      stopCluster(cl)
    } else {
      bsD <- sapply(idx, function(x){
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, pw = pw))
      }, simplify = "array")
    }
    
    # correct bias ----
    #index <- expand.grid(1:length(ps), 1:length(ps))
    #index <- index[!index[,1] == index[,2],]
    #bc <- apply(index, 1, function(i){
    #  bsO <- bsD[i[1], i[2],]
    #  act <- hrmD[i[1], i[2]]
    #  mn <- mean(bsO, na.rm = TRUE) - act
    #  bcN <- bsO - mn
    #  UCI <- quantile(bcN, prob = 0.95)
    #  LCI <- quantile(bcN, prob = 0.05)
    #  out <- c(act, LCI, UCI)
    #  names(out) <- c("ACT", "LCI", "UCI")
    #  return(out)
    #})
    # function for significant difference determination
    sigDiff <- function(x, y){
      if(x[1] < y[1] && x[2] < y[1]){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    # covert bc to 3d array
    #newBc <- array(NA, dim = c(ncol(hrmD), ncol(hrmD), 3))
    #for(i in 1:nrow(index)){
    #  newBc[index[i,1], index[i,2], 1] <- bc[1,i]
    #  newBc[index[i,1], index[i,2], 2] <- bc[2,i]
    #  newBc[index[i,1], index[i,2], 3] <- bc[3,i]
    #}
    # bootstrap means
    #bsMean <- apply(bsD, c(1,2), mean, na.rm = TRUE)
    sigMat <- matrix(NA, nrow = ncol(dRel), ncol(dRel))
    for(i in 1:ncol(pw)){
      p1 <- quantile(bsD[pw[1,i], pw[2,i],], prob = c(0.025, 0.975))
      p2 <- quantile(bsD[pw[2,i], pw[1,i],], prob = c(0.025, 0.975))
      sigMat[pw[2,i], pw[1,i]] <- sigDiff(p1, p2)
      sigMat[pw[1,i], pw[2,i]] <- sigDiff(p2, p1)
    }
    
    # test support for directional diff ----
    #weightMat <- matrix(NA, nrow = ncol(hrmD), ncol = ncol(hrmD))
    #for(i in 1:ncol(pw)){
    #  p1 <- bsD[pw[1,i], pw[2,i],]
    #  p2 <- bsD[pw[2,i], pw[1,i],]
    #  weightMat[pw[2,i], pw[1,i]] <- round(sum(p1 < p2, 
    #                                           na.rm = TRUE)/nbs, 1)
    #  weightMat[pw[1,i], pw[2,i]] <- round(sum(p2 < p1, 
    #                                           na.rm = TRUE)/nbs, 1)
    #}
    #weightMat[bsMean < filter_threshold] <- 0
    #bsMeanPlt <- bsMean
    #bsMeanPlt[bsMean < filter_threshold] <- 0
    #bsMeanSig <- dRel
    dRelPlt[!sigMat] <- 0
    #qgraph::qgraph(bsMeanPlt, nodeNames = sapply(dat$indnms, "[", 1),
    #               legend = TRUE, posCol = "blue", label.color = plot_col,
    #               edge.labels = round(bsMeanPlt, 2))
    #title(paste("Mean relative migration network (", nbs, 
    #            " bootstraps)", sep = ""))
    if(plot){
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, label.color = plot_col,
                     edge.labels = TRUE)
      title(paste("Significant relative migration network (", nbs, 
                  " bootstraps)", sep = ""))
      dev.off()
      #qgraph::qgraph(bsMeanPlt, nodeNames = sapply(dat$indnms, "[", 1),
      #               legend = TRUE, posCol = "blue", label.color = plot_col,
      #               edge.labels = round(bsMeanPlt, 2))
      #title(paste("Mean relative migration network (", nbs, 
      #            " bootstraps)", sep = ""))
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, "[", 1),
                     legend = TRUE, posCol = plot_col, label.color = plot_col,
                     edge.labels = TRUE)
      title(paste("Significant relative migration network (", nbs, 
                  " bootstraps)", sep = ""))
    }
    par(mfrow = c(1,1))
  }
  
  if(nbs != 0L){
    list(relMig = dRel,
         relMigSig = dRelPlt)
  } else {
    list(relMig = dRel)
  }
}

# code end ----
#' Bootstrapping function for use with divMigrate
#' 
#' Kevin Keenan (2014)

# Bootstrapping function definition ----
bsFun <- function(genos, idx, af, pw){
  
  nl <- length(af)
  
  # sub-sample genos ----
  sampleFun <- function(input, idx){
    return(input[idx,,])
  }
  # sub-sample genos
  genos <- mapply(sampleFun, input = genos, idx = idx,
                  SIMPLIFY = FALSE)
  
  # calculate allele frequencies ----
  #sourceCpp("src/myTab.cpp")
  
  alf <- lapply(genos, function(x){
    apply(x, 2, function(y){
      if(all(is.na(y))){
        return(NA)
      } else {
        y <- as.vector(na.omit(y))
        nms <- unique(y)[order(unique(y))]
        ot <- myTab(y)
        names(ot) <- nms
        return(ot)
      }
    })
  })
  
  # organise allele frequencies
  alf <- lapply(1:nl, function(i){
    lapply(alf, "[[", i)
  })
  
  alSort <- function(x, y){
    idx <- lapply(x, function(z){
      match(names(z), rownames(y))
    })
    for(i in 1:length(idx)){
      y[idx[[i]], i] <- x[[i]]
    }
    return(y)
  }
  
  # generate allele frequency output
  af <- mapply(alSort, x = alf, y = af, SIMPLIFY = FALSE)
  
  # calculate hths from boostrapped allele frequencies ----
  hths <- lapply(af, pwHt, pw = pw-1)
  # seperate ht and hs matrices
  ht <- lapply(hths, "[[", "ht")
  hs <- lapply(hths, "[[", "hs")
  
  # Calculate locus Jost's D ----
  
  # function for locus d
  d <- function(ht, hs){
    return(((ht-hs)/(1-hs))*2)
  }
  # locus d
  dloc <- mapply(`d`, ht = ht, hs = hs, SIMPLIFY = "array")
  
  # calculate the harmonic mean of Locus D ----
  
  # set any nan values to 1
  dloc[is.nan(dloc)] <- 1
  # calculate multilocus d
  hrmD <- apply(dloc, c(1,2), function(x){
    mn <- mean(x, na.rm = TRUE)
    vr <- var(x, na.rm = TRUE)
    return(1/((1/mn) + vr * (1/mn)^3))
  })
  
  dMig <- (1 - hrmD) / hrmD
  dMig[is.infinite(dMig)] <- NA
  dRel <- dMig/max(dMig, na.rm = TRUE)
  #rnkMat <- dMig
  #rnkMat[] <- rank(dRel, na.last = "keep")
  return(dRel)
}
# Function definition end ----
################################################################################
#################################     END ALL       ############################
################################################################################