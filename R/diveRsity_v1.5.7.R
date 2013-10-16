################################################################################
################################################################################
##                              diveRsity v1.5.7                              ##  
##                            by Kevin Keenan QUB                             ##  
##            An R package for the calculation of differentiation             ##
##              statistics and locus informativeness statistics               ##  
##                V 1.2.0 and up allows parallel computations                 ##  
##                            GPL3 Kevin Keenan 2013                          ##  
################################################################################
################################################################################

# divPart, a wrapper function for the calculation of differentiation stats.
divPart<-function(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
                  WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE, 
                  bootstraps = 0, plot = FALSE, parallel = FALSE){
                   
  ############################ Argument definitions ############################
  D <- infile
  on <- outfile
  gp <- gp
  fst <- WC_Fst
  bstrps <- bootstraps
  bsls <- bs_locus
  bspw <- bs_pairwise
  plt <- plot
  para <- parallel
  pWise <- pairwise
  
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
    accDat <- pre.divLowMemory(list(infile = D,
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
        write.xlsx(ot1,file=paste(of,"[divPart].xlsx",sep=""),
                   sheetName="Standard_stats",col.names=T,
                   row.names=F,append=F)
        # Estimated stats
        write.xlsx(ot2,file=paste(of,"[divPart].xlsx",sep=""),
                   sheetName="Estimated_stats",col.names=T,
                   row.names=F,append=T)
      } else {
        # text file alternatives
        std<-file(paste(of,"Standard-stats[divPart].txt",sep=""), "w")
        cat(paste(colnames(ot1),sep=""),"\n",sep="\t",file=std)
        for(i in 1:nrow(ot1)){
          cat(ot1[i,],"\n",file=std,sep="\t")
        }
        close(std)
        est<-file(paste(of,"Estimated-stats[divPart].txt",sep=""),"w")
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
    if (para && para_pack) {
      #count cores
      library("doParallel")
      cores <- detectCores()
      cl<-makeCluster(cores)
      registerDoParallel(cl)
    }
    
    # Used only if bootstraps is greater than zero
    if(bsls == TRUE){
      
      if (para && para_pack) {
        
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
          write.xlsx(bs_out,file=paste(of,"[divPart].xlsx",sep=""),
                     sheetName="Locus_bootstrap",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          bts<-file(paste(of,"Locus-bootstrap[divPart].txt",sep=""), "w")
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
    if(pWise || bspw){
      pw <- combn(accDat$npops,2)
      pwmat <- pw + 1
      #pw data creator
      ind_vectors <- lapply(1:accDat$npops, function(x){
        rep(x, accDat$pop_sizes[[x]])}
      )
      #      
      pre_data <- matrix(rep("", ((accDat$nloci + 1) * (accDat$nloci + 1))),
                         ncol = (accDat$nloci + 1))
      pre_data[1,] <- rep("", (accDat$nloci + 1))
      #
      for(i in 2:(accDat$nloci + 1)){
        pre_data[i, 1] <- accDat$locus_names[(i-1)]
      }
      #
      pw_data<-list()
      for (i in 1:ncol(pw)){
        pw_data[[i]]<-data.frame(rbind(pre_data,
                                       c("POP",as.vector(rep("",accDat$nloci))),
                                       cbind(ind_vectors[[pw[1,i]]],
                                             matrix(noquote(accDat$pop_list
                                                            [[pw[1,i]]]),
                                                    ncol=accDat$nloci)),
                                       c("POP",as.vector(rep("",accDat$nloci))),
                                       cbind(ind_vectors[[pw[2,i]]],
                                             matrix(noquote(accDat$pop_list
                                                            [[pw[2,i]]]),
                                                    ncol=accDat$nloci))))
      }
      true_stat_gp_in <- list()
      if(fst == TRUE){
        pw_glb <- matrix(rep(0, (8 * (ncol(pw)))), ncol = 8)
      } else {
        pw_glb <- matrix(rep(0, (6 * (ncol(pw)))), ncol = 6)
      }
      for (i in 1:ncol(pw)){
        true_stat_gp_in[[i]] <- list(infile = pw_data[[i]],
                                     gp = gp, bootstrap = FALSE,
                                     locs = FALSE, fst = fst)
      }
      if (para && para_pack) {
        
        true_stat <- parLapply(cl, true_stat_gp_in, pre.divLowMemory)
        # close core connections if not needed further
        if (bspw == FALSE){
          stopCluster(cl)
        }
      } else {
        true_stat <- lapply(true_stat_gp_in, pre.divLowMemory)
      }
      for(i in 1:ncol(pw)){
        if(fst==TRUE){
          pw_glb[i,]<-c(true_stat[[i]]$gst_all,true_stat[[i]]$gst_all_hedrick,
                        true_stat[[i]]$djost_all,true_stat[[i]]$gst_est_all,
                        true_stat[[i]]$gst_est_all_hedrick,
                        true_stat[[i]]$djost_est_all,
                        as.numeric(true_stat[[i]]$fstat[2:3]))
        } else {
          pw_glb[i,]<-c(true_stat[[i]]$gst_all,true_stat[[i]]$gst_all_hedrick,
                        true_stat[[i]]$djost_all,true_stat[[i]]$gst_est_all,
                        true_stat[[i]]$gst_est_all_hedrick,
                        true_stat[[i]]$djost_est_all)
        }
      }
      if(fst==TRUE){
        pwMatList <- lapply(1:8, function(x){
          matrix(rep("--", ((accDat$npops+1) ^ 2)), 
                 ncol = (accDat$npops + 1),
                 nrow = (accDat$npops + 1))
        })
      } else {
        pwMatList <- lapply(1:6, function(x){
          matrix(rep("--", ((accDat$npops+1)^2)),
                 ncol = (accDat$npops + 1),
                 nrow = (accDat$npops + 1))
        })
      }
      if(fst==TRUE){
        pwMatListOut <- lapply(1:8, function(x){
          matrix(rep(NA, ((accDat$npops)^2)),
                 ncol = (accDat$npops),
                 nrow = (accDat$npops))
        })
      } else {
        pwMatListOut <- lapply(1:6, function(x){
          matrix(rep(NA,((accDat$npops)^2)),
                 ncol = (accDat$npops),
                 nrow = (accDat$npops))
        })
      }
      names(pwMatList) <- namer
      names(pwMatListOut) <- namer
      #write pw res to matrices
      pnames <- c("", accDat$pop_names)
      pnamesOut <- accDat$pop_names
      if(fst==TRUE){
        for(i in 1:8){
          for(j in 1:ncol(pw)){
            pwMatList[[i]][pwmat[2, j], pwmat[1, j]] <- pw_glb[j, i]
            pwMatList[[i]][pwmat[1, j], pwmat[2, j]] <- ""
            pwMatListOut[[i]][pw[2, j], pw[1, j]] <- pw_glb[j, i]
            #pwMatListOut[[i]][pw[1,j],pw[2,j]]<-""
          }
          pwMatList[[i]][1, ] <- pnames
          pwMatList[[i]][, 1] <- pnames
          dimnames(pwMatListOut[[i]]) <- list(pnamesOut, pnamesOut)
        }
      } else {
        for(i in 1:6){
          for(j in 1:ncol(pw)){
            pwMatList[[i]][pwmat[2, j], pwmat[1, j]] <- pw_glb[j, i]
            pwMatList[[i]][pwmat[1, j], pwmat[2, j]] <- ""
            pwMatListOut[[i]][pw[2, j], pw[1, j]] <- pw_glb[j, i]
            #pwMatListOut[[i]][pw[1,j],pw[2,j]]<-""
          }
          pwMatList[[i]][1, ] <- pnames
          pwMatList[[i]][, 1] <- pnames
          dimnames(pwMatListOut[[i]]) <- list(pnamesOut, pnamesOut)
        }
      }
      
      
      # write object create
      #pnames list
      
      pwWrite <- pwMatList[[1]]
      pwWrite <- rbind(c(names(pwMatList)[1], rep("", accDat$npops)), pwWrite,
                       rep("", (accDat$npops + 1)))
      if(fst==TRUE){
        for(i in 2:8){
          pwWrite <- rbind(pwWrite, c(names(pwMatList)[i],
                                      rep("", accDat$npops)),
                           pwMatList[[i]], rep("", (accDat$npops + 1)))
        }
      } else {
        for(i in 2:6){
          pwWrite <- rbind(pwWrite, c(names(pwMatList)[i],
                                      rep("",accDat$npops)),
                           pwMatList[[i]], rep("",(accDat$npops+1)))
        }
      }
      if(!is.null(on)){
        if(write_res == TRUE){
          # write data to excel
          # Load dependencies
          # pw stats
          write.xlsx(pwWrite,file=paste(of,"[divPart].xlsx",sep=""),
                     sheetName="Pairwise-stats",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          pw_outer<-file(paste(of,"Pairwise-stats[divPart].txt",sep=""), "w")
          for(i in 1:nrow(pwWrite)){
            cat(pwWrite[i,],"\n",file=pw_outer,sep="\t")
          }
          close(std)
        }
      }
      #cleanup
      rm("pwWrite")
      ##
      zzz<-gc()
      rm(zzz)
    }
    
    
    #Bootstrap
    if(bspw == TRUE){
      
      
      # Bootstrap results data object 
      # bs_pw_glb = bootstrap pairwise global stats
      #if(fst == TRUE){
      #  bs_pw_glb <- matrix(rep(0, (8*bstrps)), ncol = 8, nrow = bstrps)
      #} else {
      #  bs_pw_glb <- matrix(rep(0, (6*bstrps)), ncol = 6, nrow = bstrps)
      #}
      # output results data object
      # pw_res = pairwise results
      if(fst==TRUE){
        pw_res <- lapply(1:8, function(x){
          matrix(nrow = ncol(pw), ncol = 3)
        })
      } else {
        pw_res <- lapply(1:6, function(x){
          matrix(nrow = ncol(pw), ncol = 3)
        })
      }
      #
      #
      
      #parallel processing option
      if (para && para_pack) {
        #create a readGenepopX list
        bs_pw_glb<-list()
        data_res<-list()
        bs_pw_para<-list()
        for(i in 1:ncol(pw)){
          input <- list(infile = pw_data[[i]], gp = gp, bootstrap = TRUE,
                        locs = FALSE, fst = fst)
          # silence for memory efficiency
          #pw_inlist<-list()
          #for(j in 1:bstrps){
          #  pw_inlist[[j]] <- input
          #}
          if(fst == TRUE){
            bs_pw_glb[[i]] <- matrix(rep(0, (8*bstrps)), ncol = 8,
                                     nrow = bstrps)
          } else {
            bs_pw_glb[[i]] <- matrix(rep(0, (6*bstrps)), ncol = 6, 
                                     nrow = bstrps)
          }
          clusterExport(cl, c("input", "pre.divLowMemory"),
                        envir = environment())
          bs_pw_para <- parLapply(cl, 1:bstrps, function(...){
            pre.divLowMemory(input)
          })
          for(j in 1:bstrps){
            if(fst == TRUE){
              bs_pw_glb[[i]][j,] <- c(bs_pw_para[[j]]$gst_all,
                                      bs_pw_para[[j]]$gst_all_hedrick,
                                      bs_pw_para[[j]]$djost_all,
                                      bs_pw_para[[j]]$gst_est_all,
                                      bs_pw_para[[j]]$gst_est_all_hedrick,
                                      bs_pw_para[[j]]$djost_est_all,
                                      as.numeric(bs_pw_para[[j]]$fstats[2:3]))
            } else {
              bs_pw_glb[[i]][j,] <- c(bs_pw_para[[j]]$gst_all,
                                      bs_pw_para[[j]]$gst_all_hedrick,
                                      bs_pw_para[[j]]$djost_all,
                                      bs_pw_para[[j]]$gst_est_all,
                                      bs_pw_para[[j]]$gst_est_all_hedrick,
                                      bs_pw_para[[j]]$djost_est_all)
              
            }
          }
        }
        #
        # confidence interval calculator function
        pwCi <- lapply(bs_pw_glb, function(x){
          res <- apply(x, 2, function(y){
            ci <- as.vector(quantile(y, probs = c(0.025, 0.975), na.rm = TRUE))
            means <- mean(y, na.rm = TRUE)
            return(c(means, ci))
          })
          mu <- res[1,]
          lci <- res[2,]
          uci <- res[3,]
          list(mu = mu, lci = lci, uci = uci)
        })
        # create easy access data structure for each
        mu <- t(sapply(1:length(pwCi), function(i){
          return(pwCi[[i]]$mu)
        }))
        lci <- t(sapply(1:length(pwCi), function(i){
          return(pwCi[[i]]$lci)
        }))
        uci <- t(sapply(1:length(pwCi), function(i){
          return(pwCi[[i]]$uci)
        }))
        
        for(i in 1:ncol(pw)){
          for(j in 1:ncol(mu)){
            pw_res[[j]][i, 1] <- round(mu[i, j], 4)
            pw_res[[j]][i, 2] <- round(lci[i, j], 4)
            pw_res[[j]][i, 3] <- round(uci[i, j], 4)
            pw_res[[j]][is.na(pw_res[[j]])] <- 0
          }
        }
        stopCluster(cl)
      } else {
        #sequential vectorized
        #pw_inlist<-list()
        #for(i in 1:ncol(pw)){
        #  input <- list(infile = pw_data[[i]],
        #                gp = gp, bootstrap = TRUE, 
        #                locs = FALSE, fst = fst)
        #  pw_inlist[[i]] <- list()
        #  for(j in 1:bstrps){
        #    pw_inlist[[i]][[j]] <- input
        #  }
        #}
        bs_pw_glb <- list()
        for(i in 1:ncol(pw)){
          if(fst == TRUE){
            bs_pw_glb[[i]] <- matrix(rep(0, (8*bstrps)), ncol = 8,
                                     nrow = bstrps)
          } else {
            bs_pw_glb[[i]] <- matrix(rep(0, (6*bstrps)), ncol = 6, 
                                     nrow = bstrps)
          }
        }
        #create a readGenepopX list
        bs_pw_glb <- list()
        data_res <- list()
        bs_pw_para <- list()
        for(i in 1:ncol(pw)){
          input <- list(infile = pw_data[[i]],
                        gp = gp, bootstrap = TRUE,
                        locs = FALSE, fst = fst)
          # silence for memory efficiency
          #pw_inlist <- list()
          #for(j in 1:bstrps){
          #  pw_inlist[[j]] <- input
          #}
          if(fst == TRUE){
            bs_pw_glb[[i]] <- matrix(rep(0, (8*bstrps)), ncol = 8,
                                     nrow = bstrps)
          } else {
            bs_pw_glb[[i]] <- matrix(rep(0, (6*bstrps)), ncol = 6,
                                     nrow = bstrps)
          }
          bs_pw_para <- lapply(1:bstrps, function(...){
            pre.divLowMemory(input)
          })
          for(j in 1:bstrps){
            if(fst == TRUE){
              bs_pw_glb[[i]][j,] <- c(bs_pw_para[[j]]$gst_all,
                                      bs_pw_para[[j]]$gst_all_hedrick,
                                      bs_pw_para[[j]]$djost_all,
                                      bs_pw_para[[j]]$gst_est_all,
                                      bs_pw_para[[j]]$gst_est_all_hedrick,
                                      bs_pw_para[[j]]$djost_est_all,
                                      as.numeric(bs_pw_para[[j]]$fstat[2:3]))
              
            } else {
              bs_pw_glb[[i]][j,] <- c(bs_pw_para[[j]]$gst_all,
                                      bs_pw_para[[j]]$gst_all_hedrick,
                                      bs_pw_para[[j]]$djost_all,
                                      bs_pw_para[[j]]$gst_est_all,
                                      bs_pw_para[[j]]$gst_est_all_hedrick,
                                      bs_pw_para[[j]]$djost_est_all)
              
            }
          }
        } 
        # confidence interval calculator function
        pwCi <- lapply(bs_pw_glb, function(x){
          res <- apply(x, 2, function(y){
            ci <- as.vector(quantile(y, probs = c(0.025, 0.975), na.rm = TRUE))
            means <- mean(y, na.rm = TRUE)
            return(c(means, ci))
          })
          mu <- res[1,]
          lci <- res[2,]
          uci <- res[3,]
          list(mu = mu, lci = lci, uci = uci)
        })
        # create easy access data structure for each
        mu <- t(sapply(1:length(pwCi), function(i){
          return(pwCi[[i]]$mu)
        }))
        lci <- t(sapply(1:length(pwCi), function(i){
          return(pwCi[[i]]$lci)
        }))
        uci <- t(sapply(1:length(pwCi), function(i){
          return(pwCi[[i]]$uci)
        }))
        
        for(i in 1:ncol(pw)){
          for(j in 1:ncol(mu)){
            pw_res[[j]][i, 1] <- round(mu[i, j], 4)
            pw_res[[j]][i, 2] <- round(lci[i, j], 4)
            pw_res[[j]][i, 3] <- round(uci[i, j], 4)
            pw_res[[j]][is.na(pw_res[[j]])] <- 0
          }
        }
        #
      }
      #
      # pairwise comparisons
      # pw_names = pairwise population names
      pw_nms <- paste(accDat$pop_names[pw[1,]],
                      accDat$pop_names[pw[2,]], sep = " vs. ")
      #
      pw_nms1 <- paste(pw[1,], pw[2,], sep = " vs. ")
      #
      names(pw_res) <- namer
      #
      pw_res1 <- pw_res
      if(fst == TRUE){
        for(i in 1:8){
          dimnames(pw_res1[[i]]) <- list(pw_nms, 
                                         c("Mean", "Lower_CI", "Upper_CI"))
        }
      } else {
        for(i in 1:6){
          dimnames(pw_res1[[i]]) <- list(pw_nms, 
                                         c("Mean", "Lower_CI", "Upper_CI"))
        }
      }
      # bs results output object header
      hdr <- matrix(c("Pairwise", "Mean", "Lower_95%CI", "Upper_95%CI"),
                    ncol = 4)
      pw_bs_out <- matrix(rbind(hdr, c(names(pw_res)[1],"" ,"" ,""),
                                cbind(pw_nms, pw_res[[1]])), ncol = 4)
      if(fst == TRUE){
        for(i in 2:8){
          pw_bs_out <- matrix(rbind(pw_bs_out, c(names(pw_res)[i], "", "", ""),
                                    cbind(pw_nms, pw_res[[i]])), ncol = 4)
        }
      } else {
        for(i in 2:6){
          pw_bs_out <- matrix(rbind(pw_bs_out, c(names(pw_res)[i], "", "", ""),
                                    cbind(pw_nms, pw_res[[i]])), ncol = 4)
        }
      }
      if(!is.null(on)){
        if(write_res==TRUE){
          write.xlsx(pw_bs_out,file=paste(of,"[divPart].xlsx",sep=""),
                     sheetName="Pairwise_bootstrap",col.names=F,
                     row.names=F,append=T)
        } else {
          # text file alternatives
          pw_bts<-file(paste(of,"Pairwise-bootstrap[divPart].txt",sep=""), "w")
          cat(paste(colnames(pw_bs_out),sep=""),"\n",sep="\t",file=pw_bts)
          for(i in 1:nrow(pw_bs_out)){
            cat(pw_bs_out[i,],"\n",file=pw_bts,sep="\t")
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
      pwso<-list()
      for(i in 1:length(pw_res)){
        pwso[[i]]<-order(pw_res[[i]][,1],decreasing=F)
        #if(length(pwso[[i]]) >= 100){
        #  pwso[[i]]<-pwso[[i]][(length(pwso[[i]])-99):length(pwso[[i]])]
        #}
      }
      names(pwso)<-namer
      # define plot parameters 
      plot.call_pw<-list()
      plot.extras_pw<-list()
      xy.labels_pw<-list()
      y.pos_pw<-list()
      x.pos_pw=1:length(pwso[[i]])
      fn_pre_pw<-list()
      direct=of
      #Plot Gst_Nei
      plot.call_pw[[1]]=c("plot(pw_res[[4]][pwso[[4]],1],
                          ylim=c(0,(max(pw_res[[4]][,3])+
                          min(pw_res[[4]][,3]))),xaxt='n',
                          ylab=names(pw_res)[4],type='n',
                          xlab='Pairwise comparisons 
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[1]]=c("points(pw_res[[4]][pwso[[4]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],2],
                            1:length(pwso[[4]]),pw_res[[4]][pwso[[4]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[5]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[1]]=data.frame(pairwise_name=pw_nms[pwso[[4]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[4]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[4]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[4]],1],4))
      
      y.pos_pw[[1]]=pw_res[[4]][pwso[[4]],1]
      fn_pre_pw[[1]]<-names(pw_res)[4]
      
      
      
      # Plot Gst_Hedrick
      plot.call_pw[[2]]=c("plot(pw_res[[5]][pwso[[5]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(pw_res)[5],type='n',
                          xlab='Pairwise comparisons
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[2]]=c("points(pw_res[[5]][pwso[[5]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[5]]),pw_res[[5]][pwso[[5]],2],
                            1:length(pwso[[5]]),pw_res[[5]][pwso[[5]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[6]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[2]]=data.frame(pairwise_name=pw_nms[pwso[[5]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[5]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[5]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[5]],1],4))
      
      y.pos_pw[[2]]=pw_res[[5]][pwso[[5]],1]
      fn_pre_pw[[2]]<-names(pw_res)[5]
      
      
      # Plot D_jost
      plot.call_pw[[3]]=c("plot(pw_res[[6]][pwso[[6]],1],
                          ylim=c(0,1),xaxt='n',ylab=names(pw_res)[6],type='n',
                          xlab='Pairwise comparisons 
                          \n (Hover over a point to see pairwise info.)',
                          cex.lab=1.2,cex.axis=1.3,las=1)")
      
      plot.extras_pw[[3]]=c("points(pw_res[[6]][pwso[[6]],1],
                            pch=15,col='black',cex=1);
                            arrows(1:length(pwso[[6]]),pw_res[[6]][pwso[[6]],2],
                            1:length(pwso[[6]]),pw_res[[6]][pwso[[6]],3],code=3,
                            angle=90,length=0.05,lwd=0.1);
                            abline(h=as.numeric(plot_data321[7]),
                            lwd=1,lty=2,col='red')")
      
      xy.labels_pw[[3]]=data.frame(pairwise_name=pw_nms[pwso[[6]]],
                                   Gst_Nei=round(pw_res[[4]][pwso[[6]],1],4),
                                   Gst_Hedrick=round(pw_res[[5]][pwso[[6]],1],4),
                                   D_jost=round(pw_res[[6]][pwso[[6]],1],4))
      
      y.pos_pw[[3]]=pw_res[[6]][pwso[[6]],1]
      fn_pre_pw[[3]]<-names(pw_res)[6]
      #plot(Fst_WC)
      if(fst==TRUE){
        plot.call_pw[[4]]=c("plot(pw_res[[8]][pwso[[8]],1],
                            ylim=c(0,(max(pw_res[[8]][,3])+
                            min(pw_res[[8]][,3]))),xaxt='n',ylab=names(pw_res)[8],type='n',
                            xlab='Pairwise comparisons 
                            \n (Hover over a point to see pairwise info.)',
                            cex.lab=1.2,cex.axis=1.3,las=1)")
        
        plot.extras_pw[[4]]=c("points(pw_res[[8]][pwso[[8]],1],
                              pch=15,col='black',cex=1);
                              arrows(1:length(pwso[[8]]),pw_res[[8]][pwso[[8]],2],
                              1:length(pwso[[8]]),pw_res[[8]][pwso[[8]],3],code=3,
                              angle=90,length=0.05,lwd=0.1);
                              abline(h=as.numeric(plot_data321[7]),
                              lwd=1,lty=2,col='red')")
        
        xy.labels_pw[[4]]=data.frame(pairwise_name=pw_nms[pwso[[8]]],
                                     Gst_Nei=round(pw_res[[4]][pwso[[8]],1],4),
                                     Gst_Hedrick=round(pw_res[[5]][pwso[[8]],1],4),
                                     D_jost=round(pw_res[[6]][pwso[[8]],1],4),
                                     Fst_WC=round(pw_res[[8]][pwso[[8]],1],4))
        
        y.pos_pw[[4]]=pw_res[[8]][pwso[[8]],1]
        fn_pre_pw[[4]]<-names(pw_res)[8]
      }
    }
    ############################### Bootstrap end ################################
    
    
    ################################# Plot resuts ################################
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
           bs_pairwise = pw_res1)
    } else if(bspw == TRUE && bsls == FALSE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_pairwise = pw_res1)
    } else if(bspw == FALSE && bsls == TRUE && pWise == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise,
           bs_locus = bs_res1)
    } else if(bspw == FALSE && bsls == FALSE && pWise == TRUE){
      list(standard = ot1out,
           estimate = ot2out,
           pairwise = pwMatListOut,
           meanPairwise = meanPairwise)
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
# divPart end                                                                  #
################################################################################
#
#
#
#
#
################################################################################
# div.part: deprecated
################################################################################
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
  gp=x$gp
  infile=x$infile
  bootstrap=x$bootstrap
  locs=x$locs
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
  npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                  which(data1[,1]=="pop")))
  pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
              which(data1[,1]=="pop"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
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
# inCalc, a wrapper function for the calculation of locus informativeness     #
################################################################################
inCalc<-function(infile, outfile=NULL, gp=3, bs_locus=FALSE, bs_pairwise=FALSE,
                 bootstraps=0, plot=FALSE, parallel=FALSE){
  D=infile
  gp=gp
  pw=bs_pairwise
  BS=bs_locus
  NBS=bootstraps
  on=outfile
  plt=plot
  para = parallel
  if(pw==T && NBS<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else if (BS==T && NBS<2){
    bs_warning<-{paste("[STOPPED]",
                       "bootsraps must be greater than 2")
    }
    cat(noquote(bs_warning))
  } else {
    write_res<-is.element("xlsx",installed.packages()[,1])
    if(write_res==TRUE && !is.null(on)){
      require("xlsx")
    } else {
      if(!is.null(on)){
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
    }
    if(!is.null(on)){
      suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                             "-[diveRsity]","/",sep="")))
    }
    
    of <- paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
    # Parallel system opti
    if(para){
      para_pack_inst<-is.element(c("parallel","doParallel","foreach",
                                   "iterators"),installed.packages()[,1])
      para_pack <- all(para_pack_inst)
    }
    if (para && para_pack == FALSE){
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
    ##
    
    #source("in.bootstrap.R")
    inls2<-list(D,gp,"FALSE",0,"FALSE")
    res_out<-in.bs(inls2)[[1]]
    if(!is.null(on)){
      if(write_res==TRUE){
        write.xlsx(res_out,file=paste(of,"[inCalc].xlsx",sep=""),
                   sheetName="In_allele_stats",col.names=T,row.names=T,append=F)
      } else {
        all_out<-file(paste(of,"Allele-In[inCalc].txt",sep=""),"w")
        cat(paste(colnames(res_out),sep=""),"\n",sep="\t",file=all_out)
        for(i in 1:nrow(res_out)){
          cat(res_out[i,],"\n",sep="\t",file=all_out)
        }
        close(all_out)
      }
    }
    ######################################################################
    # overall In
    if(BS==T){
      inls1<-list(D,gp,BS,NBS,"TRUE")
      bs_sum1<-in.bs(inls1)
      if(!is.null(on)){
        if(write_res==T){
          write.xlsx(bs_sum1,file=paste(of, "[inCalc].xlsx",sep=""),
                     sheetName="Overall_Bootstrap",col.names=T,
                     row.names=T,append=T)
        } else {
          all_bs<-file(paste(of,"Overall-bootstrap[inCalc].txt",sep=""),"w")
          cat(paste(colnames(bs_sum1),sep=""),"\n",sep="\t",file=all_bs)
          for(i in 1:nrow(bs_sum1)){
            cat(bs_sum1[i,],"\n",sep="\t",file=all_bs)
          }
          close(all_bs)
        }
      }     
      loc_nms<-rownames(bs_sum1)
      if(plt && !is.null(on)){
        lso<-order(bs_sum1[,1],decreasing=F)
        png(filename=paste(of, on,"_In_plot.png",sep=""),width=800,height=600)
        par(mar=c(6,5,1,1))
        plot(bs_sum1[lso,1],ylim=c(0,(max(bs_sum1[,3])+0.1)),xaxt='n',
             ylab=expression('Locus '*I[n]),
             xlab="",cex.lab=1.5,cex.axis=1.3,las=1,type='n')
        points(bs_sum1[lso,1],pch=15,col='black',cex=1)
        suppressWarnings(arrows(1:nrow(bs_sum1),bs_sum1[lso,2], 
                                1:nrow(bs_sum1),bs_sum1[lso,3],
                                code=3,angle=90,length=0.05,
                                lwd=0.1))
        axis(1,at=1:nrow(bs_sum1),labels=loc_nms[lso],las=3)
        dev.off()
      }
    }
    # pairwise locus In bootstrap
    if(pw==T){
      inls<-list(D, gp, FALSE, TRUE)
      names(inls)<-c("infile","gp","bootstrap","locs")
      data<-readGenepopX(inls)
      af<-data$allele_freq
      np<-data$npops
      nl<-data$nloci
      nal<-data$nalleles
      ln<-data$loci_names
      ps<-data$pop_sizes
      pl<-data$pop_list
      pwc<-combn(np,2)
      pn<-data$pop_names
      
      iv<-list()
      for(i in 1:np){
        iv[[i]]<-noquote(paste(rep(i,ps[i]),",",sep=""))
      }
      
      
      pre_data<-matrix(rep("",((nl+1)*(nl+1))),
                       ncol=(nl+1))
      pre_data[1,]<-rep("",(nl+1))
      for(i in 2:(nl+1)){
        pre_data[i,1]<-ln[(i-1)]
      }
      
      pw_data<-list()
      for (i in 1:ncol(pwc)){
        pw_data[[i]]<-data.frame(rbind(pre_data,
                                       c("POP",as.vector(rep("",nl))),
                                       cbind(iv[[pwc[1,i]]],
                                             matrix(noquote(pl[[pwc[1,i]]]),
                                                    ncol=nl)),
                                       c("POP",as.vector(rep("",nl))),
                                       cbind(iv[[pwc[2,i]]],
                                             matrix(noquote(pl[[pwc[2,i]]]),
                                                    
                                                    ncol=nl))))
      }
      pw_bs<-list()
      pw_bs_in<-list()
      pw_only<-TRUE
      pw_bs_out<-list()
      for(i in 1:ncol(pwc)){
        pw_bs_in[[i]]<-list(pw_data[[i]],gp,pw,NBS,pw_only)
      }
      if (para && para_pack){
        library("doParallel")
        cores <- detectCores()
        cl <- makeCluster(cores)
        registerDoParallel(cl)
        pw_bs<-parLapply(cl, pw_bs_in, in.bs)
        stopCluster(cl)
      } else {
        pw_bs<-lapply(pw_bs_in, in.bs)
      }
      for(i in 1:ncol(pwc)){
        #  pw_bs[[i]]<-in.bs(pw_data[[i]],gp,pw,NBS)[[2]]
        pw_bs_out[[i]]<-matrix(cbind(rownames(pw_bs[[i]]),
                                     pw_bs[[i]][,1:3]),ncol=4)
      }
      pw_nms<-paste(pn[pwc[1,]],pn[pwc[2,]],sep=" vs. ")
      names(pw_bs)<-pw_nms
      hdr<-c("Loci","Actual_In","Lower_95CI","Upper_95CI")
      pw_in_bs<-matrix(rbind(hdr,c(names(pw_bs)[1],"","",""),pw_bs_out[[1]]),
                       ncol=4)
      for(j in 2:ncol(pwc)){
        pw_in_bs<-matrix(rbind(pw_in_bs,c(names(pw_bs)[j],"","",""),
                               pw_bs_out[[j]]),ncol=4)
      }
      if(!is.null(on)){
        if(write_res==TRUE){
          write.xlsx(pw_in_bs,file=paste(of, "[inCalc].xlsx",sep=""),
                     sheetName="Pairwise_bootstraps",col.names=F,
                     row.names=F,append=T)
        } else {
          pw_bs<-file(paste(of,"Pairwise-bootstrap[inCalc].txt",sep=""),"w")
          cat(paste(colnames(pw_in_bs),sep=""),"\n",sep="\t",file=pw_bs)
          for(i in 1:nrow(pw_in_bs)){
            cat(pw_in_bs[i,],"\n",sep="\t",file=pw_bs)
          }
          close(pw_bs)
        }
      }
    }
    
    if(BS==F && pw==F){
      list(Allele_In=res_out)
    } else if (BS==T && pw==F){
      list(Allele_In=res_out,
           l_bootstrap=bs_sum1)
    } else if (BS==F && pw==T){
      list(Allele_In=res_out,
           PW_bootstrap=pw_bs)
    } else if (BS==T && pw==T){
      list(Allele_In=res_out,
           l_bootstrap=bs_sum1,
           PW_bootstrap=pw_bs)
    }      
  }
}
################################################################################
# inCalc end                                                                  #
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
#
#
#
#
################################################################################
# in.bs, a function for the bootstrap calculations of locus informativeness    #
################################################################################
in.bs<-function(x){
  D=x[[1]]
  gp=x[[2]]
  BS=x[[3]]
  NBS=x[[4]]
  pw_only=x[[5]]
  readGenepopX <- function (x) {
    gp=x$gp
    infile=x$infile
    bootstrap=x$bootstrap
    locs=x$locs
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
    npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                    which(data1[,1]=="pop")))
    pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
                which(data1[,1]=="pop"),(nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
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
    for(i in 1:npops){
      for(j in 1:nloci){
        allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
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
           bs_file=bs_data_file)
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
           locs=locs)
    }
  }
  inls<-list(D, gp, FALSE, TRUE)
  names(inls)<-c("infile","gp","bootstrap","locs")
  data<-readGenepopX(inls)
  af<-data$allele_freq
  np<-data$npops
  nl<-data$nloci
  nal<-data$nalleles
  ln<-data$loci_names
  ps<-data$pop_sizes
  pl<-data$pop_list
  ## Calc P[i]
  p<-list()
  for (i in 1:nl){
    p[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      p[[i]][j]<- sum(af[[i]][j,])/np
    }
  }
  exp1<-list()
  for (i in 1:nl){
    exp1[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      exp1[[i]][j]<- (-p[[i]][j]*log(p[[i]][j]))
    }
  }
  exp2_sep<-list()
  for(i in 1:nl){
    exp2_sep[[i]]<-matrix(rep(0,np*nal[i]))
    dim(exp2_sep[[i]])<-c(nal[i],np)
  }
  for(i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      for (z in 1:np){
        exp2_sep[[i]][j,z]<-(af[[i]][j,z]/np)*
          log(af[[i]][j,z])
      }
    }
  }
  ## Replace NaN's with 0.000
  for (i in 1:nl){
    exp2_sep[[i]][exp2_sep[[i]]=="NaN"]<-0
  }
  exp2<-list()
  for (i in 1:nl){
    exp2[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      exp2[[i]][j]<-sum(exp2_sep[[i]][j,])
    }
  }
  In<-list()
  for (i in 1:nl){
    In[[i]]<-vector()
  }
  for (i in 1:nl){
    for (j in 1:nrow(af[[i]])){
      In[[i]][j]<-sum(exp1[[i]][j]+exp2[[i]][j])
    }
  }
  In_sum<-vector()
  for (i in 1:nl){
    In_sum[i]<-sum(In[[i]])
  }
  results_out<-matrix(rep(NA,(nl*(max(nal)+1))),
                      nrow=nl,ncol=(max(nal)+1))
  for (i in 1:nl){
    results_out[i,1:nal[i]]<- round(In[[i]],4)
  }
  results_out[,(max(nal)+1)]<-round(In_sum,4)
  ##############################################################################
  if(BS==T){
    bs_sum<-matrix(rep(0,(nl*NBS)),ncol=nl)
    colnames(bs_sum)<-ln
    inls_bs<-list(D,gp,TRUE,TRUE)
    names(inls_bs)<-c("infile","gp","bootstrap","locs")
    for(w in 1:NBS){
      data_bs<-readGenepopX(inls_bs)
      af_bs<-data_bs$allele_freq
      np_bs<-data_bs$npops
      nl_bs<-data_bs$nloci
      nal_bs<-data_bs$nalleles
      ln_bs<-data_bs$loci_names
      ## Calc P[i]
      p_bs<-list()
      for (i in 1:nl_bs){
        p_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          p_bs[[i]][j]<- sum(af_bs[[i]][j,])/np_bs
        }
      }
      exp1_bs<-list()
      for (i in 1:nl_bs){
        exp1_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          exp1_bs[[i]][j]<- (-p_bs[[i]][j]*log(p_bs[[i]][j]))
        }
      }
      exp2_sep_bs<-list()
      for(i in 1:nl_bs){
        exp2_sep_bs[[i]]<-matrix(rep(0,np_bs*nal_bs[i]))
        dim(exp2_sep_bs[[i]])<-c(nal_bs[i],np_bs)
      }
      for(i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          for (z in 1:np_bs){
            exp2_sep_bs[[i]][j,z]<-(af_bs[[i]][j,z]/np_bs)*
              log(af_bs[[i]][j,z])
          }
        }
      }
      ## Replace NaN's with 0.000
      for (i in 1:nl_bs){
        exp2_sep_bs[[i]][exp2_sep_bs[[i]]=="NaN"]<-0
      }
      exp2_bs<-list()
      for (i in 1:nl_bs){
        exp2_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          exp2_bs[[i]][j]<-sum(exp2_sep_bs[[i]][j,])
        }
      }
      In_bs<-list()
      for (i in 1:nl_bs){
        In_bs[[i]]<-vector()
      }
      for (i in 1:nl_bs){
        for (j in 1:nrow(af_bs[[i]])){
          In_bs[[i]][j]<-sum(exp1_bs[[i]][j]+exp2_bs[[i]][j])
        }
      }
      In_sum_bs<-vector()
      for (i in 1:nl_bs){
        In_sum_bs[i]<-sum(In_bs[[i]])
      }
      results_out_bs<-matrix(rep(NA,(nl_bs*(max(nal_bs)+1))),
                             nrow=nl_bs,ncol=(max(nal_bs)+1))
      for (i in 1:nl_bs){
        results_out_bs[i,1:nal_bs[i]]<- In_bs[[i]]
      }
      results_out_bs[,(max(nal_bs)+1)]<-In_sum_bs
      
      bs_sum[w,]<-results_out_bs[,(max(nal_bs)+1)]
    }
    in_bs_out<-matrix(rep(0,(nl_bs*3)),ncol=3)
    colnames(in_bs_out)<-c("In","Lower_95CI","Upper_95CI")
    rownames(in_bs_out)<-ln_bs
    
    for(i in 1:nl){
      in_bs_out[i,1]<-round(mean(bs_sum[,i], na.rm = TRUE),4)
      in_bs_out[i,2] <- round(quantile(bs_sum[,i], 
                                       probs = 0.025, na.rm = TRUE), 4)
      in_bs_out[i,3] <- round(quantile(bs_sum[,i], 
                                       probs = 0.975, na.rm = TRUE), 4)
    }
  }
  colnames(results_out)<-c(paste("Allele.",1:max(nal),sep=""),"Sum")
  rownames(results_out)<-ln
  results_out[is.na(results_out)]<-""
  if(BS==T && pw_only==F){
    list(In_alleles=results_out,
         in_bs_out=in_bs_out)
  }else if(BS==F){
    list(In_alleles=results_out)
  }else if(BS==T && pw_only==T){
    return(in_bs_out)
  }
}
################################################################################
# in.bs end                                                                    #
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
# readGenepop.user, a usable function for basic population parameters          #
################################################################################
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
  npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                  which(data1[,1]=="pop")))
  pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
              which(data1[,1]=="pop"),(nrow(data1)+1))
  pop_sizes<-vector()
  for(i in 1:npops){
    pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
  }
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
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
    gp=x$gp
    bootstrap=x$bootstrap
    # define file reader
    ###########################################################################
    # Master file reader
    ###########################################################################
    fileReader <- function(infile){
      if(typeof(infile)=="list"){
        return(infile) 
      } else if (typeof(infile)=="character"){
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if(ext == "arp"){
          arp2gen(infile)
          cat("Arlequin file converted to genepop format! \n")
          infile <- paste(flForm[1], ".gen", sep = "")
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        # find number of columns
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if(popLoc[1] == 3){
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:(length(dat)-3)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x){
          x <- unlist(strsplit(x, split = "\\s+"))
          if(is.element("", x)){
            x <- x[- (which(x == ""))]
          }
          if(is.element(",", x)){
            x <- x[- (which(x ==","))]
          }
          if(length(x) != 1 && length(x) != no_col){
            x <- paste(x, collapse = "")
          }
          if(length(x) < no_col){
            tabs <- paste(rep(NA, (no_col - length(x))), sep = "\t", 
                          collapse = "\t")
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
    pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
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
        kks2<-colSums(x$indtyp[[gdData[i]]]*
                        (kkptild-rep(kkpbar,
                                     each = x$npops))^2)/((x$npops-1)*kknbar)
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
difPlot <- function (x, outfile= NULL, interactive = FALSE) {
  x=x
  on=outfile
  inta<-interactive
  #output from divPart
  #require(plotrix)
  
  if(is.null(on) == TRUE && inta == TRUE){
    of = paste(getwd(),"/", sep = "")
  } else {
    suppressWarnings(dir.create(paste(getwd(), "/", 
                                      on, "-[diveRsity]", "/", sep="")))
    of=paste(getwd(),"/",on,"-[diveRsity]","/",sep="")
  }
  
  if(!exists("inta",-1)){
    inta<-FALSE
  }
  if(inta == TRUE) {
    sp.header<-list()
    colleer<-list()
    colleer<<-colorRampPalette(c("blue","white"))
    require(sendplot)
    direct<-of
    pwc<-combn(ncol(x[[3]][[1]]),2)
    pwNames<-paste(colnames(x[[3]][[1]])[pwc[1,]],
                   colnames(x[[3]][[1]])[pwc[2,]],
                   sep=' vs ')
    
    gst_lab <- as.vector(x[[3]][[4]])
    gst_lab <- na.omit(gst_lab)
    collab111<-list()
    #
    if(length(x[[3]]) > 6){
      fst_lab <- as.vector(x[[3]][[8]])
      fst_lab<-na.omit(fst_lab)
    }
    #
    gpst_lab <- as.vector(x[[3]][[5]])
    gpst_lab<-na.omit(gpst_lab)
    #
    Dest_lab <- as.vector(x[[3]][[6]])
    Dest_lab<-na.omit(Dest_lab)
    #
    
    fl_ext<-c(".tif","Dot.png","Dot.tif")
    if (length(x[[3]]) > 6){
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
    
    plot.call <- "image(1:nrow(abx[[3]][[4]]),1:nrow(abx[[3]][[4]]),
    abx[[3]][[4]],ylab='',xlab='',main='Pairwise Gst',xaxt='n',yaxt='n',
    col = colleer(50),las=1,cex.main=3)"
    ##
    plot.extras <- "color.legend(nrow(abx[[3]][[4]])/5,
    nrow(abx[[3]][[4]])/3,
    nrow(abx[[3]][[4]])/4,
    nrow(abx[[3]][[4]])/1.2,
    collab111,
    rect.col=colleer(50),
    gradient='y',
    cex=3)"
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
    if(length(x[[3]]) > 6){
      collab111 <<- c(round(min(fst_lab),3),
                      round(mean(fst_lab),3),
                      round(max(fst_lab),3))
      plot.call <- "image(1:nrow(abx[[3]][[8]]),1:nrow(abx[[3]][[8]]),
      abx[[3]][[8]],ylab = '',xlab = '',xaxt = 'n',yaxt = 'n',
      main = 'Pairwise Fst',col = colleer(50),las = 1,cex.main = 3)"
      ##
      plot.extras <- "color.legend(nrow(abx[[3]][[8]])/5,
      nrow(abx[[3]][[8]])/3,
      nrow(abx[[3]][[8]])/4,
      nrow(abx[[3]][[8]])/1.2,
      collab111,
      rect.col=colleer(50),
      gradient='y',
      cex=3)"
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
    plot.call <- "image(1:nrow(abx[[3]][[5]]),1:nrow(abx[[3]][[5]]),
    abx[[3]][[5]],ylab='',xlab='',xaxt='n',yaxt='n',
    main='Pairwise Gst (Hedrick)',col = colleer(50),las=1,cex.main=3)"
    
    plot.extras <- "color.legend(nrow(abx[[3]][[5]])/5,nrow(abx[[3]][[5]])/3,
    nrow(abx[[3]][[5]])/4,nrow(abx[[3]][[5]])/1.2,collab111,
    rect.col=colleer(50),gradient='y',cex=3)"
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
    plot.call <- "image(1:nrow(abx[[3]][[6]]),1:nrow(abx[[3]][[6]]),
    abx[[3]][[6]],ylab='',xlab='',xaxt='n',yaxt='n',main='Pairwise D (Jost)',
    col = colleer(50),las=1,cex.main=3)"
    plot.extras <- "color.legend(nrow(abx[[3]][[6]])/5,nrow(abx[[3]][[6]])/3,
    nrow(abx[[3]][[6]])/4,nrow(abx[[3]][[6]])/1.2,collab111,
    rect.col=colleer(50),gradient='y',cex=3)"
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
    if(length(x[[3]]) > 6){
      par(mfrow=c(2,2))
    } else {
      par(mfrow=c(3,1))
    }
    colleer<-colorRampPalette(c("blue","white"))
    cols<-colleer(50)
    #Gst
    image(1:nrow(x[[3]][[4]]),
          1:nrow(x[[3]][[4]]),
          x[[3]][[4]],
          ylab="population",
          xlab="population",
          main="Pairwise Gst",
          col=cols,
          las=1)
    gst<-as.vector(x[[3]][[4]])
    gst<-as.vector(na.omit(gst))
    collab111<-c(round(min(gst),3),
                 round(mean(gst),3),
                 round(max(gst),3))
    
    color.legend(nrow(x[[3]][[4]])/5,
                 nrow(x[[3]][[4]])/3,
                 nrow(x[[3]][[4]])/4,
                 nrow(x[[3]][[4]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
    if(length(x[[3]]) > 6){
      #Fst
      image(1:nrow(x[[3]][[8]]),
            1:nrow(x[[3]][[8]]),
            x[[3]][[8]],
            ylab="population",
            xlab="population",
            main="Pairwise Theta",
            col = cols,
            las=1)
      fst<-as.vector(x[[3]][[8]])
      fst<-as.vector(na.omit(fst))
      collab111<-c(round(min(fst),3),round(mean(fst),3),round(max(fst),3))
      
      color.legend(nrow(x[[3]][[8]])/5,
                   nrow(x[[3]][[8]])/3,
                   nrow(x[[3]][[8]])/4,
                   nrow(x[[3]][[8]])/1.2,
                   collab111,
                   cols,
                   gradient="y")
    }
    #Hedrick's Gst
    image(1:nrow(x[[3]][[5]]),
          1:nrow(x[[3]][[5]]),
          x[[3]][[5]],
          ylab="population",
          xlab="population",
          main="Pairwise G'st",
          col = cols)
    gprimest<-as.vector(x[[3]][[5]])
    gprimest<-as.vector(na.omit(gprimest))
    collab111<-c(round(min(gprimest),3),
                 round(mean(gprimest),3),
                 round(max(gprimest),3))
    
    color.legend(nrow(x[[3]][[5]])/5,
                 nrow(x[[3]][[5]])/3,
                 nrow(x[[3]][[5]])/4,
                 nrow(x[[3]][[5]])/1.2,
                 collab111,
                 cols,
                 gradient="y")
    #Jost's D
    image(1:nrow(x[[3]][[6]]),
          1:nrow(x[[3]][[6]]),
          x[[3]][[6]],
          ylab="population",
          xlab="population",
          main="Pairwise Jost's D",
          col = cols,
          las=1)
    D<-as.vector(x[[3]][[6]])
    D<-as.vector(na.omit(D))
    collab111<-c(round(min(D),3),
                 round(mean(D),3),
                 round(max(D),3))
    
    color.legend(nrow(x[[3]][[6]])/5,
                 nrow(x[[3]][[6]])/3,
                 nrow(x[[3]][[6]])/4,
                 nrow(x[[3]][[6]])/1.2,
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

chiCalc <- function(infile = NULL, outfile = NULL, gp = 3, minFreq = NULL){
  inputs <- list(infile = infile, gp = gp, bootstrap = FALSE)
  minFreq <- minFreq
  
  # define read genepop function
  #############################################################################
  # readGenepopX, a function for the generation of basic population parameters #
  #############################################################################
  readGenepopX <- function (x) {
    gp=x$gp
    infile=x$infile
    bootstrap=x$bootstrap
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
    npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                    which(data1[,1]=="pop")))
    pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
                which(data1[,1]=="pop"),(nrow(data1)+1))
    pop_sizes<-vector()
    for(i in 1:npops){
      pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
    }
    pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
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
divOnline <- function(){
    shiny::runApp(system.file('diveRsity-online', package = 'diveRsity'))
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
microPlexer <- function(){
  shiny::runApp(system.file('microPlexer', package = 'diveRsity'))
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
divBasic <- function (infile = NULL, outfile = NULL, gp = 3) {
  infile =  infile
  gp = gp
  on = outfile
  
  # create a results dir
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[diveRsity]","/",sep="")))
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
  }
  
  data1 <- fileReader(infile)
  data1[data1==0]<-NA;data1[data1=="999999"]<-NA;data1[data1=="000000"]<-NA
  #raw_data<-data1
  npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                  which(data1[,1]=="pop")))
  pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
              which(data1[,1]=="pop"),(nrow(data1)+1))
  pop_sizes <- sapply(1:npops, function(i){
    pop_pos[(i+1)] - pop_pos[i]-1
  })
  minSize <- min(pop_sizes) 
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,10)
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
  ###vectorize loci_pop_sizes#################################################
  lps<-function(x){#
    lsp_count<-as.vector(colSums(!is.na(x)))#
    return(lsp_count)#
  }#
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
      actab[[x]][[y]]/(indtyppop[[x]][y]*2)
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
  posGeno <- apply(allele_names, 2, function(x){
    sapply(x, function(y){
      if(length(y) == 0){
        return(NA)
      } else {
        genos <- expand.grid(y, y)
        genos.sort <- t(apply(genos, 1, sort))
        genos <- unique(genos.sort)
        geno <- paste(genos[,1], genos[,2], sep = "") 
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
  list(locus_pop_size = locPopSize,
       Allele_number = obsAlls,
       proportion_Alleles = propAlls,
       Allelic_richness = AR,
       Ho = hetObs,
       He = hetExp,
       HWE = HWE,
       mainTab = writeOut)
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
  if(typeof(infile)=="list"){
    return(infile) 
  } else if (typeof(infile)=="character"){
    flForm <- strsplit(infile, split = "\\.")[[1]]
    ext <- flForm[[length(flForm)]]
    if(ext == "arp"){
      arp2gen(infile)
      cat("Arlequin file converted to genepop format! \n")
      infile <- paste(flForm[1], ".gen", sep = "")
    }
    dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
    # find number of columns
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    no_col <- popLoc[1] - 1
    if(popLoc[1] == 3){
      locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
      dat <- c(dat[1], locs, dat[3:(length(dat)-3)])
    }
    popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
    no_col <- popLoc[1] - 1
    dat1 <- sapply(dat, function(x){
      x <- unlist(strsplit(x, split = "\\s+"))
      if(is.element("", x)){
        x <- x[- (which(x == ""))]
      }
      if(is.element(",", x)){
        x <- x[- (which(x ==","))]
      }
      if(length(x) != 1 && length(x) != no_col){
        x <- paste(x, collapse = "")
      }
      if(length(x) < no_col){
        tabs <- paste(rep(NA, (no_col - length(x))), sep = "\t", 
                      collapse = "\t")
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
    gp=x$gp
    bootstrap=x$bootstrap
    # define file reader
    ###########################################################################
    # Master file reader
    ###########################################################################
    fileReader <- function(infile){
      if(typeof(infile)=="list"){
        return(infile) 
      } else if (typeof(infile)=="character"){
        flForm <- strsplit(infile, split = "\\.")[[1]]
        ext <- flForm[[length(flForm)]]
        if(ext == "arp"){
          arp2gen(infile)
          cat("Arlequin file converted to genepop format! \n")
          infile <- paste(flForm[1], ".gen", sep = "")
        }
        dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
        # find number of columns
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        if(popLoc[1] == 3){
          locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
          dat <- c(dat[1], locs, dat[3:(length(dat)-3)])
        }
        popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
        no_col <- popLoc[1] - 1
        dat1 <- sapply(dat, function(x){
          x <- unlist(strsplit(x, split = "\\s+"))
          if(is.element("", x)){
            x <- x[- (which(x == ""))]
          }
          if(is.element(",", x)){
            x <- x[- (which(x ==","))]
          }
          if(length(x) != 1 && length(x) != no_col){
            x <- paste(x, collapse = "")
          }
          if(length(x) < no_col){
            tabs <- paste(rep(NA, (no_col - length(x))), sep = "\t", 
                          collapse = "\t")
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
    pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,6)
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
  npops<-length(c(which(data1[,1]=="Pop"),which(data1[,1]=="POP"),
                  which(data1[,1]=="pop")))
  pop_pos<- c(which(data1[,1]=="POP"),which(data1[,1]=="Pop"),
              which(data1[,1]=="pop"),(nrow(data1)+1))
  loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
  # Calculate the minimum sample size
  pop_sizes <- sapply(1:npops, function(i){
    pop_pos[(i+1)] - pop_pos[i]-1
  })
  #minSize <- min(pop_sizes) 
  pop_names<-substr(data1[(pop_pos[1:npops]+1),1],1,10)
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
bigDivPart <- function(infile = NULL, outfile = NULL, WC_Fst = FALSE,
                       format = NULL){
  
  
  fstat = WC_Fst
  on = outfile
  if (!is.null(on) && format != "txt" && format != "xlsx"){
    stop("Please provide a valid output file format")
  }
  ############################
  
  # use file reader modified from:
  # mlt-thinks.blogspot
  
  fastScan <- function(fname){
    s <- file.info(fname)$size 
    buf <-  readChar(fname, s, useBytes = TRUE)
    return(strsplit(buf, "\n", fixed = TRUE, 
                    useBytes = TRUE)[[1]])
  }
  
  # read data
  
  dat <- fastScan(fname = infile)
  
  # remove the last line if it is blank
  
  if(length(strsplit(dat[length(dat)], split = "\\s+")[[1]]) == 1){
    dat <- dat[-(length(dat))]
  }
  
  # remove the fastScan function
  rm(fastScan)
  z <- gc()
  rm(z)
  
  ############################
  # set up parallel env
  #library("doParallel")
  #cores <- detectCores()
  
  
  # identify population locations
  popLocation <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(dat))
  # extract loci names
  
  pop_pos <- c(popLocation, (length(dat)+1))
  
  loci_names <- as.vector(sapply(dat[2:(pop_pos[1] - 1)], function(x){
    gsub(pattern = "\\s+", replacement = "", x)
  }))
  
  ########################
  # seperate populations #
  ########################
  # get population sizes
  
  popSizes <- NULL
  for(i in 1:(length(pop_pos) - 1)){
    popSizes[i] <- length((pop_pos[i]+1):(pop_pos[(i+1)] - 1))
  }
  
  # create population subsets
  
  # pop data only
  
  pops <- dat[-(c(1:(popLocation[1]-1), popLocation))]
  
  # calculate the row indexes for each population
  
  popList <- lapply(seq_along(popSizes), function(i){
    if(i == 1){
      indx <- 1:popSizes[i]
    } else {
      indx <- (sum(popSizes[1:(i-1)])+1):((sum(popSizes[1:(i-1)])) + 
                                            popSizes[i])
    }
    return(pops[indx])
  })
  
  npops <- length(popList)
  nloci <- length(loci_names)
  pop_sizes <- popSizes
  
  ####################################
  # remove dat
  rm(dat, pops)
  z <- gc(reset = TRUE)
  rm(z)
  
  
  # determine the genepop format of data
  
  testStr <- strsplit(popList[[1]][1], split = "\\s+")[[1]]
  
  gpEst <- sapply(testStr, function(x){
    if(is.character(x)){
      nchar(x)/2
    } else {
      NA
    }
  })
  
  rm(testStr)
  
  
  # take the mode of testStr 
  
  gp <- as.numeric(names(sort(-table(gpEst)))[1])
  
  ############################
  # organise data into a list of arrays
  prePopList <- lapply(popList, function(x){
    y <- array(data = NA, dim = c(length(x), (nloci+1), 2))
    colnames(y) <- c("ind", loci_names)
    for(j in 1:length(x)){
      data <- strsplit(x[j], split = "\\s+")[[1]]
      if(data[2] == ","){
        data <- data[-2]
      }
      data[data == "NANA"] <- NA
      data[data == "0"] <- NA
      data[data == "000000"] <- NA
      data[data == "999999"] <- NA
      data[data == "-9-9"] <- NA
      data[data == "0000"] <- NA
      y[j, 2:(nloci+1), 1] <- substr(data[2:(nloci+1)], 1, gp)
      y[j, 2:(nloci+1), 2] <- substr(data[2:(nloci+1)], gp + 1, gp * 2)
      y[j, 1, 1] <- data[1]
      y[j, 1, 2] <- data[1]
    }
    return(y)
  })
  rm(popList)
  
  # get individuals names
  ind_names <- lapply(prePopList, function(x){
    return(x[ , 1, 1])
  })
  
  
  # pop names
  
  pop_names <- sapply(ind_names, function(x){
    return(x[1])
  })
  
  
  # non bootstrapped stats
  nb <- bigPreDiv(prePopList, FALSE, nloci, npops, popSizes, fstat)
  
  # generate output structures
  
  #standard stats
  
  stdOut <- data.frame(loci = c(loci_names, "Global"),
                       H_st = c(nb$hst, NA),
                       D_st = c(nb$dst, NA),
                       G_st = c(nb$gst, nb$gst_all),
                       G_hed_st = c(nb$gst_hedrick, 
                                    nb$gst_all_hedrick),
                       D_Jost = c(nb$djost, nb$djost_all))
  
  # estimated stats
  
  if(fstat){
    estOut <- data.frame(loci = c(loci_names, "Global"),
                         Harmonic_N = c(nb$locus_harmonic_N, NA),
                         H_st_est = c(nb$hst_est, NA),
                         D_st_est = c(nb$dst_est, NA),
                         G_st_est = c(nb$gst_est, nb$gst_est_all),
                         G_hed_st = c(nb$gst_est_hedrick, 
                                      nb$gst_est_all_hedrick),
                         D_Jost = c(nb$djost_est, nb$djost_est_all),
                         Fst_WC = nb$fstats[,1], Fit_WC = nb$fstats[,2])
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
  
  # write the file to either excel or text
  
  # create results directory
  if(!is.null(on)){
    suppressWarnings(dir.create(path=paste(getwd(),"/",on,
                                           "-[diveRsity]","/",sep="")))
    of = paste(getwd(), "/", on, "-[diveRsity]", "/", sep = "")
  }
  
  # check if xlsx is installed
  write_res <- is.element("xlsx", installed.packages()[, 1])
  
  if(!is.null(on)){
    if(write_res && format == "xlsx"){
      
      # load xlsx package
      require("xlsx")
      
      # standard stats
      write.xlsx(stdOut, file = paste(of, "[bigDivPart].xlsx", sep = ""),
                 sheetName = "Standard_stats", col.names = TRUE,
                 row.names = FALSE, append = FALSE)
      
      # Estimated stats
      write.xlsx(estOut, file = paste(of,"[bigDivPart].xlsx", sep = ""),
                 sheetName = "Estimated_stats", col.names = TRUE,
                 row.names = FALSE, append = TRUE)
    } else {
      # text file alternatives
      
      # standard
      std <- file(paste(of, "Standard-stats[bigDivPart].txt", sep = ""), "w")
      cat(paste(colnames(stdOut), sep = ""), "\n", sep = "\t", file = std)
      stdOut <- as.matrix(stdOut)
      for(i in 1:nrow(stdOut)){
        cat(stdOut[i, ], "\n", file = std, sep = "\t")
      }
      close(std)
      
      # estimated
      est <- file(paste(of, "Estimated-stats[bigDivPart].txt", sep = ""), "w")
      cat(paste(colnames(estOut), sep = ""), "\n", sep = "\t", file = est)
      estOut <- as.matrix(estOut)
      for(i in 1:nrow(estOut)){
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
  
  # popList
  
  if(bs){
    popList <- lapply(prePopList, function(x){
      boot <- sample(1:length(x[,1,1]), replace = TRUE)
      return(x[boot,(2:(nloci+1)),])
    })
  } else {
    popList <- lapply(prePopList, function(x){
      return(x[,(2:(nloci+1)),])
    })
  }
  
  # count the numbers of individuals typed per population
  
  indtyp <- lapply(popList, function(x){
    apply(x, 2, function(y){
      length(na.omit(y[,1]))
    })
  })
  
  # get unique alleles per pop
  
  alls <- lapply(seq_along(popList), function(i){
    apply(popList[[i]], 2, function(x){
      return(unique(c(x[,1], x[,2])))
    })
  })
  
  
  # get unique alleles across pops
  
  all_alleles <- lapply(1:nloci, function(i){
    alleles <- lapply(alls, function(x){
      return(x[[i]])
    })
    return(sort(unique(unlist(alleles))))
  })
  
  # count all observed allele numbers per population
  # (parallel is slower)
  
  obsAlls <- lapply(popList, function(x){
    apply(x, 2, function(y){
      als <- unique(c(na.omit(y[,1]), na.omit(y[,2])))
      counts <- sapply(als, function(z){
        res <- length(which(y == z))
        return(res)
      })
    })
  })
  
  
  # calculate allele frequencies
  
  allele_freq <- lapply(1:nloci, function(i){
    loc <- matrix(nrow = length(all_alleles[[i]]),
                  ncol = npops)
    rownames(loc) <- all_alleles[[i]]
    for(j in 1:npops){
      o <- obsAlls[[j]][[i]]
      n <- indtyp[[j]][i]
      loc[names(o), j] <- o/(2*n)
    }
    loc[is.na(loc)] <- 0
    return(loc)
  })
  
  # generate harmonic mean pop sizes per locus
  preLoc <- lapply(indtyp, function(x){
    return(1/x)
  })
  
  loci_harm_N <- sapply(1:nloci, function(i){
    loc <- sapply(1:npops, function(j){
      return(preLoc[[j]][i])
    })
    return(npops/sum(loc))
  })
  
  loci_harm_N <- round(loci_harm_N, 2)
  
  # convert indtyp to per locus format
  indtypLoc <- lapply(1:nloci, function(i){
    res <- sapply(1:npops, function(j){
      return(indtyp[[j]][i])
    })
  })
  rm(indtyp)
  
  ###############################################
  #  Calculate Weir and Cockerham's (1984) Fst  #
  ###############################################
  if(fstat){            
    badData <- sapply(indtypLoc, function(y){
      is.element(0, y)
    })
    if(sum(badData) > 0){
      nl <- nloci - (sum(badData))
    } else{
      nl <- nloci
    }
    gdData<-which(!badData)
    badData<-which(badData)
    
    # create all genot object
    all_genot <- array(data = NA, dim = c(sum(ps), length(gdData), 1))
    for(i in 1:npops){
      if(i == 1){
        res <- apply(popList[[i]], 2, function(y){
          return(paste0(y[,1], y[,2]))
        })
        all_genot[1:ps[i], ,] <- res[,gdData]
        rm(res)
      } else {
        res <- apply(popList[[i]], 2, function(y){
          return(paste0(y[,1], y[,2]))
        })
        all_genot[(sum(ps[1:(i-1)]) + 1): sum(ps[1:i]), , ] <- res[,gdData]
        rm(res)
      }
    }
    
    all_genot[all_genot == "NANA"] <- NA
    
    # count Genotypes
    #cl <- makeCluster(cores)      
    genoCount <- apply(all_genot, 2, table)
    #stopCluster(cl)
    
    # reformat genoCount names
    nameFormat <- function(x){
      nms <- names(x)
      lgth <- nchar(nms[1])
      newNms <- sapply(nms, function(y){
        paste(substr(y, 1, lgth/2), "/", substr(y, (lgth / 2) + 1, lgth), 
              sep = "")
      })
      names(x) <- newNms
      return(x)
    }
    # run
    genoCount <- lapply(genoCount, nameFormat)
    
    
    # calculate mean heterozygosity per locus
    h_sum <- list()
    for(i in 1:length(gdData)){
      h_sum[[i]] <- vector()
      cnSplit <- strsplit(names(genoCount[[i]]), "/")
      for(j in 1:length(all_alleles[[gdData[i]]])){
        het_id1 <- lapply(cnSplit, is.element, 
                          all_alleles[[gdData[i]]][j])
        het_id2 <- lapply(het_id1, sum)
        het_id1 <- which(het_id2 == 1)
        h_sum[[i]][j] <- sum(genoCount[[i]][het_id1])
      }
    }
    indtyp_tot <- lapply(indtypLoc, sum)
    
    kk_hsum <- lapply(1:ncol(all_genot), function(i){
      list(h_sum[[i]], indtyp_tot[[gdData[i]]])
    })
    
    kk_hbar<-lapply(kk_hsum, function(x){
      return(x[[1]]/x[[2]])
    })
    
    pdat <- lapply(1:length(all_genot[1,,1]), function(i){
      list(allele_freq[[gdData[i]]], indtypLoc[[gdData[i]]])
    })
    
    kk_p <- lapply(pdat, function(x){
      if(is.null(x[[1]]) == FALSE){
        apply(x[[1]], 1, function(y){
          y*(2*x[[2]])
        })
      }
    })
    
    res <- matrix(0, (nloci+1), 2)
    colnames(res) <- c("Fst_WC","Fit_WC")
    #rownames(res) <- c(loci_names, "All")
    A <- vector()
    a <- vector()
    b <- vector()
    c <- vector()
    for(i in 1:length(gdData)){
      kknbar <- indtyp_tot[[gdData[i]]]/npops
      kknC <- (indtyp_tot[[gdData[i]]] - sum(indtypLoc[[gdData[i]]] ^ 2) / 
                 indtyp_tot[[gdData[i]]]) / (npops - 1)
      kkptild <- kk_p[[i]]/(2*indtypLoc[[gdData[i]]])
      kkptild[kkptild == "NaN"] <- NA
      kkpbar <- colSums(kk_p[[i]])/(2 * indtyp_tot[[gdData[i]]])
      kks2 <- colSums(indtypLoc[[gdData[i]]] * 
                        (kkptild - rep(kkpbar, each = npops)) ^ 2) / 
        ((npops - 1) * kknbar)
      kkA <- kkpbar * (1 - kkpbar) - (npops - 1) * kks2 / npops
      kka <- kknbar * (kks2 - (kkA - (kk_hbar[[i]] / 4)) / 
                         (kknbar - 1)) / kknC
      kkb <- kknbar * (kkA - (2 * (kknbar - 1)) * kk_hbar[[i]] / 
                         (4 * kknbar)) / (kknbar - 1)
      kkc <- kk_hbar[[i]] / 2
      A[i] <- sum(kkA, na.rm = TRUE)
      a[i] <- sum(kka, na.rm = TRUE)
      b[i] <- sum(kkb, na.rm = TRUE)
      c[i] <- sum(kkc, na.rm = TRUE)
      res[gdData[i], "Fst_WC"] <- round(sum(kka) / 
                                          sum(kka + kkb + kkc), 4)
      res[gdData[i], "Fit_WC"] <- round(1 - sum(kkc) / 
                                          sum(kka + kkb + kkc),4)
    }
    
    res[res == "NaN"] <- NA
    res[res == 0.000] <- NA
    sumA <- sum(A, na.rm = TRUE)
    suma <- sum(a, na.rm = TRUE)
    sumb <- sum(b, na.rm = TRUE)
    sumc <- sum(c, na.rm = TRUE)
    res[(nloci+1), "Fst_WC"] <- round(suma / (suma +sumb + sumc), 4)
    res[(nloci+1), "Fit_WC"] <- round(1 - sumc / 
                                        (suma + sumb + sumc), 4)
    z <- gc(reset = TRUE)
    rm(z)
    fst <- res
    rm(res)
  }
  #############
  #  end fst  #
  #############    
  
  
  # calculate observed heterozygosity
  
  ho <- lapply(popList, function(x){
    apply(x, 2, function(y){
      1 - (sum(na.omit(y[ ,1] == y[ ,2])) / length(na.omit(y[,1])))
    })
  })
  
  # calculate expected heterozygosity
  
  he <- t(sapply(allele_freq, function(x){
    apply(x, 2, function(y){
      return(1 - sum(y^2))
    })
  }))
  
  # mean frequency
  mf <- lapply(allele_freq, function(x){
    rowSums(x)/ncol(x)
  })
  
  # mean expected heterozygosity
  ht <- sapply(mf, function(x){
    1 - sum(x^2)
  })
  
  ############################
  #   calculate locus stats  #
  ############################
  
  hs <- rowSums(he) / npops
  hs_est <- hs * ((2 * loci_harm_N) / ((2 * loci_harm_N) - 1))
  ht_est <- ht + (hs_est / (2 * loci_harm_N * npops))
  # replace missing data
  ht_est[is.nan(ht_est)] <- NA
  hst <- round((ht - hs) / (1 - hs), 4)
  dst <- round(ht - hs, 4)
  gst <- round(dst / ht, 4)
  # replace missing data
  gst[is.nan(gst)] <- NA
  djost <- round((dst / (1 - hs)) * (npops / (npops - 1)), 4)
  # replace missing data
  djost[djost == 0] <- NA
  hst_est <- round((ht_est - hs_est) / (1 - hs_est), 4)
  dst_est <- round(ht_est - hs_est, 4)
  gst_est <- round(dst_est / ht_est, 4)
  # replace missing data
  gst_est[is.nan(gst_est)] <- NA
  gst_max <- ((npops - 1) * (1 - hs)) / (npops - 1 + hs)
  gst_est_max <- (((npops - 1) * (1 - hs_est)) / (npops - 1 + hs_est))
  gst_hedrick <- round(gst / gst_max, 4)
  gst_est_hedrick <- round(gst_est / gst_est_max, 4)
  gst_est_hedrick[gst_est_hedrick > 1] <- 1
  djost_est <- round((npops / (npops - 1)) * ((ht_est - hs_est) / 
                                                (1 - hs_est)), 4)
  # replace missing data
  djost_est[djost_est == 0] <- NA
  
  ############################
  #   calculate across loci  #
  ############################
  # standard
  ht_mean <- round(mean(ht, na.rm = TRUE), 4)
  hs_mean <- round(mean(hs), 4)
  gst_all <- round((ht_mean - hs_mean) / ht_mean, 4)
  gst_all_max <- round(((npops - 1) * (1 - hs_mean)) / 
                         (npops - 1 + hs_mean), 4)
  gst_all_hedrick <- round(gst_all / gst_all_max, 4)
  djost_all <- round(((ht_mean - hs_mean) / (1 - hs_mean)) * 
                       (npops / (npops - 1)), 4)
  # estimated
  hs_est_mean <- mean(hs_est, na.rm = TRUE)
  ht_est_mean <- mean(ht_est, na.rm = TRUE)
  gst_est_all <- round((ht_est_mean - hs_est_mean) / ht_est_mean, 4)
  gst_est_all_max <- round((((npops - 1) * (1 - hs_est_mean)) / 
                              (npops - 1 + hs_est_mean)), 4)
  gst_est_all_hedrick <- round(gst_est_all / gst_est_all_max, 4)
  gst_est_all_hedrick[gst_est_all_hedrick > 1] <- 1
  if (nloci == 1){
    djost_est_all <- round(djost_est, 4)
  } else {
    djost_est_all <- round(1 / (1 / mean(djost_est, na.rm = TRUE) + 
                                  (var(djost_est, na.rm = TRUE) * 
                                     (1/mean(djost_est, na.rm = TRUE)) ^ 3)), 4)
  }    
  djost_est[djost_est==0]<-NaN
  djost[djost==0]<-NaN
  # END    
  
  # return results
  if(fstat){
    list(hst = hst, dst = dst, gst = gst, gst_hedrick = gst_hedrick,
         djost = djost, locus_harmonic_N = loci_harm_N, 
         hst_est = hst_est, dst_est = dst_est, gst_est = gst_est,
         gst_est_hedrick = gst_est_hedrick, djost_est = djost_est,
         gst_all = gst_all, gst_all_hedrick = gst_all_hedrick,
         djost_all = djost_all, gst_est_all = gst_est_all,
         gst_est_all_hedrick = gst_est_all_hedrick,
         djost_est_all = djost_est_all, fstats = fst)
    
  } else{
    list(hst = hst, dst = dst, gst = gst, gst_hedrick = gst_hedrick, 
         djost = djost, locus_harmonic_N = loci_harm_N, 
         hst_est = hst_est, dst_est = dst_est, gst_est = gst_est,
         gst_est_hedrick = gst_est_hedrick, djost_est = djost_est,
         gst_all = gst_all, gst_all_hedrick = gst_all_hedrick,
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
arp2gen <- function(infile){
  # define a fastscan function
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
  if(length(outfile) > 2){
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

divMigrate <- function(infile = NULL, stat = "d_jost"){
  # check file format
  cat("Caution! The method used in this function is experimental. \n")
#   flForm <- strsplit(infile, split = "\\.")[[1]]
#   ext <- flForm[[length(flForm)]]
#   if(ext == "arp"){
#     arp2gen(infile)
#     cat("Arlequin file converted to genepop format!")
#     infile <- paste(flForm[1], ".gen", sep = "")
#   }
  dat <- fileReader(infile)
  rownames(dat) <- NULL
  dat <- as.matrix(dat)
  # determine genepop format
  p1 <- which(toupper(dat[,1]) == "POP")[1] + 1
  gp <- as.numeric(names(sort(-table(sapply(dat[p1, - 1], nchar)/2)))[1])
  dat <- as.data.frame(dat)
  rawData <- readGenepop(dat, gp = gp)
  npops <- rawData$npops
  nloci <- rawData$nloci
  # generate pairwise hypothetical matrices (use allele_freq)
  pw <- combn(npops, 2)
  # calculate ht and hs
  hths <- lapply(rawData$allele_freq, pwDivCalc, pw = pw, npops = npops)
  # seperate ht and hs matrices
  ht <- lapply(hths, "[[", 1)
  hs <- lapply(hths, "[[", 2)
  # tidy
  rm(hths)
  z <- gc()
  rm(z)
  # find the mean (use the Reduce function)
  ht_mean <- Reduce(`+`, ht)/nloci
  hs_mean <- Reduce(`+`, hs)/nloci
  # calculate Dst
  dst <- ht_mean - hs_mean
  # calculate gst
  gst <- dst/ht_mean
  gst[gst < 0.0 | is.na(gst)] <- 0
  # calculate D(Jost)
  d_jost <- ((dst)/(1-hs_mean)) * 2
  # calculate relative migration from d_jost
  d_mig <- (1 - d_jost)/d_jost
  # replace missing and negative values with 0
  d_mig[d_mig < 0 | is.na(d_mig)] <- 0
  dimnames(d_mig) <- list(paste("P", 1:npops),
                          paste("P", 1:npops))
  # standardize
  d_mig <- d_mig/max(d_mig, na.rm = TRUE)
  # test gst migration rate
  gst_mig <- 0.5 * ((1/gst) - 1)
  # fix inf
  gst_mig[is.infinite(gst_mig)] <- 0
  # standardise
  gst_mig <- gst_mig/max(gst_mig, na.rm = TRUE)
  # replace missing and negative values with 0
  gst_mig[gst_mig < 0 | is.na(gst_mig)] <- 0
  dimnames(gst_mig) <- list(paste("P", 1:npops),
                            paste("P", 1:npops))
  # test plot
  #library("qgraph")
  if(length(stat) == 2){
    par(mfrow = c(2, 1 ))
    qgraph(gst_mig, posCol = "black")
    title(expression("G"["st"]))
    qgraph(d_mig, posCol = "black")
    title(expression("D"["Jost"]))
    par(mfrow = c(1,1))
  } else if(stat == "gst"){
    qgraph(gst_mig, posCol = "black")
    title(expression("G"["st"]))
  } else if(stat == "d_jost"){
    qgraph(d_mig, posCol = "black")
    title(expression("D"["Jost"]))
  }
  list(D_mig =d_mig,
       Gst_mig = gst_mig)
}
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
    gamma <- sum(sqrt(abs(x[,1] * x[,2])))^-1 
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
#
#
#
#
#
################################################################################
#################################     END ALL       ############################
################################################################################
