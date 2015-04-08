################################################################################
# fstOnly: a memory efficient function to calculate WC Fst and Fit
################################################################################
#' @export
fstOnly <- function(infile = NULL, outfile = NULL, gp = 3, 
                    bs_locus = FALSE, bs_pairwise = FALSE, 
                    bootstraps = 0, parallel = FALSE){
  .Deprecated(new = "diffCalc", msg = "This function is no longer in use. Please use 'diffCalc' instead, \nSee ?diffCalc for usage details.", 
              old = "fstOnly")
}
# fstOnly <- function(infile = NULL, outfile = NULL, gp = 3, 
#                     bs_locus = FALSE, bs_pairwise = FALSE, 
#                     bootstraps = 0, parallel = FALSE){
#   # create a directory
#   if(!is.null(outfile)){
#     suppressWarnings(dir.create(path=paste(getwd(), "/", outfile,
#                                            "-[fstWC]", "/",sep="")))
#   }
#   
#   # define the fstWC function
#   #############################################################################
#   # fstWC: a function co calculate weir and cockerhams fis, fit, and fst
#   #############################################################################
#   fstWC<-function(x){
#     badData <- sapply(x$indtyp, function(y){
#       is.element(0, y)
#     })
#     if(sum(badData) > 0){
#       nl <- x$nloci - (sum(badData))
#     } else{
#       nl <- x$nloci
#     }
#     gdData<-which(!badData)
#     badData<-which(badData)
#     if (nl == 1) {
#       all_genot<-x$pop_list[[1]][,gdData]
#       if(x$npops > 1){
#         for(i in 2:x$npops){
#           all_genot <- c(all_genot, x$pop_list[[i]][,gdData])
#         }
#       }
#       all_genot <- matrix(all_genot, ncol = 1)
#     } else {
#       all_genot<-matrix(x$pop_list[[1]][,gdData], ncol = length(gdData))
#       if(x$npops > 1){
#         for(i in 2:x$npops){
#           all_genot<-rbind(all_genot, x$pop_list[[i]][,gdData])
#         }
#       }
#     }
#     genot<-apply(all_genot,2,unique)
#     genot<-lapply(genot, function(x){
#       if (sum(is.na(x))>0){
#         y<-which(is.na(x)==TRUE)
#         x_new<-x[-y]
#         return(x_new)
#       } else {
#         return(x)
#       }
#     })
#     #count genotypes
#     
#     genoCount<-list()
#     for(i in 1:ncol(all_genot)){
#       genoCount[[i]]<-matrix(0,ncol=length(genot[[i]]))
#       for(j in 1:length(genot[[i]])){
#         genoCount[[i]][,j]<-length(which(all_genot[,i] == genot[[i]][j]))
#       }
#       if (x$gp==3){
#         colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,3),"/",
#                                         substr(genot[[i]],4,6),sep="")
#       } else if (x$gp==2){
#         colnames(genoCount[[i]])<-paste(substr(genot[[i]],1,2),"/",
#                                         substr(genot[[i]],3,4),sep="")
#       }
#     }
#     
#     h_sum<-list()
#     for(i in 1:ncol(all_genot)){
#       h_sum[[i]]<-vector()
#       cnSplit<-strsplit(colnames(genoCount[[i]]),"/")
#       for(j in 1:length(x$all_alleles[[gdData[i]]])){
#         het_id1<-lapply(cnSplit, is.element, x$all_alleles[[gdData[i]]][j])
#         het_id2<-lapply(het_id1, sum)
#         het_id2<-as.vector(het_id2)
#         het_id3<-which(het_id2==1)
#         h_sum[[i]][j]<-sum(genoCount[[i]][1,het_id3])
#       }
#     }
#     indtyp_tot<-lapply(x$indtyp, sum)
#     kk_hsum <- lapply(1:ncol(all_genot), function(i){
#       list(h_sum[[i]], indtyp_tot[[gdData[i]]])
#     })
#     kk_hbar<-lapply(kk_hsum, function(x){
#       return(x[[1]]/x[[2]])
#     })
#     
#     pdat <- lapply(1:ncol(all_genot), function(i){
#       list(x$allele_freq[[gdData[i]]], x$indtyp[[gdData[i]]])
#     })
#     
#     kk_p<-lapply(pdat, function(x){
#       if(is.null(x[[1]])==FALSE){
#         apply(x[[1]], 1, function(y){
#           y*(2*x[[2]])
#         })
#       }
#     })
#     res<-matrix(0,(x$nloci+1),3)
#     colnames(res)<-c("Fis_WC","Fst_WC","Fit_WC")
#     rownames(res)<-c(x$loci_names, "All")
#     A<-vector()
#     a<-vector()
#     b<-vector()
#     c<-vector()
#     for(i in 1:ncol(all_genot)){
#       kknbar<-indtyp_tot[[gdData[i]]]/x$npops
#       kknC<-(indtyp_tot[[gdData[i]]]-sum(x$indtyp[[gdData[i]]]^2)/
#                indtyp_tot[[gdData[i]]])/(x$npops-1)
#       kkptild<-kk_p[[i]]/(2*x$indtyp[[gdData[i]]])
#       kkptild[kkptild=="NaN"]<-NA
#       kkpbar<-colSums(kk_p[[i]])/(2*indtyp_tot[[gdData[i]]])
#       kks2<-colSums(x$indtyp[[gdData[i]]]*
#                       (kkptild-rep(kkpbar,each = x$npops))^2)/((x$npops-1)*
#                                                                  kknbar)
#       kkA<-kkpbar*(1-kkpbar)-(x$npops-1)*kks2/x$npops
#       kka<-kknbar*(kks2-(kkA-(kk_hbar[[i]]/4))/(kknbar-1))/kknC
#       kkb<-kknbar*(kkA-(2*(kknbar-1))*kk_hbar[[i]]/(4*kknbar))/(kknbar-1)
#       kkc<-kk_hbar[[i]]/2
#       A[i]<-sum(kkA)
#       a[i]<-sum(kka)
#       b[i]<-sum(kkb)
#       c[i]<-sum(kkc)
#       res[gdData[i],"Fis_WC"]<- round(1-sum(kkc)/sum(kkb+kkc),4)
#       res[gdData[i],"Fst_WC"]<- round(sum(kka)/sum(kka+kkb+kkc),4)
#       res[gdData[i],"Fit_WC"]<- round(1-sum(kkc)/sum(kka+kkb+kkc),4)
#     }
#     res[res=="NaN"]<-NA
#     res[res==0.000]<-NA
#     sumA<-sum(na.omit(A))
#     suma<-sum(na.omit(a))
#     sumb<-sum(na.omit(b))
#     sumc<-sum(na.omit(c))
#     res[(x$nloci+1),"Fis_WC"]<-round(1-sumc/(sumb+sumc),4)
#     res[(x$nloci+1),"Fst_WC"]<-round(suma/(suma+sumb+sumc),4)
#     res[(x$nloci+1),"Fit_WC"]<-round(1-sumc/(suma+sumb+sumc),4)
#     #res[is.na(res)]<-NaN
#     list(Fstats=res,
#          multiLoc<-res[(x$nloci+1),])
#   }
#   #############################################################################
#   # end fstWC
#   #############################################################################
#   #
#   #
#   #
#   # define the readGenepopX function
#   readGenepopX <- function (x) {
#     infile=x$infile
#     #gp=x$gp
#     bootstrap=x$bootstrap
#     # define file reader
#     ###########################################################################
#     # Master file reader
#     ###########################################################################
#     fileReader <- function(infile){
#       if (typeof(infile) == "list") {
#         return(infile)
#       } else if (typeof(infile) == "character") {
#         flForm <- strsplit(infile, split = "\\.")[[1]]
#         ext <- flForm[[length(flForm)]]
#         if (ext == "arp") {
#           convRes <- arp2gen(infile)
#           if (!is.null(convRes)) {
#             cat("Arlequin file converted to genepop format! \n")
#             infile <- paste(flForm[1], ".gen", sep = "")
#           } else {
#             infile <- paste(flForm[1], ".gen", sep = "")
#           }
#         }
#         dat <- scan(infile, sep = "\n", what = "character", quiet = TRUE)
#         if(length(strsplit(dat[4], split = "\\s+")[[1]][-1]) > 1){
#           locs <- strsplit(dat[2], split = "\\s+")[[1]]
#           if(length(locs != 1)){
#             locs <- strsplit(dat[2], split = ",")[[1]]
#           }
#           locs <- as.character(sapply(locs, function(x){
#             x <- strsplit(x, split = "")[[1]]
#             if(is.element(",", x)){
#               x <- x[-(which(x == ","))]
#             }
#             return(paste(x, collapse = ""))
#           }))
#           dat <- c(dat[1], locs, dat[-(1:2)])
#         }
#         
#         
#         popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
#         no_col <- popLoc[1] - 1
#         if (popLoc[1] == 3) {
#           locs <- unlist(strsplit(dat[2], split = c("\\,", "\\s+")))
#           dat <- c(dat[1], locs, dat[3:length(dat)])
#         }
#         popLoc <- grep("^([[:space:]]*)pop([[:space:]]*)$", tolower(dat))
#         no_col <- popLoc[1] - 1
#         dat1 <- sapply(dat, function(x) {
#           x <- unlist(strsplit(x, split = "\\s+"))
#           if (is.element("", x)) {
#             x <- x[-(which(x == ""))]
#           }
#           if (is.element(",", x)) {
#             x <- x[-(which(x == ","))]
#           }
#           if (length(x) != 1 && length(x) != no_col) {
#             x <- paste(x, collapse = "")
#           }
#           if (length(x) < no_col) {
#             tabs <- paste(rep(NA, (no_col - length(x))), 
#                           sep = "\t", collapse = "\t")
#             line <- paste(x, tabs, sep = "\t")
#             line <- unlist(strsplit(line, split = "\t"))
#             return(line)
#           } else {
#             return(x)
#           }
#         })
#       }
#       out <- as.data.frame(t(dat1))
#       rownames(out) <- NULL
#       return(out)
#     }
#     data1 <- fileReader(infile)
#     if(is.null(x$gp)){
#       rownames(data1) <- NULL
#       data1 <- as.matrix(data1)
#       # determine genepop format
#       p1 <- which(toupper(data1[,1]) == "POP")[1] + 1
#       gp <- as.numeric(names(sort(-table(sapply(data1[p1,(-1)], nchar)/2)))[1])
#       data1 <- as.data.frame(data1)
#     } else {
#       gp <- x$gp
#     }
#     if(gp == 3){
#       data1[data1==0]<-NA
#       data1[data1=="999999"]<-NA
#       data1[data1=="000000"]<-NA
#       data1[data1=="NANA"]<-NA
#     } else if(gp == 2){
#       data1[data1==0]<-NA
#       data1[data1=="9999"]<-NA
#       data1[data1=="0000"]<-NA
#       data1[data1=="NA"]<-NA
#     }    
#     raw_data<-data1
#     npops<-length(which(toupper(data1[,1]) == "POP"))
#     pop_pos<- c(which(toupper(data1[,1]) == "POP"), (nrow(data1)+1))
#     pop_sizes<-vector()
#     for(i in 1:npops){
#       pop_sizes[i]<- pop_pos[(i+1)] - pop_pos[i]-1
#     }
#     pop_names<-as.character(data1[(pop_pos[1:npops]+1),1])
#     pop_weights<- 1/pop_sizes
#     
#     n_harmonic<-npops/sum(pop_weights)
#     
#     N<-pop_sizes
#     
#     nloci<- (pop_pos[1]-2)
#     loci_names<-as.vector(data1[2:(pop_pos[1]-1),1])
#     pop_list<-list()
#     for (i in 1:npops){
#       pop_list[[i]]<-as.matrix(data1[(pop_pos[i]+1):(pop_pos[(i+1)]-1),
#                                      2:(nloci+1)])
#     }
#     # check if all populations have at least some data at loci
#     extCheck <- sapply(1:length(pop_list), function(i){
#       sum(is.na(pop_list[[i]])) == nloci * pop_sizes[i]
#     })
#     if (sum(extCheck) > 0){
#       npops <- npops - sum(extCheck)
#       pop_list <- pop_list[-(which(extCheck == TRUE))]
#       pop_sizes <- pop_sizes[-(which(extCheck == TRUE))]
#       pop_names <- pop_names[-(which(extCheck == TRUE))]
#       pop_weights <- pop_weights[-(which(extCheck == TRUE))]
#       N <- N[-(which(extCheck == TRUE))]
#       #raw_data fix
#       noPop <- which(extCheck == TRUE)
#       indexer <- lapply(noPop, function(i){
#         (pop_pos[i] + 1):(pop_pos[(i+1)])
#       })
#       indexer <- unlist(indexer)
#       raw_data <- raw_data[-(indexer), ]    
#     }  
#     if (gp==3) {
#       plMake<-function(x){
#         out <- matrix(sprintf("%06g",as.numeric(x)),
#                       nrow = nrow(x), ncol = ncol(x))
#         if (Sys.info()["sysname"] == "Darwin"){
#           out[out == "0000NA"] <- "    NA"
#         }
#         return(out)
#       }
#     } else if (gp==2) {
#       plMake<-function(x){
#         out <- matrix(sprintf("%04g",as.numeric(x)),
#                       nrow = nrow(x), ncol = ncol(x))
#         if (Sys.info()["sysname"] == "Darwin"){
#           out[out == "00NA"] <- "  NA"
#         }
#         return(out)
#       }
#     }
#     suppressWarnings(pop_list<-lapply(pop_list, plMake))
#     
#     if (gp == 3){
#       for(i in 1:npops){
#         pop_list[[i]][pop_list[[i]] == "    NA"]<-NA
#       }
#     } else if (gp == 2){
#       for(i in 1:npops){
#         pop_list[[i]][pop_list[[i]] == "  NA"] <-NA
#       }
#     }
#     
#     if(bootstrap == T){
#       bs<-function(x){
#         return(matrix(x[sample(nrow(x),replace=TRUE), ],ncol=ncol(x)))
#       }
#       pop_list<-lapply(pop_list, bs)
#     }  
#     
#     ###vectorize loci_pop_sizes###############################################
#     
#     lps<-function(x){#
#       lsp_count<-as.vector(colSums(!is.na(x)))#
#       return(lsp_count)#
#     }#
#     pre_loci_pop_sizes<-lapply(pop_list,lps)#
#     pls<-matrix(ncol=nloci,nrow=npops)#
#     for(i in 1:length(pre_loci_pop_sizes)){#
#       pls[i,]<-pre_loci_pop_sizes[[i]]#
#     }#
#     #convert pls to loci_pop_sizes format
#     loci_pop_sizes<-split(pls,col(pls))
#     
#     
#     #vectorized loci_pop_weights##############################################
#     
#     pre_loc_weights<- 1/pls
#     loci_pop_weights1<-split(pre_loc_weights,col(pre_loc_weights))
#     loci_harm_N<-npops/colSums(pre_loc_weights)
#     
#     #end vectorized loci_pop_weights##########################################
#     
#     ###vectorize pop_alleles##################################################
#     if (gp==3){
#       pl_ss<-function(x){  # where x is object pop_list
#         pl<-list()
#         pl[[1]]<-matrix(substr(x,1,3),ncol=nloci)
#         pl[[2]]<-matrix(substr(x,4,6),ncol=nloci)
#         return(pl)
#       }
#     } else {
#       pl_ss<-function(x){  # where x is object pop_list
#         pl<-list()
#         pl[[1]]<-matrix(substr(x,1,2),ncol=nloci)
#         pl[[2]]<-matrix(substr(x,3,4),ncol=nloci)
#         return(pl)
#       }
#     }
#     pop_alleles<-lapply(pop_list,pl_ss)
#     #end vectorize pop_alleles################################################
#     
#     #vectorize allele_names###################################################
#     
#     alln<-function(x){ # where x is the object pop_alleles (returned by pl_ss())
#       res<-list()
#       for(i in 1:ncol(x[[1]])){
#         res[i]<-list(sort(unique(c(x[[1]][,i],x[[2]][,i])),decreasing=F))
#       }
#       return(res)
#     }
#     
#     allele_names<-lapply(pop_alleles,alln)
#     
#     
#     loci_combi<-allele_names[[1]]
#     for(j in 1:nloci){
#       for(i in 2:npops){
#         loci_combi[[j]]<-c(loci_combi[[j]],allele_names[[i]][[j]])
#       }
#     }
#     
#     #all_alleles vectorized###################################################
#     
#     aaList<-function(x){
#       return(sort(unique(x,decreasing=FALSE)))
#     }
#     all_alleles<-lapply(loci_combi,aaList)
#     
#     #end all_alleles vectorized###############################################
#     
#     aa<-all_alleles
#     aa<-lapply(aa, FUN=`list`, npops)
#     afMatrix<-function(x){
#       np<-x[[2]]
#       z<-matrix(rep(0,(np*length(x[[1]]))),ncol=np, nrow=length(x[[1]]))
#       rownames(z)<-x[[1]]
#       return(z)
#     }
#     allele_freq<-lapply(aa,afMatrix)
#     
#     
#     #combine pop_alleles
#     parbind<-function(x){
#       rbind(x[[1]],x[[2]])
#     }
#     pa1<-lapply(pop_alleles, parbind)
#     #create a function to tabulate the occurance of each allele
#     afTab<-function(x){
#       lapply(1:ncol(x), function(i){
#         return(table(x[,i]))
#       })
#     }
#     actab<-lapply(pa1, afTab)
#     
#     afs<-function(x){
#       afsint<-function(y){
#         length(na.omit(y))/2
#       }
#       apply(x,2,afsint)
#     }
#     indtyppop<-lapply(pa1,afs)
#     #calculate allele frequencies
#     afCalcpop<-lapply(1:length(actab), function(x){
#       lapply(1:length(actab[[x]]),function(y){
#         actab[[x]][[y]]/(indtyppop[[x]][y]*2)
#       })
#     })
#     #assign allele freqs to frequency matrices
#     obs_count<-allele_freq
#     for(i in 1:npops){
#       for(j in 1:nloci){
#         allele_freq[[j]][names(afCalcpop[[i]][[j]]),i]<-afCalcpop[[i]][[j]]
#         obs_count[[j]][names(actab[[i]][[j]]),i]<-actab[[i]][[j]]
#       }
#     }
#     
#     indtyp<-list()
#     for(i in 1:nloci){
#       indtyp[[i]]<-vector()
#     }
#     for(i in 1:npops){
#       for(j in 1:nloci){
#         indtyp[[j]][i]<-indtyppop[[i]][j]
#       }
#     }
#     
#     if(bootstrap==T){
#       ind_vectors<-list()
#       for(i in 1:npops){
#         ind_vectors[[i]]<-noquote(paste(rep(i,pop_sizes[i]),",",sep=""))
#       }
#       pre_data<-matrix(rep("\t",((nloci+1)*(nloci+1))),
#                        ncol=(nloci+1))
#       pre_data[1,]<-c("Title",rep("\t",nloci))
#       for(i in 2:(nloci+1)){
#         pre_data[i,1]<-loci_names[(i-1)]
#       }
#       pop_data<-list()
#       for(i in 1:npops){
#         pop_data[[i]]<-matrix(rbind(c("POP",as.vector(rep("\t",nloci))),
#                                     cbind(ind_vectors[[i]],pop_list[[i]])),
#                               ncol=(nloci+1))
#       }
#       bs_data_file<-matrix(rbind(pre_data,pop_data[[1]]),ncol=(nloci+1))
#       for(i in 2:npops){
#         bs_data_file<-matrix(rbind(bs_data_file,pop_data[[i]]),ncol=(nloci+1))
#       }
#       bs_data_file<-data.frame(bs_data_file)
#     }
#     nalleles<-vector()
#     for(i in 1:nloci){
#       nalleles[i]<- nrow(allele_freq[[i]])
#     }
#     ##########################################################################
#     list(pop_list = pop_list,
#          npops = npops,
#          nloci = nloci,
#          pop_sizes = pop_sizes,
#          pop_alleles = pop_alleles,
#          all_alleles = all_alleles,
#          allele_freq = allele_freq,
#          loci_harm_N = loci_harm_N,
#          loci_names = loci_names,
#          pop_names = pop_names,
#          indtyp = indtyp,
#          gp = gp)
#   }
#   ############################################################################
#   # readGenepopX end                                                        #
#   ############################################################################
#   
#   #setup a parallel cluster if parallel == TRUE
#   if(parallel){
#     para_pack <- is.element(c("parallel","doParallel","foreach","iterators"),
#                             installed.packages()[,1])
#     if(sum(para_pack) != 4){
#       stop("Please install all required parallel packages")
#     } else {
#       library("doParallel")
#       cores <- detectCores()
#       cl <- makeCluster(cores)
#       registerDoParallel(cl)
#     }
#   }
#   
#   
#   # create the baseline
#   accData <- readGenepopX(list(infile = infile, gp = gp, bootstrap = FALSE))
#   glbFst <- fstWC(accData)
#   if(bs_locus){
#     # calculate locus bootstraps
#     if(bs_locus && parallel){
#       input <- list(infile = infile, gp = gp, bootstrap = TRUE)
#       clusterExport(cl, c("fstWC", "readGenepopX", "input"), 
#                     envir = environment())
#       loc_stats <- parLapply(cl, 1:bootstraps, function(...){
#         gps <- readGenepopX(input)
#         fst <- fstWC(gps)$Fstats
#         return(fst)
#       })
#     } else if(bs_locus && !parallel){
#       input <- list(infile = infile, gp = gp, bootstrap = TRUE)
#       loc_stats <- lapply(1:bootstraps, function(...){
#         gps <- readGenepopX(input)
#         fst <- fstWC(gps)$Fstats
#         return(fst)
#       })
#     }
#     # compile bs_locus results
#     loc_fst <- sapply(loc_stats, function(x){
#       return(x[,"Fst_WC"])
#     })
#     loc_fit <- sapply(loc_stats, function(x){
#       return(x[,"Fit_WC"])
#     })
#     locBS <- list(loc_fst = loc_fst,
#                   loc_fit = loc_fit)
#     locBSci <- lapply(locBS, function(x){
#       apply(x, 1, function(y){
#         return(as.vector(((sd(y)/sqrt(length(y))) * 1.96)))
#       })
#     })
#     locBSout <- lapply(1:2, function(i){
#       if(i == 1){
#         cbind(actual = round(glbFst[[1]][,"Fst_WC"], 4),
#               lower_CI = round(glbFst[[1]][,"Fst_WC"] - locBSci[[i]], 4),
#               upper_CI = round(glbFst[[1]][,"Fst_WC"] + locBSci[[i]], 4))
#         
#       } else if(i == 2){
#         cbind(actual = round(glbFst[[1]][,"Fit_WC"], 4),
#               lower_CI = round(glbFst[[1]][,"Fit_WC"] - locBSci[[i]], 4),
#               upper_CI = round(glbFst[[1]][,"Fit_WC"] + locBSci[[i]], 4))
#       }
#     })
#     # write the locus results
#     locOut <- rbind(c("actual", "lower_CI", "upper_CI"),
#                     locBSout[[1]],
#                     c("--", "--", "--"),
#                     c("","",""),
#                     c("actual", "lower_CI", "upper_CI"),
#                     locBSout[[2]])
#     locNames <- c("Fst", rownames(locBSout[[1]]),
#                   "--", "", "Fit", rownames(locBSout[[1]]))
#     locOut <- cbind(locNames, locOut)
#     dimnames(locOut) <- NULL
#     
#     if (!is.null(outfile)){
#       of = paste(getwd(), "/", outfile, "-[fstWC]", "/", sep = "")
#       if(is.element("xlsx", installed.packages()[, 1])){
#         # write data to excel
#         # Load dependencies
#         library("xlsx")
#         # standard stats
#         write.xlsx(locOut, file = paste(of, "[fstWC].xlsx", sep = ""),
#                    sheetName = "Locus_stats", col.names = FALSE,
#                    row.names = FALSE, append = FALSE)
#       } else {
#         # text file alternatives
#         std <- file(paste(of, "Locus_stats-[fstWC].txt", sep = ""), "w")
#         for(i in 1:nrow(locOut)){
#           cat(locOut[i,], "\n", file = std, sep = "\t")
#         }
#         close(std)
#       }
#     }
#     names(locBSout) <- c("Fst", "Fit")
#   }
#   
#   
#   ##########################################################################
#   ##                              PAIRWISE                                ##
#   ##########################################################################
#   
#   
#   
#   # calculate pairwise bootstraps
#   if(bs_pairwise){
#     pw <- combn(accData$npops,2)
#     pwmat <- pw + 1
#     ind_vectors <- lapply(1:accData$npops, function(x){
#       rep(x, accData$pop_sizes[[x]])}
#     )
#     #      
#     pre_data <- matrix(rep("", ((accData$nloci + 1) * (accData$nloci + 1))),
#                        ncol = (accData$nloci + 1))
#     pre_data[1,] <- rep("", (accData$nloci + 1))
#     #
#     for(i in 2:(accData$nloci + 1)){
#       pre_data[i, 1] <- accData$loci_names[(i-1)]
#     }
#     #
#     pw_data<-list()
#     for (i in 1:ncol(pw)){
#       pw_data[[i]]<-data.frame(rbind(pre_data,
#                                      c("POP",as.vector(rep("",accData$nloci))),
#                                      cbind(ind_vectors[[pw[1,i]]],
#                                            matrix(noquote(accData$pop_list
#                                                           [[pw[1,i]]]),
#                                                   ncol=accData$nloci)),
#                                      c("POP",as.vector(rep("",accData$nloci))),
#                                      cbind(ind_vectors[[pw[2,i]]],
#                                            matrix(noquote(accData$pop_list
#                                                           [[pw[2,i]]]),
#                                                   ncol=accData$nloci))))
#     }
#     # define true stat res obj
#     pw_glb <- matrix(rep(0, (2 * (ncol(pw)))), ncol = 2)
#     
#     # true stat input
#     trueStatIn <- lapply(pw_data, function(x){
#       list(infile = x, gp = gp, bootstrap = FALSE)
#     })
#     
#     # calculate true stats
#     if(parallel){
#       clusterExport(cl, c("readGenepopX", "fstWC"), envir = environment())
#       tparSapply <- function(...) t(parSapply(...))
#       trueStat <- tparSapply(cl, trueStatIn, function(x){
#         input <- readGenepopX(x)
#         return(fstWC(input)[[2]][2:3])
#       })
#     } else {
#       tsapply <- function(...) t(sapply(...))
#       trueStat <- tsapply(cl, trueStatIn, function(x){
#         input <- readGenepopX(x)
#         return(fstWC(input)[[2]][2:3])
#       })
#     }
#     
#     fstBS <- function(x){
#       fstIn <- readGenepopX(x)
#       fstOut <- fstWC(fstIn)
#       return(fstOut[[2]][2:3])
#     }
#     
#     if(parallel){
#       # bootstrap pairwise populations
#       clusterExport(cl, c("pw_data", "gp", "fstWC", "bootstraps", 
#                           "readGenepopX", "fstBS"), envir = environment())
#       pwRES <- parLapply(cl, 1:ncol(pw), function(i){
#         input <- list(infile = pw_data[[i]], gp = gp, bootstrap = TRUE)
#         out <- replicate(bootstraps, fstBS(input))
#         fstCI <- (sd(out[1, ])/sqrt(length(out[1,]))) * 1.96
#         fitCI <- (sd(out[2, ])/sqrt(length(out[2,]))) * 1.96
#         return(c(fst = fstCI, fit = fitCI))
#       })
#       stopCluster(cl)
#     } else {
#       # bootstrap pairwise populations
#       pwRES <- lapply(1:ncol(pw), function(i){
#         input <- list(infile = pw_data[[i]], gp = gp, bootstrap = TRUE)
#         out <- replicate(bootstraps, fstBS(input))
#         fstCI <- (sd(out[1, ])/sqrt(length(out[1,]))) * 1.96
#         fitCI <- (sd(out[2, ])/sqrt(length(out[2,]))) * 1.96
#         return(c(fst = fstCI, fit = fitCI))
#       })
#     }
#     pwOut <- lapply(1:2, function(i){      
#       lapply(1:ncol(pw), function(j){
#         if(i == 1){
#           return(round(c(trueStat[j,"Fst_WC"], 
#                          trueStat[j,"Fst_WC"] - pwRES[[j]]["fst"],
#                          trueStat[j,"Fst_WC"] + pwRES[[j]]["fst"]), 4))
#         } else if(i == 2){
#           return(round(c(trueStat[j,"Fit_WC"], 
#                          trueStat[j,"Fit_WC"] - pwRES[[j]]["fit"],
#                          trueStat[j,"Fit_WC"] + pwRES[[j]]["fit"]),4))
#         }
#       })
#     })
#     pwOut1 <- lapply(pwOut, function(x){
#       out <- as.data.frame(do.call("rbind", x))
#       colnames(out) <- c("actual", "lower_CI", "upper_CI")
#       rownames(out) <- paste(accData$pop_names[pw[1,]], " v ",
#                              accData$pop_names[pw[2,]], sep = "")
#       return(out)
#     })
#     
#     # write pw bootstrap results
#     pwWriteOut <- rbind(c("actual", "lower_CI", "upper_CI"),
#                         pwOut1[[1]],
#                         c("--", "--", "--"),
#                         c("","",""),
#                         c("actual", "lower_CI", "upper_CI"),
#                         pwOut1[[2]])
#     pwNames <- c("Fst", rownames(pwOut1[[1]]),
#                  "--", "", "Fit", rownames(pwOut1[[1]]))
#     pwOut <- as.matrix(cbind(pwNames, pwWriteOut))
#     dimnames(pwOut) <- NULL
#     
#     if (!is.null(outfile)){
#       of = paste(getwd(), "/", outfile, "-[fstWC]", "/", sep = "")
#       if(is.element("xlsx", installed.packages()[, 1])){
#         # write data to excel
#         # Load dependencies
#         library("xlsx")
#         # standard stats
#         if(bs_locus){
#           write.xlsx(pwOut, file = paste(of, "[fstWC].xlsx", sep = ""),
#                      sheetName = "Pairwise_stats", col.names = FALSE,
#                      row.names = FALSE, append = TRUE)
#         } else {
#           write.xlsx(pwOut, file = paste(of, "[fstWC].xlsx", sep = ""),
#                      sheetName = "Pairwise_stats", col.names = FALSE,
#                      row.names = FALSE, append = FALSE)
#         }                   
#       } else {
#         # text file alternatives
#         std <- file(paste(of, "Pairwise_stats-[fstWC].txt", sep = ""), "w")
#         for(i in 1:nrow(pwOut)){
#           cat(pwOut[i,], "\n", file = std, sep = "\t")
#         }
#         close(std)
#       }
#     }
#     names(pwOut1) <- c("Fst", "Fit")
#   }
#   # return results to the R enviroment
#   if(bs_locus && bs_pairwise){
#     list(locus = locBSout,
#          pairwise = pwOut1)
#   } else if (bs_locus && !bs_pairwise){
#     return(locus = locBSout)
#   } else if (bs_pairwise && !bs_locus){
#     return(pairwise = pwOut1)
#   }
# }
################################################################################
# END
################################################################################