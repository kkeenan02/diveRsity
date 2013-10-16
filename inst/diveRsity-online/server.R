# Load packages
library("shiny")
library("plotrix")#, lib.loc = "/home/kkeenan/depends/")
library("diveRsity")#, lib.loc = "/home/kkeenan/depends/")
#library("shinyIncubator")

shinyServer(function(input, output, session) {

  # calculate std and est stats
  stdest <- reactive ({
    if(input$goButton==0) return(NULL)
    isolate({
      infile <- input$file$datapath
      divPart(infile = infile,
              outfile = NULL,
              gp = input$gp,
              pairwise = FALSE,
              WC_Fst = input$WC_Fst,
              bs_locus = FALSE,
              bs_pairwise = FALSE,
              bootstraps = FALSE,
              plot = FALSE,
              parallel = FALSE)
    })
  })
  
  # Calculate pairwise matrix
  pwOut <- reactive ({
    if(input$goButton==0) return(NULL)
    isolate({
      infile <- input$file$datapath
      divPart(infile = infile,
              outfile = NULL,
              gp = input$gp,
              pairwise = input$pairwise,
              WC_Fst = input$WC_Fst,
              bs_locus = FALSE,
              bs_pairwise = FALSE,
              bootstraps = FALSE,
              plot = FALSE,
              parallel = input$parallel)
    })
  })
  
  # Calculate locus bs
  lbsOut <- reactive ({
    if(input$goButton==0) return(NULL)
    isolate({
      infile <- input$file$datapath
      divPart(infile = infile,
              outfile = NULL,
              gp = input$gp,
              pairwise = FALSE,
              WC_Fst = input$WC_Fst,
              bs_locus = input$bs_locus,
              bs_pairwise = FALSE,
              bootstraps = input$bootstraps,
              plot = FALSE,
              parallel = input$parallel)    
    })
  })
  
  # Calculate pairwise bs
  pwbsOut <- reactive ({
    if(input$goButton==0) return(NULL)
    isolate({
      infile <- input$file$datapath
      divPart(infile = infile,
              outfile = NULL,
              gp = input$gp,
              pairwise = FALSE,
              WC_Fst = input$WC_Fst,
              bs_locus = FALSE,
              bs_pairwise = input$bs_pairwise,
              bootstraps = input$bootstraps,
              plot = FALSE,
              parallel = input$parallel)
    })
  })
  
  
  
  divBout <- reactive({
    if(input$goButton==0) return(NULL)
    isolate({
      if(!is.null(input$divBasic)){
        divBasic(infile = input$file$datapath,
                 outfile = NULL,
                 gp = input$gp)
      }
    })  
  })
  
  #############################################################################
  # divBasic output
  #############################################################################
  output$divB <- renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      if(input$divBasic && !is.null(input$file)){
        res <- divBout()
        return(as.data.frame(res$mainTab))
      }
    })    
  })
  
  #############################################################################
  # Download divBasic results
  #############################################################################
  output$divBdl <- downloadHandler(
    filenames <- function(file){
      paste("divBasic", Sys.Date(), "-[diveRsity-online].txt", sep = "")
    },
    content <- function(file){
      res <- divBout()
      outer <- res$mainTab
      write.table(outer, file, append = FALSE, quote = FALSE,
                  sep = "\t", eol = "\r\n", row.names = FALSE,
                  col.names = FALSE)
    }
  )
                  
  
  
  #############################################################################
  # Standard stats
  #############################################################################  
  output$std <-  renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      if(!is.null(input$file)) {
        out <- stdest()
        return(as.data.frame(out$standard))
      }    
    })
  })
  
  
  
  #Download standard data
  output$dlstd <- downloadHandler(
    #if(!is.null(input$file)) {
      filename <- function() {
        paste("standard_", Sys.Date(), "_[diveRsity-online].txt", sep = "")
      },
      content <- function(file) {
        out <- stdest()
        prestd <- out$standard
        std <- cbind(rownames(prestd), prestd)
        colnames(std) <- c("Loci", colnames(prestd))
        write.table(std, file, append = FALSE, quote = FALSE,
                    sep = "\t", eol = "\r\n", row.names = FALSE)
      }
    #} else {
     # print("No file specified!")
    #}
  )   
  #############################################################################
  # Estimated stats
  #############################################################################
  output$est <- renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      if(!is.null(input$file)) {
        out <- stdest()
        return(as.data.frame(out$estimate))
      }
    })
  })
  
  #Download standard data
  output$dlest <- downloadHandler(
    #if(!is.null(input$file)) {
      filename <- function() {
        paste("estimate_", Sys.Date(), "_[diveRsity-online].txt", sep = "")
      },
      content <- function(file) {
        out <- stdest()
        preest <- out$estimate
        est <- cbind(rownames(preest), preest)
        colnames(est) <- c("Loci", colnames(preest))
        write.table(est, file, append = FALSE, quote = FALSE,
                    sep = "\t", eol = "\r\n", row.names = FALSE)
      }
  )
  #############################################################################
  # Pairwise matrices
  #############################################################################
  output$pw <- renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      out <- pwOut()
      if(is.element("pairwise", names(out))){
        pw_fix <- lapply(out$pairwise, function(x){
          matrix(x, ncol = ncol(x), nrow = nrow(x))
        })
        for(i in 1:length(pw_fix)){
          pw_fix[[i]][is.na(pw_fix[[i]])] <- ""
        }
        spltr <- matrix(rep("", (ncol(pw_fix[[1]]))+1), nrow = 1, 
                        ncol = (ncol(pw_fix[[1]])+1))
        rownames(spltr) <- NULL
        rowcol <- c("",colnames(out$pairwise[[1]]))
        dimnames(rowcol) <- NULL
        spltr_nm <- matrix(c("Gst_est", rep("", (length(spltr)-1))), 
                           ncol = length(spltr), nrow = 1)
        rownames(spltr_nm) <- NULL
        pre_pw <- rbind(rowcol[-1], pw_fix[[4]])
        pw <- rbind(spltr_nm, cbind(rowcol, pre_pw))
        if(!input$WC_Fst){
          for(i in 5:6){
            spltr_nm <- matrix(c(names(out$pairwise)[i], 
                                 rep("", (length(spltr)-1))), 
                               ncol = length(spltr), nrow = 1)
            rownames(spltr_nm) <- NULL
            pre_pw <- rbind(rowcol[-1], pw_fix[[i]])
            pw <- rbind(pw, spltr, spltr_nm, cbind(rowcol, pre_pw))
          }
        } else {
          for (i in c(5,6,7,8)){
            spltr_nm <- matrix(c(names(out$pairwise)[i], 
                                 rep("", (length(spltr)-1))), 
                               ncol = length(spltr), nrow = 1)
            rownames(spltr_nm) <- NULL
            pre_pw <- rbind(rowcol[-1], pw_fix[[i]])
            pw <- rbind(pw, spltr, spltr_nm, cbind(rowcol, pre_pw))
          }
        }
        dimnames(pw) <- NULL
        return(pw)
      }
    })  
  })
  
  #Download pairwise matrix data
  output$dlpw <- downloadHandler(
    filename <- function() {
      paste("pairwise_matrix_", Sys.Date(), "_[diveRsity-online].txt", 
            sep = "")
    },
    content <- function(file) {
      out <- pwOut()
      pw_fix <- lapply(out$pairwise, function(x){
        matrix(x, ncol = ncol(x), nrow = nrow(x))
      })
      for(i in 1:length(pw_fix)){
        pw_fix[[i]][is.na(pw_fix[[i]])] <- ""
      }
      spltr <- matrix(rep("", (ncol(pw_fix[[1]]))+1), nrow = 1, 
                      ncol = (ncol(pw_fix[[1]])+1))
      rownames(spltr) <- NULL
      rowcol <- c("",colnames(out$pairwise[[1]]))
      dimnames(rowcol) <- NULL
      spltr_nm <- matrix(c("Gst_est", rep("", (length(spltr)-1))), 
                         ncol = length(spltr), nrow = 1)
      rownames(spltr_nm) <- NULL
      pre_pw <- rbind(rowcol[-1], pw_fix[[4]])
      pw <- rbind(spltr_nm, cbind(rowcol, pre_pw))
      if(!input$WC_Fst){
        for(i in 5:6){
          spltr_nm <- matrix(c(names(out$pairwise)[i], 
                               rep("", (length(spltr)-1))), 
                             ncol = length(spltr), nrow = 1)
          rownames(spltr_nm) <- NULL
          pre_pw <- rbind(rowcol[-1], pw_fix[[i]])
          pw <- rbind(pw, spltr, spltr_nm, cbind(rowcol, pre_pw))
        }
      } else {
        for (i in c(5,6,7,8)){
          spltr_nm <- matrix(c(names(out$pairwise)[i], 
                               rep("", (length(spltr)-1))), 
                             ncol = length(spltr), nrow = 1)
          rownames(spltr_nm) <- NULL
          pre_pw <- rbind(rowcol[-1], pw_fix[[i]])
          pw <- rbind(pw, spltr, spltr_nm, cbind(rowcol, pre_pw))
        }
      }
      dimnames(pw) <- NULL  
      write.table(pw, file, append = FALSE, quote = FALSE,
                  sep = "\t", eol = "\r\n", row.names = FALSE,
                  col.names = FALSE)
    }
  )
  #############################################################################
  # Locus bootstraps
  #############################################################################
  output$bs_loc <- renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      if(input$bs_locus == TRUE){
        out <- lbsOut()
        splt <- c("--","--","--","--")
        rownames(splt) <-  NULL
        splt_nm <- c("Gst_est","","","")
        rownames(splt_nm) <- NULL
        bs_loc <- rbind(splt_nm, cbind(rownames(out$bs_locus$Gst_est),
                                       out$bs_locus$Gst_est))
        if(!input$WC_Fst){
          for (i in 5:6){
            splt_nm <- c(names(out$bs_locus)[i],"","","")
            adder <- cbind(rownames(out$bs_locus[[i]]),
                           out$bs_locus[[i]])
            suppressWarnings(bs_loc <- rbind(bs_loc, splt, splt_nm, adder))
          }
        } else {
          for (i in c(5,6,7,8)){
            splt_nm <- c(names(out$bs_locus)[i],"","","")
            adder <- cbind(rownames(out$bs_locus[[i]]),
                           out$bs_locus[[i]])
            suppressWarnings(bs_loc <- rbind(bs_loc, splt, splt_nm, adder))
          }
        }
        rownames(bs_loc) <- NULL
        colnames(bs_loc) <- c("Loci", "Actual", "Lower", "Upper")
        return(bs_loc)
      }
    })
  })
  
  #Download bs _pw data
  output$dllcbs <- downloadHandler(
    filename <- function() {
      paste("Locus_bootstrap_", Sys.Date(), "_[diveRsity-online].txt", sep = "")
    },
    content <- function(file) {
      if(input$bs_locus == TRUE){
        out <- lbsOut()
        splt <- c("--","--","--","--")
        rownames(splt) <-  NULL
        splt_nm <- c("Gst_est","","","")
        rownames(splt_nm) <- NULL
        bs_loc <- rbind(splt_nm, cbind(rownames(out$bs_locus$Gst_est),
                                       out$bs_locus$Gst_est))
        if(!input$WC_Fst){
          for (i in 5:6){
            splt_nm <- c(names(out$bs_locus)[i],"","","")
            adder <- cbind(rownames(out$bs_locus[[i]]),
                           out$bs_locus[[i]])
            suppressWarnings(bs_loc <- rbind(bs_loc, splt, splt_nm, adder))
          }
        } else {
          for (i in c(5,6,7,8)){
            splt_nm <- c(names(out$bs_locus)[i],"","","")
            adder <- cbind(rownames(out$bs_locus[[i]]),
                           out$bs_locus[[i]])
            suppressWarnings(bs_loc <- rbind(bs_loc, splt, splt_nm, adder))
          }
        }
        rownames(bs_loc) <- NULL
        colnames(bs_loc) <- c("Loci", "Actual", "Lower", "Upper")
      }
      write.table(bs_loc, file, append = FALSE, quote = FALSE,
                  sep = "\t", eol = "\r\n", row.names = FALSE)
    }
  )
  
  ############################################################################
  # Pairwise bootstrap
  ############################################################################  
  output$pw_bs <- renderTable({
    if(input$goButton==0) return(NULL)
    isolate({
      if(input$bs_pairwise == TRUE){
        out <- pwbsOut()
        splt <- c("--","--","--","--")
        splt_nm <- c("Gst_est", "","", "")
        pw <- rbind(splt_nm, cbind(rownames(out$bs_pairwise$Gst_est),
                                   out$bs_pairwise$Gst_est))
        if(!input$WC_Fst){
          for(i in 5:6){
            splt_nm <- c(names(out$bs_pairwise)[i], "", "", "")
            adder <- cbind(rownames(out$bs_pairwise[[i]]),
                           out$bs_pairwise[[i]])
            suppressWarnings(pw <- rbind(pw, splt, splt_nm, adder))                                               
          }
        } else {
          for(i in c(5,6,7,8)){
            splt_nm <- c(names(out$bs_pairwise)[i], "", "", "")
            adder <- cbind(rownames(out$bs_pairwise[[i]]),
                           out$bs_pairwise[[i]])
            suppressWarnings(pw <- rbind(pw, splt, splt_nm, adder))
          }
        }
        rownames(pw) <- NULL
        colnames(pw) <- c("POPS", "Actual", "Lower", "Upper")
        return(pw)    
      }
    })
  })
  
  #Download bs _pw data
  output$dlpwbs <- downloadHandler(
    filename <- function() {
      paste("Pairwise_bootstrap_", Sys.Date(), "_[diveRsity-online].txt",
            sep = "")
    },
    content <- function(file) {
      if(input$bs_pairwise == TRUE){
        out <- pwbsOut()
        splt <- c("--","--","--","--")
        splt_nm <- c("Gst_est", "","", "")
        pw <- rbind(splt_nm, cbind(rownames(out$bs_pairwise$Gst_est),
                                   out$bs_pairwise$Gst_est))
        if(!input$WC_Fst){
          for(i in 5:6){
            splt_nm <- c(names(out$bs_pairwise)[i], "", "", "")
            adder <- cbind(rownames(out$bs_pairwise[[i]]),
                           out$bs_pairwise[[i]])
            suppressWarnings(pw <- rbind(pw, splt, splt_nm, adder))                                               
          }
        } else {
          for(i in c(5,6,7,8)){
            splt_nm <- c(names(out$bs_pairwise)[i], "", "", "")
            adder <- cbind(rownames(out$bs_pairwise[[i]]),
                           out$bs_pairwise[[i]])
            suppressWarnings(pw <- rbind(pw, splt, splt_nm, adder))
          }
        }
        rownames(pw) <- NULL
        colnames(pw) <- c("POPS", "Actual", "Lower", "Upper")    
      }
      
      write.table(pw, file, append = FALSE, quote = FALSE,
                  sep = "\t", eol = "\r\n", row.names = FALSE)
    }
  )
  
  #plot attempt
  output$cor <- renderPlot({
    if(input$goButton==0) return(NULL)
    isolate({
      if(input$corplot == TRUE){
        infile <- input$file$datapath
        x <- readGenepop(infile, input$gp, FALSE)
        y <- divPart(infile = infile,
                     outfile = NULL,
                     gp = input$gp,
                     pairwise = FALSE,
                     WC_Fst = TRUE,
                     bs_locus = FALSE,
                     bs_pairwise = FALSE,
                     bootstraps = 0,
                     plot = FALSE,
                     parallel = FALSE)
        par(mfrow = c(2, 2))
        par(mar = c(4, 5, 2, 2))
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
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 8] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression(hat(theta)), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 8] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor1 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 8], x[[16]])
        sig <- sigStar(cor1)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor1$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 4] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression(G[st]), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 4] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor2 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 4], x[[16]])
        sig <- sigStar(cor2)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor2$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 5] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression("G'"[st]), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 5] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor3 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 5], x[[16]])
        sig <- sigStar(cor3)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor3$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 6] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression(D[est]), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 6] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor4 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 6], x[[16]])
        sig <- sigStar(cor4)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor4$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
      }   
    })
  })
  output$corplt <- downloadHandler(
    filename = function() {
      paste("corPlot_", Sys.Date(), "_[diveRsity-online].pdf",
            sep = "")
    },
    content = function(file) {
      temp <- tempfile()
      on.exit(unlink(temp))
      if(input$corplot == TRUE){
        if(is.null(input$file)) {
          infile <- "./Test_data.txt"
        } else {
          infile <- input$file$datapath
        }
        x <- readGenepop(infile, input$gp, FALSE)
        y <- divPart(infile = infile,
                     outfile = NULL,
                     gp = input$gp,
                     pairwise = FALSE,
                     WC_Fst = TRUE,
                     bs_locus = FALSE,
                     bs_pairwise = FALSE,
                     bootstraps = 0,
                     plot = FALSE,
                     parallel = FALSE)
        par(mfrow = c(2, 2))
        par(mar = c(4, 5, 2, 2))
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
        pdf(file = temp)
        par(mfrow = c(2,2), mar = c(5,5,2,2))
        #par(mfrow = c(2,2))
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 8] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression(hat(theta)), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 8] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor1 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 8], x[[16]])
        sig <- sigStar(cor1)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor1$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 4] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression(G[st]), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 4] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor2 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 4], x[[16]])
        sig <- sigStar(cor2)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor2$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 5] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression("G'"[st]), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 5] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor3 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 5], x[[16]])
        sig <- sigStar(cor3)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor3$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        plot(y[[2]][1:(nrow(y[[2]]) - 1), 6] ~ x[[16]], pch = 16, 
             xlab = "Number of alleles", ylab = expression(D[est]), 
             ylim = c(0, 1), las = 1, cex.lab = 1.5)
        abline(lm(y[[2]][1:(nrow(y[[2]]) - 1), 6] ~ x[[16]]), col = "red", 
               lwd = 2)
        cor4 <- cor.test(y[[2]][1:(nrow(y[[2]]) - 1), 6], x[[16]])
        sig <- sigStar(cor4)
        text(x = max(x[[16]])/1.5, y = 0.8, 
             labels = paste("r = ", round(cor4$estimate[[1]], 3), " ", sig, 
                            sep = ""), cex = 2)
        dev.off()
        bytes <- readBin(temp, "raw", file.info(temp)$size)
        writeBin(bytes, file)
      }
    }
  )
  
  
})
