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