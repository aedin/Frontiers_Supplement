###############################################################################
# TOY DATA HISTOGRAMS
###############################################################################
raw <- rnorm(10000, mean=1.5, sd=1.5)

df = list("raw" = raw,
"Scale"=scale(raw, scale=TRUE, center=FALSE),
"Center"= scale(raw, scale=FALSE, center=TRUE),
"Scale+Center" = scale(raw, scale=TRUE, center=TRUE),
"Log2 *" = log2(raw))


plotHist<-function(x, n, xlim=c(-5,5)){
  svg(file=paste0(n, ".svg"), width = 4, height = 3)
  hist(x,main=n, xlab="",col="grey90",freq = FALSE, ylab="", xlim=xlim, yaxt='n',breaks=seq(-20, 20, 0.5))
  abline(v=mean(x,na.rm=TRUE), col="red", lwd=2, lty=2, cex=2)
  #legend("topleft", legend=c(paste0("mean=", round(mean(x, na.rm=TRUE),2)), paste0("stdev=",round(sd(x, na.rm=TRUE),2))), bty = "n")
  dev.off()
}
lapply(seq_along(df), function(i) plotHist(df[[i]], names(df)[i])) 

###############################################################################
# SCMIXOLOGY(10X) HISTOGRAMS
###############################################################################

load('scmix_3lines_intersected.rda')
raw_10x <- counts(sce_sc_10x_qc)

df_10x = list("raw (10x)" = raw_10x,
              "Scale"=scale(raw_10x, scale=TRUE, center=FALSE),
              "Center"= scale(raw_10x, scale=FALSE, center=TRUE),
              "Scale+Center" = scale(raw_10x, scale=TRUE, center=TRUE),
              "Log2 *" = log2(raw_10x))

jpgHist<-function(x, n){
  jpeg(file=paste0(n, "1.jpg"), width = 400, height = 300)
  hist(x,main=n, xlab="",col="grey90",freq = FALSE, ylab="", yaxt='n', cex.axis = 1.5)
  abline(v=mean(x,na.rm=TRUE), col="red", lwd=2, lty=2, cex=2)
  #legend("topleft", legend=c(paste0("mean=", round(mean(x, na.rm=TRUE),2)), paste0("stdev=",round(sd(x, na.rm=TRUE),2))), bty = "n")
  dev.off()
}

for(ind in 1:5){
  jpgHist(df_10x_bigsf[[ind]], names(df_10x_bigsf)[ind])
}
