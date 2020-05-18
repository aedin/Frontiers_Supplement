library("ggplot2")
library("ggpubr")
library("reshape2")

m= rgamma(100,0.5)
#m[sample(1:1000, 200)]<-0
x<-data.frame(raw=m, scaled=scale(m, center=FALSE, scale=TRUE),
              centered=scale(m, center=TRUE, scale=FALSE), 
              centered_scaled=scale(m, center=TRUE, scale=TRUE), 
              log2= log2(m+1))

steps= rbind(names(x)[c(1,2)], names(x)[c(1,3)], names(x)[c(1,4)], names(x)[c(1,5)])
for(i in 1:4) {
  ifile=paste0(paste(steps[i,], collapse = "_"), ".png")
  p<-ggdensity(melt(x[,steps[i,]], variable.name = "data"), x = "value",  ylab="", xlab="", add="mean",
            fill =   "data", facet.by = "data",
            palette = c("#E7B800","#00AFBB"))
            
  ggpar(p,ytickslab=FALSE, legend.title = "")
  ggsave(ifile)
}

# par(mfrow=c(2,3))
# plot(density(m), main="unprocessed")
# plot(density(scale(m, center=TRUE, scale=FALSE)), main="centered")
# plot(density(scale(m, center=TRUE, scale=TRUE)), main="centered and scaled")
# 
# plot(density(log2(m+1)), main="log2 +1",xlab="")
# plot(density(scale(log2(m+1), center=TRUE, scale=FALSE)), main="centered, log2+1")
# plot(density(scale(log2(m+1), center=TRUE, scale=TRUE)), main="centered and scaled, log 2+1")
