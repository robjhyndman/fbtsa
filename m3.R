rm(list=ls())

library(Mcomp)
library(gridExtra)
library(ggplot2)
library(GGally)
library(forecast)

savepdf <- function(file, width=16, height=10)
{
  fname <<- paste("figs/",file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54, pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}
endpdf <- function() {crop::dev.off.crop(fname)}

load("M3Features.Rdata")
M3Features <- as.data.frame(M3Features)
colnames(M3Features) <- c("SpectralEntropy","Trend","Seasonality",
                          "Period","ResidACF","Lambda")
M3Features = subset(as.data.frame(M3Features), select = c(4,1,2,3,5,6))
feaNames = c('Period', "Entropy",'Trend','Seasonality','ACF','Lambda')
M3Features$Period= as.factor(M3Features$Period)

# Histograms
gghist <- function(data, mapping, ...)
{
  x <- data[,as.character(mapping$x[2])]
  bw <- 0.2*bw.nrd0(x) + 0.8*bw.SJ(x)
  p <- ggplot(data, mapping) +
    geom_density(col=NA, fill="#cc5900", bw=bw)
  return(p)
}

savepdf("PairwisePlot", height=15,width=16)
yk_ggally_densityDiag <- wrap(gghist, adjust = 0.5)
yk_ggally_barDiag <-  wrap(ggally_barDiag, colour="#cc5900",
                           fill ="#cc5900", width = 0.2)
ggpairs(M3Features,
          diag = list(continuous = yk_ggally_densityDiag,
                      discrete = yk_ggally_barDiag),
          columnLabels = feaNames,
          axisLabels = "none",
          lower=list(continuous = wrap("points", alpha = 0.5,  size=0.2)))
endpdf()

for(i in 1:NCOL(M3Features))
{
  fname <- paste("kde",colnames(M3Features)[i],sep="")
  savepdf(fname)
  if(is.factor(M3Features[,i]))
    print(ggally_barDiag(M3Features,
                   mapping = aes_string(colnames(M3Features)[i]),
                   colour="#cc5900", fill="#cc5900", width=0.2))
  else
    print(gghist(M3Features,
             mapping = aes_string(colnames(M3Features)[i])))
  endpdf()
}

savepdf("histPeriod")
ggally_barDiag(M3Features,
               mapping = aes(Period), width=0.2,
               colour="#cc5900", fill="#cc5900")
endpdf()

savepdf("TrendSE")
ggplot(M3Features, aes(x=SpectralEntropy,y=Trend)) + geom_point()
endpdf()

savepdf("ACF1SE")
ggplot(M3Features, aes(x=SpectralEntropy,y=ResidACF)) + geom_point()
endpdf()



#Consider only long series
n <- unlist(lapply(M3,function(x){x$n}))
M3Featureslong <- M3Features[n>50,]
M3long <- M3[names(M3)[n>50]]

fnames <- c("M3Freq","M3spec","M3trend","M3season","M3acf","M3lambda")
titles <- c("Period","Spectral Entropy","Trend","Seasonality","ACF1","Box Cox")

k <- NROW(M3Featureslong)

for(i in 1:6)
{
  j <- order(M3Featureslong[,i])
  savepdf(paste(fnames[i],"Lo",sep=""), width=20, height=7)
  print(autoplot(M3long[[j[1]]]$x) +
    ylab(M3long[[j[1]]]$sn) + xlab(""))
  endpdf()
  savepdf(paste(fnames[i],"Hi",sep=""), width=20, height=7)
  print(autoplot(M3long[[j[k]]]$x) +
    ylab(M3long[[j[k]]]$sn) + xlab(""))
  endpdf()
}


x <- M3long[["N2096"]]
plot(x)
fit <- stl(x$x, s.window=11)

savepdf("N2096")
autoplot(fit) + xlab("") +
  scale_x_continuous(breaks=seq(1982,1992,by=1), minor_breaks = NULL)
endpdf()


## INSTANCE SPACE

load("M3Features.Rdata")
feaNames = c("Entropy",'Trend','Seasonality','Period','ACF','Lambda')
scaledM3Features = scale(as.data.frame(M3Features))
M3.pc.cr <- princomp((scaledM3Features))
DataMean = colMeans(M3Features)
DataSD = apply(M3Features,2,sd)
PC1 <- M3.pc.cr$scores[, 1]
PC2 <- M3.pc.cr$scores[, 2]
targetcol = '#D95F02'
targetsize = 10
egs = c(2963, 970,139,626,815,731, 781, 882, 2786,2794,2612,2442)
egspc = data.frame(pc1 = PC1[egs], pc2 = PC2[egs])

# plot the instance space
require(ggplot2)
savepdf("InstanceSpace", height=10,width=14)
ggplot() +
  geom_point(data = data.frame(x = PC1, y = PC2),
             aes(x = PC1, y = PC2), colour = "#939393", shape = 19, size = 2) +
  ggtitle("First two PCs explain 60% of variation")
endpdf()

mycols = c("#1B9E77", "#FE4902")
myggplotFig3 <- function(d, i) {
  ggplot(data = d, mapping = aes(x = PC1, y = PC2)) +
      geom_point(aes(colour = value),  size = 2, shape=19) +
      ggtitle(feaNames[i]) + scale_colour_gradientn(colours = mycols)
}
for(i in 1:6)
{
  savepdf(paste("Feature",i,sep=""))
  print(myggplotFig3(d = data.frame(x = PC1, y = PC2,
                            value = M3Features[, i]),i))
  endpdf()
}

## MIN MASE

load("M3MASEout.Rdata")
#library(splus2R)
minMASE = apply(M3MASEout,1,min)
d = data.frame(x=PC1,y=PC2)
a = quantile(minMASE,probs=c(0.2,0.5))[1]
b = quantile(minMASE,probs=c(0.2,0.5))[2]
d1 = data.frame(x=PC1[minMASE<a],y=PC2[minMASE<a])

savepdf("MASElo", height=10,width=14)
ggplot() +
  geom_point(data=d,aes(x=PC1,y=PC2),colour = '#939393', shape=19,size=2)+
  geom_point(data=d1,aes(x=x,y=y),colour ='#1B9E77', shape=19,size=2)+
  ggtitle('Low')
endpdf()

savepdf("MASEmid", height=10,width=14)
d2 = data.frame(x=PC1[(minMASE>=a)&(minMASE<=b)],y=PC2[(minMASE>=a)&(minMASE<=b)])
ggplot() +
  geom_point(data=d,aes(x=PC1,y=PC2),colour = '#939393', shape=19,size=2)+
  geom_point(data=d2,aes(x=x,y=y),colour ='#D95F02', shape=19,size=2)+
  ggtitle('Middle')
endpdf()

savepdf("MASEhi", height=10,width=14)
d3 = data.frame(x=PC1[minMASE>b],y=PC2[minMASE>b])
ggplot() +
  geom_point(data=d,aes(x=PC1,y=PC2),colour = '#939393', shape=19,size=2)+
  geom_point(data=d3,aes(x=x,y=y),colour ="#7570B3", shape=19,size=2)+
  ggtitle('High')
endpdf()


# MASE RAINBOW

titles = c('Naive','Seasonal naive', 'Theta', 'ETS','ARIMA','STL-AR')

myggplotFig8 <- function(d, i) {
  ggplot(data = d, mapping = aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = MASE), size = 2) +
    ggtitle(titles[i]) +
    scale_colour_gradientn(colours = mycols, limits = c(0, 4))
}

for(i in 1:6)
{
  savepdf(paste("MASE",i,sep=""))
  print(myggplotFig8(
    d = data.frame(PC1 = PC1,PC2 = PC2,MASE = pmin(4,M3MASEout[,i])), i))
  endpdf()
}

d = data.frame(PC1 = PC1,PC2 = PC2,MASE = M3MASEout[,2])
d$MASE[d$MASE > 4] = 4
p2 <- myggplotFig8(d, 2) + theme(axis.title = element_text(color = 'white'))

d = data.frame(PC1 = PC1,PC2 = PC2,MASE = M3MASEout[,3])
d$MASE[d$MASE > 4] = 4
p3 <- myggplotFig8(d, 3) + theme(axis.title.x = element_text(color = 'white'))

d = data.frame(PC1 = PC1,PC2 = PC2,MASE = M3MASEout[,4])
d$MASE[d$MASE > 4] = 4
p4 <- myggplotFig8(d, 4) + theme(axis.title = element_text(color = 'white'))

d = data.frame(PC1 = PC1,PC2 = PC2,MASE = M3MASEout[,5])
d$MASE[d$MASE > 4] = 4
p5 <- myggplotFig8(d, 5)

d = data.frame(PC1 = PC1,PC2 = PC2,MASE = M3MASEout[,6])
d$MASE[d$MASE > 4] = 4
p6 <- myggplotFig8(d, 6) + theme(axis.title.y = element_text(color = 'white'))

grid.arrange(p1, p2, p3, p4, p5, p6,ncol = 2)


