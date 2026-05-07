
library(ggplot2)
methods(fortify)
library(latticeExtra)
library(sp)
data(meuse)
class(meuse)
library(gstat)
library(lattice)
# library(maptools)  # removed: package unavailable on CRAN
library(MASS)
# library(contoureR)  # removed: package unavailable on CRAN
library(plotly)
# library(anndata)  # removed: package unavailable on CRAN

# setwd("/Users/aliamiryousefi/Desktop/Tumor boundry/")  # removed: hardcoded local path
all<- (read_h5ad("adata_e24.h5ad"))
tumor<- read_h5ad("adata_e24_tumor.h5ad")

All<- cbind(all$X, all$obs)
T<- cbind(tumor$X, tumor$obs)

plot_slide <- function(patient_ID) {
  # placeholder
}



G<-read.table(file = "unmicst-LSP11347_cell.csv", header = T, sep = ",")
plot(G$X_centroid[log(G$SOX10)<8.2], G$Y_centroid[log(G$SOX10)<8.2], pch=".", col= "Green")
points(G$X_centroid[log(G$SOX10)>8.2], G$Y_centroid[log(G$SOX10)>8.2], pch=".", col= "Red")

#pixel size
ps<- 0.325
#boundary length in micron (between 10-40)
beta<- 30
#obtaining the number of bins
Xbins<- floor(((max(G$X_centroid)-min(G$X_centroid))*ps+1)/beta)
Ybins<- floor(((max(G$Y_centroid)-min(G$Y_centroid))*ps+1)/beta)

Xoffset <- (max(G$X_centroid)-min(G$X_centroid))/40
Yoffset <- (max(G$Y_centroid)-min(G$Y_centroid))/40

x<- G$X_centroid
y<- G$Y_centroid

x.grid<- seq(min(x)-Xoffset, max(x)+Xoffset, (max(x)-min(x)+2*Xoffset)/Xbins)
y.grid<- seq(min(y)-Yoffset, max(y)+Yoffset, (max(y)-min(y)+2*Yoffset)/Ybins)


vol<- matrix(0, Xbins, Ybins)

for (i in 1:Xbins){
  for (j in 1:Ybins){
    sec<- which(G$X_centroid < x.grid[i+1] & G$X_centroid >= x.grid[i] & G$Y_centroid < y.grid[j+1] & G$Y_centroid >= y.grid[j])
    if (length(sec)>0){
      vol[i,j]<- sum(log(G[sec, "SOX10"]))/length(sec)
    }
  }
}

filled.contour(t(vol), plot.axes = {
  axis(1)
  axis(2)
  contour(t(vol), add = TRUE, lwd = 0.5)
}
)

fig <- plot_ly(
       type = 'contour',
       z = vol, line = list(smoothing = 0.5), colorscale = 'Jet',
       contours = list(
           coloring = 'heatmap'
       )
   )
fig


#simple trick to get the inner tumor or stroma borders


#in tumor
biVol<- 1*(vol < 8.2)


#in stroma 
#biVol<- 1*(vol > 8.2)

#it seems the inner tumor part is a bit more sampled so to be checked!

library(ggplot2)
x = 1:nrow(biVol)
y = 1:ncol(biVol)
z = expand.grid(x=x,y=y); z$z = apply(z,1,function(xx){ biVol[ xx[1],xx[2] ]} )
df = getContourLines(z, levels = c(0,1))
ggplot(df,aes(x,y,group=Group,colour=z)) + geom_path()
plot(df$x, df$y, pch=".")

xy<-unique(cbind(df$x,df$y))
xy<- xy+1 



all<- (read_h5ad("adata_e24.h5ad"))
tumor<- read_h5ad("adata_e24_tumor.h5ad")

All<- cbind(all$X, all$obs)
T<- cbind(tumor$X, tumor$obs)

#MART1, SOX9, MITF, 
gMART1<- allgates$gates["MART1","LSP11347"]

gSOX9<- allgates$gates["SOX9","LSP11347"]

gMITF<- allgates$gates["MITF","LSP11347"]

##agregating all the cells on the boundry for 30 micron

inboundTumorcells<- numeric(dim(G)[1])
inboundStromcells<- numeric(dim(G)[1])


#run this only in conjunction with #in tumor
for (k in 1:dim(xy)[1]){
  i<- xy[k,1]
  j<- xy[k,2]
  sec<- which(G$X_centroid < x.grid[i+1] & G$X_centroid >= x.grid[i] & G$Y_centroid < y.grid[j+1] & G$Y_centroid >= y.grid[j])
  inboundTumorcells[sec]<-1
}


#run this only in conjunction with #in stroma
for (k in 1:dim(xy)[1]){
  i<- xy[k,1]
  j<- xy[k,2]
  sec<- which(G$X_centroid < x.grid[i+1] & G$X_centroid >= x.grid[i] & G$Y_centroid < y.grid[j+1] & G$Y_centroid >= y.grid[j])
  inboundStromcells[sec]<-1
}


gMART1<- exp(gMART1)
gSOX9<- exp(gSOX9)
gMITF<- exp(gMITF)
#MSM lineages

ppp<- G$MART1 > gMART1 & G$SOX9 > gSOX9 & G$MITF > gMITF
ppn<- G$MART1 > gMART1 & G$SOX9 > gSOX9 & G$MITF < gMITF
pnp<- G$MART1 > gMART1 & G$SOX9 < gSOX9 & G$MITF > gMITF
pnn<- G$MART1 > gMART1 & G$SOX9 < gSOX9 & G$MITF < gMITF
npp<- G$MART1 < gMART1 & G$SOX9 > gSOX9 & G$MITF > gMITF
npn<- G$MART1 < gMART1 & G$SOX9 > gSOX9 & G$MITF < gMITF
nnp<- G$MART1 < gMART1 & G$SOX9 < gSOX9 & G$MITF > gMITF
nnn<- G$MART1 < gMART1 & G$SOX9 < gSOX9 & G$MITF < gMITF

GC<- cbind(G, ppp, ppn, pnp, pnn, npp, npn, nnp, nnn, inboundStromcells, inboundTumorcells)





GCinbt<- subset(GC, GC$inboundTumorcells==1)
GCint<- subset(GC, GC$inboundTumorcells==0 & log(GC$SOX10) > 8.2 )

GCinbs<- subset(GC, GC$inboundStromcells==1)
GCins<- subset(GC, GC$inboundStromcells==0 & log(GC$SOX10) < 8.2 )


Lineage_porp_inbt<- c(sum(GCinbt$ppp), sum(GCinbt$ppn), sum(GCinbt$pnp), sum(GCinbt$pnn), sum(GCinbt$npp), sum(GCinbt$npn), sum(GCinbt$nnp), sum(GCinbt$nnn))/ dim(GCinbt)[1]

Lineage_porp_int<- c(sum(GCint$ppp), sum(GCint$ppn), sum(GCint$pnp), sum(GCint$pnn), sum(GCint$npp), sum(GCint$npn), sum(GCint$nnp), sum(GCint$nnn))/ dim(GCint)[1]

Lineage_porp_inbs<- c(sum(GCinbs$ppp), sum(GCinbs$ppn), sum(GCinbs$pnp), sum(GCinbs$pnn), sum(GCinbs$npp), sum(GCinbs$npn), sum(GCinbs$nnp), sum(GCinbs$nnn))/ dim(GCinbs)[1]

Lineage_porp_ins<- c(sum(GCins$ppp), sum(GCins$ppn), sum(GCins$pnp), sum(GCins$pnn), sum(GCins$npp), sum(GCins$npn), sum(GCins$nnp), sum(GCins$nnn))/ dim(GCins)[1]


bardata<-rbind(Lineage_porp_int, Lineage_porp_inbt, -Lineage_porp_inbs, -Lineage_porp_ins)
rownames(bardata)<- c("INS-Tumor", "INB-Tumor", "INB-Stroma", "INS-Stroma")
colnames(bardata)<- c("+++", "++-", "+-+", "+--", "-++", "-+-", "--+", "---")
barplot(height = bardata, beside = TRUE, col=c("darkred", "red", "orange", "yellow"), width = 2, ylim = c(-0.6,0.6), las=1, legend.text = rownames(bardata))
title(main = "Tumor Boundary Lineage Profile", font.main = 4, sub="MART1/SOX9/MITF")
#text(30.2,-0.4, "MART1", cex=1.2, col = "blue")

