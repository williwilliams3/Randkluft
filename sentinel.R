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

set.seed(1)

# setwd("~/Desktop/sentinel")  # removed: hardcoded local path

all<- (read_h5ad("adata_e24.h5ad"))
tumor<- read_h5ad("adata_e24_tumor_v3.h5ad")

All<- cbind(all$X, all$obs)
T<- cbind(tumor$X, tumor$obs)

Tumor<- numeric(dim(All)[1])
lineage2<- character(dim(All)[1])
proliferation2<- character(dim(All)[1])
prolif2_lineage2<- character(dim(All)[1])
SOX9_prolif<- character(dim(All)[1])

All<- cbind(All, Tumor, lineage2, proliferation2, prolif2_lineage2, SOX9_prolif)

n_patients<- length(as.character(unique(All$imageid)))

for (i in 1:n_patients){
  image<- as.character(unique(All$imageid))[i]
  Ac<- subset(All, All$imageid==image)
  Tc<- subset(T, T$imageid==image)
  Ac$Tumor[which(Ac$CellID %in% Tc$CellID)]<- 1
  Ac$lineage2[which(Ac$CellID %in% Tc$CellID)] <- as.character(Tc$lineage2)
  Ac$proliferation2[which(Ac$CellID %in% Tc$CellID)] <- as.character(Tc$proliferation2)
  Ac$prolif2_lineage2[which(Ac$CellID %in% Tc$CellID)] <- as.character(Tc$prolif2_lineage2)
  Ac$SOX9_prolif[which(Ac$CellID %in% Tc$CellID)] <- as.character(Tc$SOX9_prolif)
  All[which(All$imageid==image), ] <- Ac
}

# > as.character(unique(All$imageid))
# [1] "LSP11563" "LSP11523" #"LSP11587" "LSP11747" "LSP11411" "LSP11691" "LSP11483" "LSP11595" "LSP11627" "LSP11315" "LSP11419" "LSP11355" 
#"LSP11387" "LSP11403" "LSP11339" "LSP11667" "LSP11347"
# [18] "LSP11395" "LSP11467" "LSP12422" "LSP11643" "LSP11515"
#for a each imageid, forming a Tumor vector

image<- "LSP11691"
G<- subset(All, All$imageid==image)
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


vol <- matrix(0, Xbins, Ybins)

for (i in 1:Xbins){
  for (j in 1:Ybins){
    sec<- which(G$X_centroid < x.grid[i+1] & G$X_centroid >= x.grid[i] & G$Y_centroid < y.grid[j+1] & G$Y_centroid >= y.grid[j])
    if (length(sec)>0){
      vol[i,j]<- sum(G[sec, "SOX10"])/length(sec)
    }
  }
}

fig <- plot_ly(
  type = 'contour',
  z = vol, line = list(smoothing = 0.9), colorscale = 'Jet',
  contours = list(
    coloring = 'heatmap'
  )
)
fig

filled.contour(t(vol), plot.axes = {
  axis(1)
  axis(2)
  contour(t(vol), add = TRUE, lwd = 0.5)
}
)

maxcap = max(vol)

print(maxcap)

biVol<- 1*(vol < 0.5)

x = 1:nrow(biVol)
y = 1:ncol(biVol)
z = expand.grid(x=x,y=y); z$z = apply(z,1,function(xx){ biVol[ xx[1],xx[2] ]} )

#new method based on grDevices
cl<- contourLines(x=1:nrow(biVol), y= 1:ncol(biVol), biVol, levels =c(0,1))

# Convert the list of lists to a dataframe
dfn <- do.call(rbind, lapply(cl, data.frame))

xy<-unique(cbind(dfn$x,dfn$y))
xy<- xy+1 

inboundTumorcells<- numeric(dim(G)[1])

#run this only in conjunction with #in tumor
for (k in 1:dim(xy)[1]){
  i<- xy[k,1]
  j<- xy[k,2]
  sec<- which(G$X_centroid < x.grid[i+1] & G$X_centroid >= x.grid[i] & G$Y_centroid < y.grid[j+1] & G$Y_centroid >= y.grid[j])
  inboundTumorcells[sec]<-1
}

biVol<- 1*(vol > 0.5)

x = 1:nrow(biVol)
y = 1:ncol(biVol)
z = expand.grid(x=x,y=y); z$z = apply(z,1,function(xx){ biVol[ xx[1],xx[2] ]} )

cl<- contourLines(x=1:nrow(biVol), y= 1:ncol(biVol), biVol, levels =c(0,1))

# Convert the list of lists to a dataframe
dfn <- do.call(rbind, lapply(cl, data.frame))

xy<-unique(cbind(dfn$x,dfn$y))
xy<- xy+1 

inboundStromcells<- numeric(dim(G)[1])
#run this only in conjunction with #in stroma
for (k in 1:dim(xy)[1]){
  i<- xy[k,1]
  j<- xy[k,2]
  sec<- which(G$X_centroid < x.grid[i+1] & G$X_centroid >= x.grid[i] & G$Y_centroid < y.grid[j+1] & G$Y_centroid >= y.grid[j])
  inboundStromcells[sec]<-1
}

GC<- cbind(G, inboundTumorcells, inboundStromcells)

par(mfrow = c(2,3))
######1
lineages<- unique(levels(T$lineage2))
l<- length(lineages)
Lineage_porp_inbt<- numeric(l)  
Lineage_porp_int<-  numeric(l)  
Lineage_porp_inbs<- numeric(l)   
Lineage_porp_ins<-  numeric(l)   
for (i in 1:l){
  Lineage_porp_inbt[i]<- sum(GC$inboundTumorcells==1 & GC$lineage2 == lineages[i])/ sum(GC$inboundTumorcells==1 & GC$lineage2 != "")
  Lineage_porp_int[i]<-  sum(GC$inboundTumorcells==0 & GC$Tumor==1 & GC$lineage2 == lineages[i])/ sum(GC$inboundTumorcells==0 & GC$Tumor ==1 & GC$lineage2 != "")  
  Lineage_porp_inbs[i]<- sum(GC$inboundStromcells==1 & GC$lineage2 == lineages[i])/ sum(GC$inboundStromcells==1 & GC$lineage2 != "")   
  Lineage_porp_ins[i]<-  sum(GC$inboundStromcells==0 & GC$Tumor==0 & GC$lineage2 == lineages[i])/ sum(GC$inboundStromcells==0 & GC$Tumor ==0 & GC$lineage2 != "")    
}
bardata<-rbind(Lineage_porp_int, Lineage_porp_inbt, -Lineage_porp_inbs, -Lineage_porp_ins)
rownames(bardata)<- c("Tumor", "Tumor border", "Stroma border", "Stroma")
colnames(bardata)<- c("+++", "++-", "+-+", "+--", "-++", "-+-", "--+", "---")
barplot(height = bardata, beside = TRUE, col=c("darkred", "red", "orange", "yellow"), width = 2, ylim = c(-1,1), las=1)
title(main = paste("Lineage Profile", image, " "), font.main = 4, sub="MART1/MITF/SOX9")

######2
lineages<- unique(levels(T$ROI_major_category))
l<- length(lineages)
Lineage_porp_inbt<- numeric(l)  
Lineage_porp_int<-  numeric(l)  
Lineage_porp_inbs<- numeric(l)   
Lineage_porp_ins<-  numeric(l)   
for (i in 1:l){
  Lineage_porp_inbt[i]<- sum(GC$inboundTumorcells==1 & GC$ROI_major_category == lineages[i])/ sum(GC$inboundTumorcells==1 & GC$ROI_major_category != "")
  Lineage_porp_int[i]<-  sum(GC$inboundTumorcells==0 & GC$Tumor==1 & GC$ROI_major_category == lineages[i])/ sum(GC$inboundTumorcells==0 & GC$Tumor ==1 & GC$ROI_major_category != "")  
  Lineage_porp_inbs[i]<- sum(GC$inboundStromcells==1 & GC$ROI_major_category == lineages[i])/ sum(GC$inboundStromcells==1 & GC$ROI_major_category != "")   
  Lineage_porp_ins[i]<-  sum(GC$inboundStromcells==0 & GC$Tumor==0 & GC$ROI_major_category == lineages[i])/ sum(GC$inboundStromcells==0 & GC$Tumor ==0 & GC$ROI_major_category != "")    
}
bardata<-rbind(Lineage_porp_int, Lineage_porp_inbt, -Lineage_porp_inbs, -Lineage_porp_ins)
rownames(bardata)<- c("Tumor", "Tumor border", "Stroma border", "Stroma")
lineages[5]<- "comp-reg"
colnames(bardata)<- lineages
barplot(height = bardata, beside = TRUE, col=c("darkred", "red", "orange", "yellow"), width = 2, ylim = c(-1,1), las=2)
title(main = paste("Major Category"), font.main = 4, sub="")

######3
lineages<- unique(levels(T$ROI_minor_category))
l<- length(lineages)
Lineage_porp_inbt<- numeric(l)  
Lineage_porp_int<-  numeric(l)  
Lineage_porp_inbs<- numeric(l)   
Lineage_porp_ins<-  numeric(l)   
for (i in 1:l){
  Lineage_porp_inbt[i]<- sum(GC$inboundTumorcells==1 & GC$ROI_minor_category == lineages[i])/ sum(GC$inboundTumorcells==1 & GC$ROI_minor_category != "")
  Lineage_porp_int[i]<-  sum(GC$inboundTumorcells==0 & GC$Tumor==1 & GC$ROI_minor_category == lineages[i])/ sum(GC$inboundTumorcells==0 & GC$Tumor ==1 & GC$ROI_minor_category != "")  
  Lineage_porp_inbs[i]<- sum(GC$inboundStromcells==1 & GC$ROI_minor_category == lineages[i])/ sum(GC$inboundStromcells==1 & GC$ROI_minor_category != "")   
  Lineage_porp_ins[i]<-  sum(GC$inboundStromcells==0 & GC$Tumor==0 & GC$ROI_minor_category == lineages[i])/ sum(GC$inboundStromcells==0 & GC$Tumor ==0 & GC$ROI_minor_category != "")    
}
bardata<-rbind(Lineage_porp_int, Lineage_porp_inbt, -Lineage_porp_inbs, -Lineage_porp_ins)
rownames(bardata)<- c("Tumor", "Tumor border", "Stroma border", "Stroma")
colnames(bardata)<- lineages
barplot(height = bardata, beside = TRUE, col=c("darkred", "red", "orange", "yellow"), width = 2, ylim = c(-1,1), las=2)
title(main = paste("minor Category"), font.main = 4, sub="")

####4
lineages<- unique(levels(T$proliferation2))
l<- length(lineages)
Lineage_porp_inbt<- numeric(l)  
Lineage_porp_int<-  numeric(l)  
Lineage_porp_inbs<- numeric(l)   
Lineage_porp_ins<-  numeric(l)   
for (i in 1:l){
  Lineage_porp_inbt[i]<- sum(GC$inboundTumorcells==1 & GC$proliferation2 == lineages[i])/ sum(GC$inboundTumorcells==1 & GC$proliferation2 != "")
  Lineage_porp_int[i]<-  sum(GC$inboundTumorcells==0 & GC$Tumor==1 & GC$proliferation2 == lineages[i])/ sum(GC$inboundTumorcells==0 & GC$Tumor ==1 & GC$proliferation2 != "")  
  Lineage_porp_inbs[i]<- sum(GC$inboundStromcells==1 & GC$proliferation2 == lineages[i])/ sum(GC$inboundStromcells==1 & GC$proliferation2 != "")   
  Lineage_porp_ins[i]<-  sum(GC$inboundStromcells==0 & GC$Tumor==0 & GC$proliferation2 == lineages[i])/ sum(GC$inboundStromcells==0 & GC$Tumor ==0 & GC$proliferation2 != "")    
}
bardata<-rbind(Lineage_porp_int, Lineage_porp_inbt, -Lineage_porp_inbs, -Lineage_porp_ins)
rownames(bardata)<- c("Tumor", "Tumor border", "Stroma border", "Stroma")
colnames(bardata)<- c("nP", "P")
barplot(height = bardata, beside = TRUE, col=c("darkred", "red", "orange", "yellow"), width = 2, ylim = c(-1,1), las=2)
title(main = paste("Proliferation"), font.main = 4, sub="")

####5
lineages<- unique(levels(T$SOX9_prolif))
l<- length(lineages)
Lineage_porp_inbt<- numeric(l)  
Lineage_porp_int<-  numeric(l)  
Lineage_porp_inbs<- numeric(l)   
Lineage_porp_ins<-  numeric(l)   
for (i in 1:l){
  Lineage_porp_inbt[i]<- sum(GC$inboundTumorcells==1 & GC$SOX9_prolif == lineages[i])/ sum(GC$inboundTumorcells==1 & GC$SOX9_prolif != "")
  Lineage_porp_int[i]<-  sum(GC$inboundTumorcells==0 & GC$Tumor==1 & GC$SOX9_prolif == lineages[i])/ sum(GC$inboundTumorcells==0 & GC$Tumor ==1 & GC$SOX9_prolif != "")  
  Lineage_porp_inbs[i]<- sum(GC$inboundStromcells==1 & GC$SOX9_prolif == lineages[i])/ sum(GC$inboundStromcells==1 & GC$SOX9_prolif != "")   
  Lineage_porp_ins[i]<-  sum(GC$inboundStromcells==0 & GC$Tumor==0 & GC$SOX9_prolif == lineages[i])/ sum(GC$inboundStromcells==0 & GC$Tumor ==0 & GC$SOX9_prolif != "")    
}
bardata<-rbind(Lineage_porp_int, Lineage_porp_inbt, -Lineage_porp_inbs, -Lineage_porp_ins)
rownames(bardata)<- c("Tumor", "Tumor border", "Stroma border", "Stroma")
colnames(bardata)<- c("SOX9+ nP", "SOX9- nP", "SOX9+ P", "SOX9- P")
barplot(height = bardata, beside = TRUE, col=c("darkred", "red", "orange", "yellow"), width = 2, ylim = c(-1,1), las=2)
title(main = paste("SOX9 Proliferation"), font.main = 4, sub="")

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
legend("center", legend = rownames(bardata), fill = c("darkred", "red", "orange", "yellow"), title = "Legend", cex=1.2)

