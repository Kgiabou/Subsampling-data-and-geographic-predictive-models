### SUBSAMPLING ANALYSIS ####

setwd("C:/Users/kgiab/Desktop/Michela review analysis")
load("Horse_climatic_data2.RData")
library(ecospat)
library(biomod2)
library(ade4)
library(adehabitatHS)
library(MASS)

bins <- c(seq(3,22,1), seq(24, 44,2))
occs_Asia <- read.table("DBAsia_variables.csv", h=T, sep=";", stringsAsFactors = FALSE)
occs_Asia$Lab_Code <- rownames(occs_Asia)
rownames(occs_Asia) <- c(1:nrow(occs_Asia))
occs_Asia <- occs_Asia[, c(1,2,4)]
colnames(occs_Asia) <- c("x","y","Interval")

occs_Europe <- read.table("DBEurope_variables.csv", h=T, sep=";", stringsAsFactors = FALSE)
occs_Europe$Lab_Code <- rownames(occs_Europe)
rownames(occs_Europe) <- c(1:nrow(occs_Europe))
occs_Europe <- occs_Europe[, c(1,2,4)]
colnames(occs_Europe) <- c("x","y","Interval")

## Set iteration number ##
repl <- 100
Overlap <- matrix(NA, nrow=repl, ncol=1)
dimnames(Overlap) <- list(paste("round_", c(1:repl), sep=""), "Overlap")

## For each iteration subset the data and perform the niche overlap analysis## 
for (samples in 1:repl)
{
for(r in seq_along(bins))
{
## subset cliamtic data for each time bin ##
cl_As <- subset(clim_all_Asia, Interval==bins[r])
cl_Eur <- subset(clim_all_Europe, Interval==bins[r])
## Subset Asian occurrences per time bin ##
df_A <- subset(occs_Asia, Interval==bins[r])
sams <- nrow(df_A)
## Subset European occurrences per time bin ##
df_E <- subset(occs_Europe, Interval==bins[r])
if(sams>0)
{
## Check if European data are more than Asian data ##
if(nrow(df_E)>sams)
{
## Sample European occurrences at the same level as Asian occurrences ##
vec_sam <- as.numeric(sample(1:nrow(df_E), sams))
df_E_sam <- df_E[vec_sam,]
}
else if(nrow(df_E)<sams)
## IF Asian data are more samplle in the opposite direction##
{
df_E_sam <- df_E
vec_As <- as.numeric(sample(1:nrow(df_A), nrow(df_E)))
df_A <- df_A[vec_As,]
}
if(nrow(df_A)>0)
{
occ_A <- na.exclude(ecospat.sample.envar(dfsp=df_A, colspxy=2:1,colspkept=NULL,dfvar=cl_As,colvarxy=1:2,colvar="all",resolution=0.1666565))
assign(paste("occ_Asia_", bins[r], "k", sep=""), occ_A)
}
if(nrow(df_E)>0)
{
occ_E <- na.exclude(ecospat.sample.envar(dfsp=df_E_sam, colspxy=2:1,colspkept=NULL,dfvar=cl_Eur, colvarxy=1:2,colvar="all",resolution=0.1666565))
assign(paste("occ_Europe_", bins[r], "k", sep=""), occ_E)
}
}
}

As <- ls(pattern = "occ_Asia")[c(22, 26:30, 1:21, 23:25)]
As_all <- lapply(As, get)
occ_Asia_all <- do.call(rbind, As_all)

Es <- ls(pattern = "occ_Europe_")[c(22, 26:31, 1:21, 23:25)]
Es_all <- lapply(Es, get)
occ_Europe_all <- do.call(rbind, Es_all)

n_bins_As <- as.vector(unique(occ_Asia_all$Interval))
occ_Europe_all2 <- occ_Europe_all

clim_all_Europe <- clim_all_Europe[clim_all_Europe$Interval %in% n_bins_As,]
clim_all_Asia <- clim_all_Asia[clim_all_Asia$Interval %in% n_bins_As,]
clim_all <- rbind(clim_all_Europe, clim_all_Asia)

 ## Perfomr Ecospat SDM overlap ##
PROJ = F
Xvar<-c(3:6)
nvar<-length(Xvar)
iterations<-100
R=100

# if PROJ = F
row.w.Eur.occ<-1-(nrow(occ_Europe_all2)/nrow(rbind(occ_Europe_all2, occ_Asia_all))) # prevalence of occ1
row.w.Asia.occ<-1-(nrow(occ_Asia_all)/nrow(rbind(occ_Europe_all2, occ_Asia_all))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim_all_Europe)),rep(0, nrow(clim_all_Asia)),rep(row.w.Eur.occ, nrow(occ_Europe_all2)),rep(row.w.Asia.occ, nrow(occ_Asia_all)))

row.w.Eur.env<-1-(nrow(clim_all_Europe)/nrow(clim_all))  # prevalence of clim1
row.w.Asia.env<-1-(nrow(clim_all_Asia)/nrow(clim_all))  # prevalence of clim2
row.w.env<-c(rep(row.w.Eur.env, nrow(clim_all_Europe)),rep(row.w.Asia.env, nrow(clim_all_Asia)),rep(0, nrow(occ_Europe_all2)),rep(0, nrow(occ_Asia_all)))

fac<-as.factor(c(rep(1, nrow(clim_all_Europe)),rep(2, nrow(clim_all_Asia)),rep(1, nrow(occ_Europe_all2)),rep(2, nrow(occ_Asia_all))))

# global dataset for the analysis and rows for each sub dataset
data.env.occ<-rbind(clim_all_Europe,clim_all_Asia,occ_Europe_all2,occ_Asia_all)[Xvar]
row.clim_Eur<-1:nrow(clim_all_Europe)
row.clim_As<-(nrow(clim_all_Europe)+1):(nrow(clim_all_Europe)+nrow(clim_all_Asia))
row.clim_all<-1:(nrow(clim_all_Europe)+nrow(clim_all_Asia))
row.sp_Eur<-(nrow(clim_all_Europe)+nrow(clim_all_Asia)+1):(nrow(clim_all_Europe)+nrow(clim_all_Asia)+nrow(occ_Europe_all2))
row.sp_As<-(nrow(clim_all_Europe)+nrow(clim_all_Asia)+nrow(occ_Europe_all2)+1):(nrow(clim_all_Europe)+nrow(clim_all_Asia)+nrow(occ_Europe_all2)+nrow(occ_Asia_all))

rm(list=c("clim_all"))

if(PROJ == F){	#fit of the analyse using occurences from both ranges		
	pca.cal <-dudi.pca(data.env.occ, row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
}

rm(list=c("data.env.occ"))
			# predict the scores on the axes
scores.clim_all<- pca.cal$li[row.clim_all,]
scores.clim_Eur<- pca.cal$li[row.clim_Eur,]
scores.clim_As<- pca.cal$li[row.clim_As,]
scores.sp_Eur<- pca.cal$li[row.sp_Eur,]
scores.sp_As<- pca.cal$li[row.sp_As,]


			# calculation of occurence density and test of niche equivalency and similarity 
z_Eur<- ecospat.grid.clim.dyn(scores.clim_all,scores.clim_Eur,scores.sp_Eur,R)
z_As<- ecospat.grid.clim.dyn(scores.clim_all,scores.clim_As,scores.sp_As,R)

# Calculate Overlap	for all iterations		
Over <-round(as.numeric(ecospat.niche.overlap(z_Eur,z_As,cor=T)[1]),3)
Overlap[samples,1] <- Over
rm(list=c("z_As", "z_Eur","pca.cal","scores.clim_all", "scores.clim_Eur","scores.clim_As","scores.sp_Eur", "scores.sp_As"))
write.table(Overlap, file="Horses_overlap_subsampling.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
}

## Plot Histogram of Overlap values ##
h <-hist.default(bb$Overlap, probability = TRUE, col="steelblue",
                 border = "grey", axes=F, xaxs="i", yaxs="i", main="Resampling Overlap Values")
				
abline(v=0.123, col="orange", lwd=2, lty=2)
axis(2, at=seq(0,65,5))
axis(1, at=h$mids, labels=round(h$mids, 3))
mtext(side=1, line = 2.5, cex=1.15, "Overlap")
mtext(side=2, line = 2.5, cex=1.15, "Counts")

## Plot niche changes in climatic space ##
trellis.device(device="pdf", file="Niche_Overlap_all.pdf", width=14, height=10)
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE))
ecospat.plot.niche(z_Eur,title="PCA-env - EU niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.niche(z_As,title="PCA-env - AS niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(ecospat.niche.overlap(z_Eur,z_As,cor=T)[1]),3)))
ecospat.plot.overlap.test(a,"D","Equivalency")
ecospat.plot.overlap.test(b,"D","Similarity 2->1")
ecospat.plot.overlap.test(b2,"D","Similarity 1->2")
dev.off()
