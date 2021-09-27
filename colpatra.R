# Phylogenetic analyses for Trujillo et al. in prep "Palm functional trait responses to local environmental factors in the Colombian Amazon", submitted to JTE
# Wolf Eiserhardt, 27/9/2021
# wolf.eiserhardt@bio.au.dk
# https://github.com/pebgroup

# load packages
library("phangorn")
library("picante")
library("PCPS")
library("HDInterval")

# load modified functions (derived from phangorn and PCPS)
source("pcps.curve_WE.R")
source("densiTree_WE.R")

# load community data
load("Data.RData")

# load trees
name_matching <- read.table("names.csv", sep=";", header=T)
rownames(name_matching) <- name_matching$tree_name

trees1 <- read.nexus("Phylogeny_Con_Checklist_1.nex")
trees1_pruned <- lapply(trees1, keep.tip, tip=as.vector(name_matching$tree_name))

trees2 <- read.nexus("Phylogeny_Con_Checklist_2.nex")
trees2_pruned <- lapply(trees2, keep.tip, tip=as.vector(name_matching$tree_name))

trees <- c(trees1_pruned, trees2_pruned)
class(trees) <- "multiPhylo"

mcctree <- maxCladeCred(trees)
tree_sample <- trees[sample(1:1000, 150)]

#######################################
# tree figure (Supplement 4 Figure 1) #
#######################################

pdf("phylotraits.pdf", width=8.3, height = 11.7)
densiTree_WE(tree_sample, col="blue", width=3, tip.color="black", cex=1, label.offset = 0.45, plotwidth = 3.3)

order <- c("maurarma", "iriaseti", "socrexor", "manisacc", "geondeve", "geonmaxi", "geonmacr", "hyoseleg", "euteprec", "oenobata", "oenobaca", "attamari", "attamicr", "attabuty", "desmgiga", "desmpoly", "desmmiti", "astrgyna", "bactmara", "bacthirt", "bactcoro", "bactacan", "bactbron", "bactmajo", "bactsimp")
scf <- 2

# LF
symbol <- rep(19, 25)
char <- as.character(traits$LF)
names(char) <- traits$spp
char <- char[rev(order)]
symbol[char=="sol"] <- 1
points(x=rep(1.05,25), y=1:25, pch=symbol, cex=.6)

text(1.05, 26, "LF", srt=90, adj=0)

# GF
symbol <- rep(19, 25)
char <- as.character(traits$GF)
names(char) <- traits$spp
char <- char[rev(order)]
symbol[char=="aca"] <- 1
symbol[char=="cli"] <- 0
points(x=rep(1.1,25), y=1:25, pch=symbol, cex=.6)

text(1.1, 26, "GF", srt=90, adj=0)

# StH
char <- traits$StH
names(char) <- traits$spp
char <- char[rev(order)]
points(x=rep(1.15,25), y=1:25, cex=scf*char/max(char, na.rm=T))

text(1.15, 26, "StH", srt=90, adj=0)

# LN
char <- traits$LN
names(char) <- traits$spp
char <- char[rev(order)]
points(x=rep(1.2,25), y=1:25, cex=scf*char/max(char, na.rm=T))

text(1.2, 26, "LN", srt=90, adj=0)

# PeL
char <- traits$PeL
names(char) <- traits$spp
char <- char[rev(order)]
points(x=rep(1.25,25), y=1:25, cex=scf*char/max(char, na.rm=T))

text(1.25, 26, "PeL", srt=90, adj=0)

# RL
char <- traits$RL
names(char) <- traits$spp
char <- char[rev(order)]
points(x=rep(1.3,25), y=1:25, cex=scf*char/max(char, na.rm=T))

text(1.3, 26, "RL", srt=90, adj=0)

# FD
char <- traits$FD
names(char) <- traits$spp
char <- char[rev(order)]
points(x=rep(1.35,25), y=1:25, cex=scf*char/max(char, na.rm=T))

text(1.35, 26, "FD", srt=90, adj=0)

# SN
char <- traits$SN
names(char) <- traits$spp
char <- char[rev(order)]
points(x=rep(1.4,25), y=1:25, cex=scf*char/max(char, na.rm=T))

text(1.4, 26, "SN", srt=90, adj=0)

dev.off()


#########################################
# PCPS analysis (Supplement 4 Figure 2) #
#########################################

# format species data correctly
spe2 <- spe
rownames(spe2) <- spe2$spp
spe2 <- spe2[2:26]

pcps(spe2, cophenetic(mcctree)) -> P

pdf("comphylosignal.pdf", width=6.3, height = 11.7)
par(mfrow=c(4,2))

#LF
#ces
LF4pcps <- traits$LF
LF4pcps[LF4pcps=="sol"] <- 0
LF4pcps[LF4pcps=="ces"] <- 1
LF4pcps <- as.numeric(LF4pcps)
dim(LF4pcps) <- c(length(LF4pcps), 1)

rownames(LF4pcps) <- traits$spp
colnames(LF4pcps) <- "LF"

pcps.curve_WE(spe2, cophenetic(mcctree), LF4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_LF
plot(Pcurve_LF, draw.model="ts", model.col="grey")
title("LF:cespitose")

#GF
#cli
GF4pcps <- traits$GF
GF4pcps[GF4pcps=="aca"] <- 0
GF4pcps[GF4pcps=="ere"] <- 0
GF4pcps[GF4pcps=="cli"] <- 1
GF4pcps <- as.numeric(GF4pcps)
dim(GF4pcps) <- c(length(GF4pcps), 1)

rownames(GF4pcps) <- traits$spp
colnames(GF4pcps) <- "GF"

pcps.curve_WE(spe2, cophenetic(mcctree), GF4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_GF
plot(Pcurve_GF, draw.model="ts", model.col="grey")
title("GF:climbing")

#GF
#aca
GF4pcps <- traits$GF
GF4pcps[GF4pcps=="aca"] <- 1
GF4pcps[GF4pcps=="ere"] <- 0
GF4pcps[GF4pcps=="cli"] <- 0
GF4pcps <- as.numeric(GF4pcps)
dim(GF4pcps) <- c(length(GF4pcps), 1)

rownames(GF4pcps) <- traits$spp
colnames(GF4pcps) <- "GF"

pcps.curve_WE(spe2, cophenetic(mcctree), GF4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_GF
plot(Pcurve_GF, draw.model="ts", model.col="grey")
title("GF:acaulescent")

#StH
StH4pcps <- traits$StH
dim(StH4pcps) <- c(length(StH4pcps), 1)
rownames(StH4pcps) <- traits$spp
colnames(StH4pcps) <- "StH"

pcps.curve_WE(spe2, cophenetic(mcctree), StH4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_StH
plot(Pcurve_StH, draw.model="ts", model.col="grey")
title("StH")

#LN
LN4pcps <- traits$LN
dim(LN4pcps) <- c(length(LN4pcps), 1)
rownames(LN4pcps) <- traits$spp
colnames(LN4pcps) <- "LN"

pcps.curve_WE(spe2, cophenetic(mcctree), LN4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_LN
plot(Pcurve_LN, draw.model="ts", model.col="grey")
title("LN")

#PeL
PeL4pcps <- traits$PeL
dim(PeL4pcps) <- c(length(PeL4pcps), 1)
rownames(PeL4pcps) <- traits$spp
colnames(PeL4pcps) <- "PeL"

pcps.curve_WE(spe2, cophenetic(mcctree), PeL4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_PeL
plot(Pcurve_PeL, draw.model="ts", model.col="grey")
title("PeL")

#RL
RL4pcps <- traits$RL
dim(RL4pcps) <- c(length(RL4pcps), 1)
rownames(RL4pcps) <- traits$spp
colnames(RL4pcps) <- "RL"

pcps.curve_WE(spe2, cophenetic(mcctree), RL4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_RL
plot(Pcurve_RL, draw.model="ts", model.col="grey")
title("RL")

#FD
FD4pcps <- traits$FD
dim(FD4pcps) <- c(length(FD4pcps), 1)
rownames(FD4pcps) <- traits$spp
colnames(FD4pcps) <- "FD"

pcps.curve_WE(spe2, cophenetic(mcctree), FD4pcps, null.model.ts = TRUE, runs=999) -> Pcurve_FD
plot(Pcurve_FD, draw.model="ts", model.col="grey")
title("FD")

dev.off()

#####################################################
# trait phylogenetic signal (Supplement 4, Table 1) #
#####################################################

# Sorry, this part is a bit messy (assembling table by hand rather than generating it here)

trees_renamed <- trees
for(i in 1:1000){
  trees_renamed[[i]]$tip.label <- name_matching[trees_renamed[[i]]$tip.label,"spp"]
}

# categorical traits

# life form: 
samp <- matrix(nrow=nrow(traits), ncol=2)
samp[,] <- 0
rownames(samp) <- traits$spp
colnames(samp) <- c("sol", "ces")
samp[traits$LF=="ces","ces"] <- 1
samp[traits$LF=="sol","sol"] <- 1
samp <- t(samp)

mpd.obs.z <- rep(NA, 1000)
p <- rep(NA, 1000)
for(j in 1:1000){
  ps <- ses.mpd(samp, cophenetic(trees_renamed[[i]]), "taxa.labels")
  mpd.obs.z[j] <- ps["ces","mpd.obs.z"]
  p[j] <- ps["ces","mpd.obs.p"]
}  

round(mean(mpd.obs.z),2)
round(hdi(mpd.obs.z)[1],2)
round(hdi(mpd.obs.z)[2],2)
round(mean(p),3)
round(hdi(p)[1],3)
round(hdi(p)[2],3)

# growth form: 
samp <- matrix(nrow=nrow(traits), ncol=3)
samp[,] <- 0
rownames(samp) <- traits$spp
colnames(samp) <- c("ere", "aca", "cli")
samp[traits$GF=="ere","ere"] <- 1
samp[traits$GF=="aca","aca"] <- 1
samp[traits$GF=="cli","cli"] <- 1
samp <- t(samp)

mpd.obs.z_aca <- rep(NA, 1000)
p_aca <- rep(NA, 1000)
mpd.obs.z_cli <- rep(NA, 1000)
p_cli <- rep(NA, 1000)
for(j in 1:1000){
  ps <- ses.mpd(samp, cophenetic(trees_renamed[[i]]), "taxa.labels")
  mpd.obs.z_aca[j] <- ps["aca","mpd.obs.z"]
  p_aca[j] <- ps["aca","mpd.obs.p"]
  mpd.obs.z_cli[j] <- ps["cli","mpd.obs.z"]
  p_cli[j] <- ps["cli","mpd.obs.p"]
}  

round(mean(mpd.obs.z_aca),2)
round(hdi(mpd.obs.z_aca)[1],2)
round(hdi(mpd.obs.z_aca)[2],2)
round(mean(p_aca),3)
round(hdi(p_aca)[1],3)
round(hdi(p_aca)[2],3)

round(mean(mpd.obs.z_cli),2)
round(hdi(mpd.obs.z_cli)[1],2)
round(hdi(mpd.obs.z_cli)[2],2)
round(mean(p_cli),3)
round(hdi(p_cli)[1],3)
round(hdi(p_cli)[2],3)

# continuous traits
phylosignal_results <- list()

#for(trait in colnames(traits)[4:4]){ #9
for(i in 1:6){ #9
    
  char <- traits[,i+3]
  names(char) <- traits$spp
  
  k <- rep(NA, 1000)
  p <- rep(NA, 1000)
  for(j in 1:1000){
    ps <- phylosignal(char, trees_renamed[[j]])
    k[j] <- ps$K
    p[j] <- ps$PIC.variance.P
  }
  phylosignal_results[[i]] <- cbind(k, p)
}

# print means + HPDs for table 1
#par(mfrow=c(3,2))
for(i in 1:6){
  #K
  #print(paste(round(mean(phylosignal_results[[i]][,1]),2), "[", round(hdi(phylosignal_results[[i]][,1])[1],2), round(hdi(phylosignal_results[[i]][,1])[2],2),  "]"))
  #P
  #print(paste(round(mean(phylosignal_results[[i]][,2]),3), "[", paste(round(quantile(phylosignal_results[[i]][,2], c(0.025,0.975)),3), collapse=" - "), "]"))
  print(paste(round(mean(phylosignal_results[[i]][,2]),3), "[", round(hdi(phylosignal_results[[i]][,2])[1],3), round(hdi(phylosignal_results[[i]][,2])[2],3),  "]"))
  #hist(phylosignal_results[[i]][,1], xlim=c(0,2.5))
}


           