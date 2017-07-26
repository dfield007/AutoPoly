########################################
# Power testing Cluster (BaseContinue) #
########################################

# Linux cluster setup
setwd("/data/david.field/")
library(gtools)                            

############################
#  Forward simulation code #
############################

source("getMotherGameteInfo.R")
source("importDataFile.r")            
source("mutation.r")             
source("forwardSimBaseStart.r")
source("forwardSimBaseContinue.r")
source("forwardSimMetaPopStart.r")            
source("forwardSimMetaPopContinue.r")

# loading previous workspace..... t=1
load("simAncestral_x4_t1_Workspace.RData")

# Parameters (continuing from initial startup - forwardSimBase)
lociSim<-24
ploidy<-4
DRR<-"min"
popSizeBase<-10000
n1<-500
mutationRate<-0.001
mutationModel<-"SMM"

# Saving data at..... t=2
simAncestral<-fowardSimBaseContinue(simAncestral,lociSim,ploidy,DRR,popSizeBase,n1,mutationRate,mutationModel)
write.csv(simAncestral[[7]],"simAncestral_x4_t2_genSoFar.csv",row.names=FALSE)
write.csv(simAncestral[[3]],"simAncestral_x4_t2_MeanAlleles.csv",row.names=FALSE)
write.csv(simAncestral[[1]],"simAncestral_x4_t2.csv",row.names=FALSE)
save.image("simAncestral_x4_t2_Workspace.RData") 
