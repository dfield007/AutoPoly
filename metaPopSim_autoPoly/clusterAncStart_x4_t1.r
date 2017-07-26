#####################################
# Power testing Cluster (BaseStart) #
#####################################

# MAC setup
setwd("/Users/dfield/Documents/DavidWorkLaptop/PostdocCSIRO/final2017/package/SimulationIslandModel/")
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

# Parameters
lociSim<-24
ploidy<-4
DRR<-"min"
popSizeBase<-10000
n1<-500
mutationRate<-0.001
mutationModel<-"SMM"

# Saving data at t=1
simAncestral<- fowardSimBase(lociSim,ploidy,DRR,popSizeBase,n1,mutationRate,mutationModel)
write.csv(simAncestral[[2]],"simAncestral_x4_t1_locusBounds.csv",row.names=FALSE)
write.csv(simAncestral[[3]],"simAncestral_x4_t1_MeanAlleles.csv",row.names=FALSE)
write.csv(simAncestral[[7]],"simAncestral_x4_t1_genSoFar.csv",row.names=FALSE)
write.csv(simAncestral[[1]],"simAncestral_x4_t1.csv",row.names=FALSE)
save.image("simAncestral_x4_t1_Workspace.RData")
