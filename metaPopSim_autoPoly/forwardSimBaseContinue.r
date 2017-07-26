###################################
#  Base population (Ancestral)    #
# Forward simulations, continue...#
###################################

fowardSimBaseContinue<- function(simPop,lociSim,ploidy,DRR,popSizeBase,n1,mutationRate,mutationModel) {

cat("\n ********************************************** \n")
        cat(" *          PolyAssign (1/1/2010)             * \n")
        cat(" *   Authors: David Field, Alec Zwart,        * \n")
        cat(" *            Andrew Young, Linda Broadhurst  * \n")
        cat(" ********************************************** \n")
        cat("\n ** Beginning Forward Simulation.... **\n")

  flush.console()
    #############################################
    # a) Directing required hexaploid functions #
    #############################################
  if (ploidy==6){
      #getProbGenotypeGivenPopn <- genotypeFreqPopHex           # freq of a genotype, given gamete freqs of population
      #mixedGenotypeFreqGivenPopns <- mixedGenotypeFreqPopHex   # freq of a mixed F1 genotype, given gamete freqs of 2 populations
      #getSeedlingGenotypes <- getSeedlingGenotypesHex          # list of possible genotypes for a phenotype
      #gameteFreqGivenPopn <- gameteFrequencyPopHex             # frequency of a given gamete
      if (DRR == "min") {
          getMotherGameteInfo <- getMotherGameteInfoHexGenoMin # segregation ratios: gametes from individual
      } else if (DRR == "max") {
          getMotherGameteInfo <- getMotherGameteInfoHexGenoMax # segregation ratios: gametes from individual
      }
  }
    ################################################
    # b) Directing required tetraploid functions   #
    ################################################

  if (ploidy==4) {
      #getProbGenotypeGivenPopn <- genotypeFreqPopTetr              # freq of a genotype, given gamete freqs of population
      #mixedGenotypeFreqGivenPopns <- mixedGenotypeFreqPopTetr      # freq of a mixed F1 genotype, given gamete freqs of 2 populations
      #getSeedlingGenotypes <- getSeedlingGenotypesTetr             # list of possible genotypes for a phenotype
      #gameteFreqGivenPopn <- gameteFrequencyPopTetr                # frequency of a given gamete
      if (DRR == "min") {
          getMotherGameteInfo <- getMotherGameteInfoTetrGenoMin # segregation ratios: gametes from individual
      } else if (DRR == "max") {
          getMotherGameteInfo <- getMotherGameteInfoTetrGenoMax # segregation ratios: gametes from individual
      }
  }
  strip.NAs <- function(vv) {
      return(vv[!is.na(vv)])
  }
  # Making loci column header
  if (ploidy==4) {
      locusValues<-c("a","b","c","d")

  }
  if (ploidy==6) {
      locusValues<-c("a","b","c","d","e","f")
  }
  lociLabels<-NULL
  for (thisLocus in 1:lociSim) {
      # thisLocus<-3
      lociLabels<-c(lociLabels,paste("Locus",thisLocus,locusValues,sep=""))
  }
  
  ##############################
  # Pull out previous run data #
  ############################## 
  # Genotypes
  matrixT1<- simPop[[1]]
  # Locus bounds information
  locusBounds<-simPop[[2]]
  # Generation Tally
  GenerationsSoFar<- simPop[[7]]

  # Empty matrix ready for initial allele data
  matrixT2<-matrixT1

  ##############################
  # Starting at generation t=1 #
  ##############################
  for (thisGen in 1:n1) {
       total<-n1
       # Clear out previous T2
       matrixT2[,4:(length(lociLabels)+3)]<-0
      #  thisGen<-1
          if (thisGen==as.integer((total)*0.1)) {
              cat("\n         10%..") }
          if (thisGen==as.integer((total)*0.2)) {
              cat("20%..") }
          if (thisGen==as.integer((total)*0.3)) {
              cat("30%..") }
          if (thisGen==as.integer((total)*0.4)) {
              cat("40%..") }
          if (thisGen==as.integer((total)*0.5)) {
              cat("50%..") }
          if (thisGen==as.integer((total)*0.6)) {
              cat("60%..") }
          if (thisGen==as.integer((total)*0.7)) {
              cat("70%..") }
          if (thisGen==as.integer((total)*0.8)) {
              cat("80%..") }
          if (thisGen==as.integer((total)*0.9)) {
              cat("90%..") }
          if (thisGen==as.integer((total)*1)) {
              cat("Complete \n") }
          flush.console()

      #################
      # Random Mating #
      #################
      # Randomly drawn parent for progeny in t + 1
      parents1<-sample(1:popSizeBase,popSizeBase,replace=TRUE)
      parents2<-sample(1:popSizeBase,popSizeBase,replace=TRUE)
      newProgeny<-matrix(0,popSizeBase,length(lociLabels))
      
      ####################################
      #  Begin loop over new individuals #
      ####################################
      for (thisProgeny in 1:popSizeBase) {
          # thisProgeny<-1
          ##################
          # Loop over loci #
          ##################
          for (thisLocus in 1:lociSim) {
              # thisLocus<-1
              # Current locus column range
              locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
              locusRangeNewProgeny<-(thisLocus-1)*ploidy + 1:ploidy
              # Parent1 gamete
              Parent1Genotype <- matrixT1[parents1[thisProgeny],locusRange]
              names(Parent1Genotype)<-NULL
              Parent1GametesList <-getMotherGameteInfo(Parent1Genotype)
              Parent1Gametes<- names(Parent1GametesList$prob)
              Parent1Probs<- Parent1GametesList$prob
              names(Parent1Probs)<-NULL
              thisParent1.sample<- sample(Parent1Gametes,1,replace=TRUE,Parent1Probs)
              thisParent1.sample<- t(matrix(unlist(strsplit(thisParent1.sample," ")),ploidy/2,))
              # Parent2 gamete
              Parent2Genotype <- matrixT1[parents2[thisProgeny],locusRange]
              names(Parent2Genotype)<-NULL
              Parent2GametesList <-getMotherGameteInfo(Parent2Genotype)
              Parent2Gametes<- names(Parent2GametesList$prob)
              Parent2Probs<- Parent2GametesList$prob
              names(Parent2Probs)<-NULL
              thisParent2.sample<- sample(Parent2Gametes,1,replace=TRUE,Parent2Probs)
              thisParent2.sample<- t(matrix(unlist(strsplit(thisParent2.sample," ")),ploidy/2,))
              # New individual for t+1
              newProgenyOne<-cbind(thisParent1.sample,thisParent2.sample)
              # Mutation
              newProgeny[thisProgeny,locusRangeNewProgeny]<-mutation(newProgenyOne,thisLocus,mutationRate,mutationModel,locusBounds)
          }  # End locus loop
      } # End Progeny loop
      
      ######################
      # Send to Matrix T2  #
      ######################
      matrixT2[,4:(length(lociLabels)+3)]<-newProgeny
      matrixT1<-matrixT2

    } #End Generation loop

  ##############################
  # extracting no. of alleles  #
  ##############################
  popGenStatsMatrix<-matrix(0,1,lociSim)
  rownames(popGenStatsMatrix)<-1
  for (thisLocus in 1:lociSim) {
      # thisLocus<-2
      locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
      thisLocusData<-matrixT1[,locusRange]
      popnAlleles <- sort(strip.NAs(unique(unlist(as.data.frame(thisLocusData)))))
      numAlleles<-length(popnAlleles)
      popGenStatsMatrix[1,thisLocus]<-numAlleles
      }

      meanAlleles<-apply(popGenStatsMatrix,1,mean)
      SDAlleles<-apply(popGenStatsMatrix,1,sd)
      maxNumAlleles<-apply(popGenStatsMatrix,1,max)
      minNumAlleles<-apply(popGenStatsMatrix,1,min)
      # Send data to list
      resultsForwardSim<-list()
      resultsForwardSim[["SimData"]]<-list()
      resultsForwardSim[["LocusBounds"]]<-list()
      resultsForwardSim[["Mean No. Alleles"]]<-list()
      resultsForwardSim[["Standard deviation alleles number"]]<-list()
      resultsForwardSim[["Max No. Alleles"]]<-list()
      resultsForwardSim[["Min No. Alleles"]]<-list()
      resultsForwardSim[["Generations"]]<-list()
      resultsForwardSim[["SimData"]]<-matrixT1
      resultsForwardSim[["LocusBounds"]]<-locusBounds
      resultsForwardSim[["Mean No. Alleles"]]<-meanAlleles
      resultsForwardSim[["Standard deviation alleles number"]]<-SDAlleles
      resultsForwardSim[["Max No. Alleles"]]<-maxNumAlleles
      resultsForwardSim[["Min No. Alleles"]]<-minNumAlleles
      resultsForwardSim[["Generations"]]<-n1+GenerationsSoFar
      return(resultsForwardSim)                               
      } # close function

