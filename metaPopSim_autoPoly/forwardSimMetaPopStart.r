###############################################
#  Forward simulations Metapopulation (Start) #
###############################################

fowardSimMeta<- function(lociSim,ploidy,DRR,popSizeSims,numPopsSims,n1,
                      n2,mutationRate,mutationModel,minVariability,maxAlleles,migration) {

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

  # Empty matrix ready for initial allele data
  matrixT1<-matrix(0,popSizeSims*numPopsSims,3+length(lociLabels))
  colnames(matrixT1)<-c("ID","pop","mother",lociLabels)
  rownames(matrixT1)<-1:(popSizeSims*numPopsSims)
  matrixT1[,"ID"]<-1:(popSizeSims*numPopsSims)
  popLabels<-NULL
  for (thisPop in 1:numPopsSims) {
      popLabels<-c(popLabels,rep(thisPop,popSizeSims))
  }
  matrixT1[,"pop"]<-popLabels
  matrixT1[,"mother"]<-NA
  #matrixT1<-as.data.frame(matrixT1)
  matrixT2<-matrixT1
  # head(matrixT1)
  # unique(matrixT1[,"pop"])

  ########################
  # Starting variability #
  ########################
  # Minimum genetic variability
  if (minVariability==TRUE) {
      beginAlleles<-as.integer(runif(lociSim,150,380))
      for (thisLocus in 1:lociSim) {
          # thisLocus<-1
          # Current locus column range
          locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
          matrixT1[,locusRange]<-beginAlleles[thisLocus]
      }

      # determining boundaries for allele sizes per locus
      locusBounds<-matrix(0,lociSim,2)
      colnames(locusBounds)<-c("min","max")
      for (thisLocus in 1:lociSim) {
          locusBounds[thisLocus,"min"]<-beginAlleles[thisLocus]-30
          locusBounds[thisLocus,"max"]<-beginAlleles[thisLocus]+30
      }
  allelicStates <- matrix(0,lociSim,maxAlleles)
  }
  # Maximum genetic variability
  if (minVariability==FALSE) {
      beginAlleles<-as.integer(runif(lociSim,150,380))
      # determining boundaries for allele sizes per locus
      locusBounds<-matrix(0,lociSim,2)
      colnames(locusBounds)<-c("min","max")
      for (thisLocus in 1:lociSim) {
          locusBounds[thisLocus,"min"]<-beginAlleles[thisLocus]-30
          locusBounds[thisLocus,"max"]<-beginAlleles[thisLocus]+30
      }
      # Determining the allelic states
      allelicStates <- matrix(0,lociSim,maxAlleles)
      for (thisLocus in 1:lociSim) {
           # thisLocus<-1
           allelesPossible <- as.vector(locusBounds[thisLocus,1]):as.vector(locusBounds[thisLocus,2])
           allelesThisLocus <- sample(allelesPossible,maxAlleles,replace=FALSE)
           allelicStates[thisLocus,] <- allelesThisLocus
      }
      # Fill matrixT1 with random alleles from possible allelic states
      for (thisInd in 1:nrow(matrixT1)) {
          # thisInd <- 1
          for (thisLocus in 1:lociSim) {
              # thisLocus<-1
              # Current locus column range
              locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
              thisIndGenotype <- sample(allelicStates[thisLocus,],ploidy)
              matrixT1[thisInd,locusRange]<-thisIndGenotype
          }
      }
  }

  ##############################
  # Starting at generation t=1 #
  ##############################
  for (thisGen in 1:(n1+n2)) {
       total<-n1+n2
      # thisGen<-2
      # for (thisGen in 1:total) {
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
       #}

      #################
      # Random Mating #
      #################
      # Empty matrix for next generation t + 1
      #matrixT2<-matrix(0,popSizeSims*numPopsSims,3+length(lociLabels))
      #colnames(matrixT2)<-c("ID","pop","mother",lociLabels)
      #rownames(matrixT2)<-1:(popSizeSims*numPopsSims)
      #matrixT2[,"ID"]<-1:(popSizeSims*numPopsSims)
      #popLabels<-NULL
      #for (thisPop in 1:numPopsSims) {
      #    popLabels<-c(popLabels,rep(thisPop,popSizeSims))
      #}
      #matrixT2[,"pop"]<-popLabels
      #matrixT2[,"mother"]<-NA
      #matrixT2<-as.data.frame(matrixT2)
      # Tally matrix for individuals filled in T2
      statusT2 <- matrix(0,2,numPopsSims)
      rownames(statusT2) <- c("Full","Tally")
      colnames(statusT2) <- 1:numPopsSims

      ####################################
      #  Begin loop over new individuals #
      ####################################
      while (sum(statusT2["Full",]) < numPopsSims) {

          # Drawing parents from matrix t for random mating
          parent1 <- sample(1:(popSizeSims*numPopsSims),1,replace=TRUE)
          thisRandomPop <- matrixT1[parent1,"pop"]
          thisRandomPop <- as.numeric(thisRandomPop)
          thisRandomPopInds <- as.logical(matrixT1[,"pop"] == thisRandomPop)
          thisRandomPopInds <- as.vector(matrixT1[thisRandomPopInds,"ID"])
          parent2 <- sample(thisRandomPopInds,1,replace=TRUE)
          newProgeny <- matrix(0,1,length(lociLabels))

          ##################
          # Loop over loci #
          ##################
          for (thisLocus in 1:lociSim) {
              # thisLocus<-1
              # Current locus column range
              locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
              locusRangeNewProgeny <- (thisLocus-1)*ploidy + 1:ploidy
              # Parent1 gamete
              Parent1Genotype <- matrixT1[parent1,locusRange]
              names(Parent1Genotype)<-NULL
              Parent1GametesList <- getMotherGameteInfo(Parent1Genotype)
              Parent1Gametes<- names(Parent1GametesList$prob)
              Parent1Probs<- Parent1GametesList$prob
              names(Parent1Probs)<-NULL
              thisParent1.sample<- sample(Parent1Gametes,1,replace=TRUE,Parent1Probs)
              thisParent1.sample<- t(matrix(unlist(strsplit(thisParent1.sample," ")),ploidy/2,))
              # Parent2 gamete
              Parent2Genotype <- matrixT1[parent2,locusRange]
              names(Parent2Genotype)<-NULL
              Parent2GametesList <-getMotherGameteInfo(Parent2Genotype)
              Parent2Gametes<- names(Parent2GametesList$prob)
              Parent2Probs<- Parent2GametesList$prob
              names(Parent2Probs)<-NULL
              thisParent2.sample<- sample(Parent2Gametes,1,replace=TRUE,Parent2Probs)
              thisParent2.sample<- t(matrix(unlist(strsplit(thisParent2.sample," ")),ploidy/2,))
              # New individual for t+1
              newProgenyOne <- cbind(thisParent1.sample,thisParent2.sample)
              # Mutation
              newProgeny[1,locusRangeNewProgeny] <- mutation(newProgenyOne,thisLocus,mutationRate,mutationModel,
                                                            locusBounds,allelicStates)
          }

          ########################
          # Send to dataframe T2 #
          ########################

          # Migration odds for n1 generations
          if (thisGen <= n1)  {
              MigOdds <- 1
              migOccur <- sample(1:MigOdds,1,replace=TRUE)
              # send to another population
              # randomly pick across pops other than the curent one that are not full
              Pops <- (1:numPopsSims)
              otherPops <- Pops[Pops!=thisRandomPop]
              otherPopsNotFull <- statusT2["Full",otherPops]!=1

              if (all(otherPopsNotFull == FALSE)) {
                  # Increase Tally
                  statusT2["Tally",thisRandomPop] <- statusT2["Tally",thisRandomPop]+1
                  # Check if this makes it now at max, if so, put flag=1
                  if (as.logical(statusT2["Tally",thisRandomPop] == popSizeSims)) {
                      statusT2["Full",thisRandomPop] <- 1
                  }
                  # send data to T2 matrix
                  minValue <-(thisRandomPop*popSizeSims)-popSizeSims+1
                  maxValue <- thisRandomPop*popSizeSims
                  rangeValues <- minValue:maxValue
                  currentPos <- as.vector(statusT2["Tally",thisRandomPop])
                  thisProgeny <- rangeValues[currentPos]
                  matrixT2[thisProgeny,4:(length(lociLabels)+3)] <- newProgeny
              }
              if (any(otherPopsNotFull != FALSE)) {
                  remainPops <- otherPops[otherPopsNotFull]
                  thatRandomPop <- sample(remainPops,1,replace=TRUE)
                  if (as.vector(statusT2["Full",thatRandomPop] != 1)) {
                      # Increase Tally
                      statusT2["Tally",thatRandomPop] <- statusT2["Tally",thatRandomPop]+1
                      # Check if this makes it now at max, if so, put flag=1
                      if (as.logical(statusT2["Tally",thatRandomPop]==popSizeSims)) {
                          statusT2["Full",thatRandomPop] <- 1
                      }
                      # send data to T2 matrix
                      minValue <- (thatRandomPop*popSizeSims)-popSizeSims+1
                      maxValue <- thatRandomPop*popSizeSims
                      rangeValues <- minValue:maxValue
                      currentPos <- as.vector(statusT2["Tally",thatRandomPop])
                      thisProgeny <- rangeValues[currentPos]
                      matrixT2[thisProgeny,4:(length(lociLabels)+3)] <- newProgeny
                  }
              }
          } # close thisgen <n1

          if (thisGen > n1)  {
              # Migration odds for > n2 generations
              MigOdds <- as.integer(1/migration)
              migOccur <- sample(1:MigOdds,1,replace=TRUE)
              #############
              # Migration #
              #############
              if (migOccur == 1) {
                  # send to another population
                  # randomly pick across pops other than the curent one that are not full
                  Pops <- (1:numPopsSims)
                  otherPops <- Pops[Pops!=thisRandomPop]
                  otherPopsNotFull <- statusT2["Full",otherPops]!=1

                  if (all(otherPopsNotFull==FALSE)) {
                      # Increase Tally
                      statusT2["Tally",thisRandomPop] <- statusT2["Tally",thisRandomPop]+1
                      # Check if this makes it now at max, if so, put flag=1
                      if (as.logical(statusT2["Tally",thisRandomPop]==popSizeSims)) {
                          statusT2["Full",thisRandomPop] <- 1
                      }
                      # send data to T2 matrix
                      minValue <- (thisRandomPop*popSizeSims)- popSizeSims+1
                      maxValue <- thisRandomPop*popSizeSims
                      rangeValues <- minValue:maxValue
                      currentPos <- as.vector(statusT2["Tally",thisRandomPop])
                      thisProgeny <- rangeValues[currentPos]
                      matrixT2[thisProgeny,4:(length(lociLabels)+3)] <- newProgeny
                  }
                  if (any(otherPopsNotFull!=FALSE)) {
                      remainPops<-otherPops[otherPopsNotFull]
                      thatRandomPop<-sample(remainPops,1,replace=TRUE)
                      if (as.vector(statusT2["Full",thatRandomPop] != 1)) {
                          # Increase Tally
                          statusT2["Tally",thatRandomPop] <- statusT2["Tally",thatRandomPop]+1
                          # Check if this makes it now at max, if so, put flag=1
                          if (as.logical(statusT2["Tally",thatRandomPop] == popSizeSims)) {
                              statusT2["Full",thatRandomPop] <- 1
                          }
                          # send data to T2 matrix
                          minValue <-(thatRandomPop*popSizeSims) - popSizeSims+1
                          maxValue <- thatRandomPop*popSizeSims
                          rangeValues <- minValue:maxValue
                          currentPos <- as.vector(statusT2["Tally",thatRandomPop])
                          thisProgeny <- rangeValues[currentPos]
                          matrixT2[thisProgeny,4:(length(lociLabels)+3)] <- newProgeny
                       }
                  }
              } #end IF migration

              ################
              # NO Migration #
              ################
              if (migOccur > 1) {
                  # send to same population (if not full)
                  if (as.numeric(statusT2["Full", as.numeric(thisRandomPop)]) != 1) {
                      # Increase tally
                      statusT2["Tally",thisRandomPop] <- statusT2["Tally",thisRandomPop]+1
                      # If now reached max, put full flag (=1)
                      if (as.logical(statusT2["Tally",thisRandomPop] == popSizeSims)) {
                          statusT2["Full",thisRandomPop] <- 1
                      }
                      # send data to T2 matrix
                      minValue <- (thisRandomPop*popSizeSims)-popSizeSims+1
                      maxValue <- thisRandomPop*popSizeSims
                      rangeValues <- minValue:maxValue
                      currentPos <- as.vector(statusT2["Tally",thisRandomPop])
                      thisProgeny <- rangeValues[currentPos]
                      matrixT2[thisProgeny,4:(length(lociLabels)+3)] <- newProgeny
                    }
                  } # close home pop
              } # close thisgen >n1  (second phase migration rates)

           } # close while  (if matrix is empty)
         matrixT1 <-matrixT2
         } # close generation loops

      ##############################
      # extracting no. of alleles  #
      ##############################
      popGenStatsMatrix<-matrix(0,numPopsSims,lociSim)
      rownames(popGenStatsMatrix)<-1:numPopsSims
      for (thisPopn in 1:numPopsSims) {
          # thisPopn <- 1
          for (thisLocus in 1:lociSim) {
              # thisLocus <- 1
              locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
              thisPopnThisLocus <- matrixT1[matrixT1[,"pop"] == thisPopn,locusRange]
              popnAlleles <- unique(sort(strip.NAs(unlist(thisPopnThisLocus))))
              numAlleles <- length(popnAlleles)
              popGenStatsMatrix[thisPopn,thisLocus] <- numAlleles
          }
      }
      meanAlleles<-apply(popGenStatsMatrix,1,mean)
      SDAlleles<-apply(popGenStatsMatrix,1,sd)
      maxNumAlleles<-apply(popGenStatsMatrix,1,max)
      minNumAlleles<-apply(popGenStatsMatrix,1,min)
      # Send data to list
      resultsForwardSim <- list()
      resultsForwardSim[["SimData"]] <- list()
      resultsForwardSim[["LocusBounds"]] <- list()
      resultsForwardSim[["Allelic States"]] <- list()
      resultsForwardSim[["Mean No. Alleles"]] <- list()
      resultsForwardSim[["Standard deviation alleles number"]] <- list()
      resultsForwardSim[["Max No. Alleles"]] <- list()
      resultsForwardSim[["Min No. Alleles"]] <- list()
      resultsForwardSim[["Generations"]] <- list()
      resultsForwardSim[["SimData"]] <- matrixT1
      resultsForwardSim[["LocusBounds"]] <- locusBounds
      resultsForwardSim[["Allelic States"]] <- allelicStates
      resultsForwardSim[["Mean No. Alleles"]] <- meanAlleles
      resultsForwardSim[["Standard deviation alleles number"]] <- SDAlleles
      resultsForwardSim[["Max No. Alleles"]] <- maxNumAlleles
      resultsForwardSim[["Min No. Alleles"]] <- minNumAlleles
      resultsForwardSim[["Generations"]] <- n1 + n2
      return(resultsForwardSim)
      } # close function



















