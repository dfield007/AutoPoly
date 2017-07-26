#############################
#       F'ST                #
#    genotypes/phenotypes   #
#     Author: David Field   #
#       Date: 1/2/2010      #
#############################

Fst<- function(inData,numLoci,ploidy,samples)  {

  cat("\n ********************************************** \n")
  cat(" *          PolyAssign (1/1/2010)             * \n")
  cat(" *   Authors: David Field, Alec Zwart,        * \n")
  cat(" *            Andrew Young, Linda Broadhurst  * \n")
  cat(" ********************************************** \n")
  cat("\n ** Calculating Pair-wise Fst.... **\n")
  flush.console()
  strip.NAs <- function(vv) {
      return(vv[!is.na(vv)])
  }
  
  ##################################
  #  Setup objects to receive data #
  ##################################
   # inData<-simData1

  allPopns <- as.vector(sort(unique(inData$pop)))
  numPopns <- length(allPopns)
  ##Remove seedlings, leaving adults and mothers
  inDataAdults <- subset(inData,subset=is.na(mother))
  numAdults<-nrow(inDataAdults)

  ###################
  # Fis Calculation #
  ###################
  # Fis - all populations
  FisMatrix<-matrix(0,numPopns,numLoci)
  rownames(FisMatrix)<-allPopns
  colnames(FisMatrix)<-c(1:numLoci)
  AdultSubSample<-list()
  
  #################################
  ## Begin loop over populations ##
  #################################
  for (thisPopn in allPopns) {
       # test  thisPopn<-"GA"
       thisPopnData<-inDataAdults[inDataAdults$pop==thisPopn,]
       if (nrow(thisPopnData) <= samples) {
          AdultSubSample[[thisPopn]]<-thisPopnData
       }
       if (nrow(thisPopnData) > samples) {
          thisPopnSample<-sample(1:nrow(thisPopnData),samples,replace=FALSE)
          thisPopnData<-thisPopnData[thisPopnSample,]
          AdultSubSample[[thisPopn]]<-thisPopnData                 
       }
       
       #########################
       # Pairwise comparisons  #
       #########################
       Adults<-as.vector(rownames(AdultSubSample[[thisPopn]]))
       numAdults<-length(Adults)
       numAdultCombinations<-(numAdults*(numAdults-1))/2
       AdultCombinations<-combinations(numAdults,2,Adults,repeats=FALSE)

       ##########################
       ## Begin loop over loci ##
       ##########################
       for (thisLocus in 1:numLoci) {
          # test thisLocus<-1
          # Current locus column range
          locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
          thisPopnDataTemp<-AdultSubSample[[thisPopn]]
          thisPopnDatathisLocus<-thisPopnDataTemp[,locusRange]
          temp<-matrix(0,numAdultCombinations,1)
          FisThisLocus<-cbind(AdultCombinations,temp)
        
          ###############################
          ## Loop over each adult pair ##
          ###############################
          for (thisAdultPair in 1:numAdultCombinations) {
              # test thisAdultPair<-91
              Adult1<-FisThisLocus[thisAdultPair,1]   
              Adult2<-FisThisLocus[thisAdultPair,2]
              Adult1Phenotype<-unique(strip.NAs(thisPopnDatathisLocus[Adult1,]))
              Adult2Phenotype<-unique(strip.NAs(thisPopnDatathisLocus[Adult2,]))
              # Go to next pair if either adult has no data
              if (length(Adult1Phenotype) < 1) {
                next
              }
              if (length(Adult2Phenotype) < 1) {
                next
              }
              # unique alleles in adult1
              adultAlleleMatch1 <- match(Adult1Phenotype, Adult2Phenotype)
              # unique alleles in adult2
              adultAlleleMatch2 <- match(Adult2Phenotype, Adult1Phenotype)
              # No of different bands
              Total<- c(adultAlleleMatch1,adultAlleleMatch2)
              FisThisLocus[thisAdultPair,3]<-sum(is.na(Total))
          } # end loop over adult combinations
          
          # Send the mean for current locus to FisMatrix
          FisMatrix[thisPopn,thisLocus]<-mean(as.numeric(FisThisLocus[,3]))
        } # end loop over loci
      } # end loop over populations
   
  ####################
  # Fit Calculations #
  ####################
      
  # Pairwise population comparisons
  noPopCombinations<-(numPopns*(numPopns-1))/2
  popCombos<-combinations(numPopns,2,allPopns,repeats=FALSE)
  popComboMatrix<-matrix(0,1,noPopCombinations)
  for (thisCombo in 1:nrow(popCombos)) {
      Combo<- popCombos[thisCombo,]
      popComboMatrix[1,thisCombo]<-paste(Combo,collapse=" ")
  }
 
  # Fit - all population pairs
  popComboMatrix<-t(popComboMatrix)
  FitMatrix<-matrix(0,nrow(popComboMatrix),numLoci)
  rownames(FitMatrix)<-popComboMatrix[,1]
  colnames(FitMatrix)<-c(1:numLoci)

  # Copy same matrix ready for Fst
  FstMatrix<-FitMatrix

  #############################################
  ## Begin loop over population combinations ##
  #############################################
  for (thisPopnCombo in 1:noPopCombinations) {
       # test thisPopnCombo<-1 
       Pop1<- popCombos[thisPopnCombo,1]
       Pop2<- popCombos[thisPopnCombo,2]
       # Pull out population data
       thisPopn1Data<-AdultSubSample[[Pop1]]
       thisPopn2Data<-AdultSubSample[[Pop2]]
       # Combine populations
       combinedPopData<-rbind(thisPopn1Data,thisPopn2Data)
       # test thisPopnCombo<-1 
          if (thisPopnCombo==as.integer((noPopCombinations)*0.1)) { 
              cat("\n         10%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.2)) { 
              cat("20%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.3)) { 
              cat("30%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.4)) { 
              cat("40%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.5)) { 
              cat("50%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.6)) { 
              cat("60%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.7)) { 
              cat("70%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*0.8)) { 
              cat("80%..") }    
          if (thisPopnCombo==as.integer((noPopCombinations)*0.9)) { 
              cat("90%..") }
          if (thisPopnCombo==as.integer((noPopCombinations)*1)) { 
              cat("Complete \n") }
          flush.console()    
     
       #########################
       # Pairwise comparisons  #
       #########################
       Adults<- rownames(combinedPopData)
       numAdults<-length(Adults)
       numAdultCombinations<-(numAdults*(numAdults-1))/2
       AdultCombinations<-combinations(numAdults,2,Adults,repeats=FALSE)
     
       ##########################
       ## Begin loop over loci ##
       ##########################
       for (thisLocus in 1:numLoci) {
          #test thisLocus<-1
          # Current locus column range
          locusRange <- 3 + (thisLocus-1)*ploidy + 1:ploidy
          combinedPopDataThisLocus<-combinedPopData[,locusRange]
          temp<-matrix(0,numAdultCombinations,1)
          FitThisLocus<-cbind(AdultCombinations,temp)
          
          ###############################
          ## Loop over each adult pair ##
          ###############################
          for (thisAdultPair in 1:numAdultCombinations) {
              # test thisAdultPair<-1
              Adult1<-FitThisLocus[thisAdultPair,1]   
              Adult2<-FitThisLocus[thisAdultPair,2]
              Adult1Phenotype<-unique(strip.NAs(combinedPopDataThisLocus[Adult1,]))
              Adult2Phenotype<-unique(strip.NAs(combinedPopDataThisLocus[Adult2,]))
              # Go to next pair if either adult has no data
              if (length(Adult1Phenotype) < 1) {
                next
              }
              if (length(Adult2Phenotype) < 1) {
                next
              }
              # unique alleles in adult1
              adultAlleleMatch1 <- match(Adult1Phenotype, Adult2Phenotype)
              # unique alleles in adult2
              adultAlleleMatch2 <- match(Adult2Phenotype, Adult1Phenotype)
              # No of different bands
              Total<- c(adultAlleleMatch1,adultAlleleMatch2)
              FitThisLocus[thisAdultPair,3]<-sum(is.na(Total))
          } # end loop over adult combinations
          
          # Send the mean for current locus to FisMatrix
          FitMatrix[thisPopnCombo,thisLocus]<-mean(as.numeric(FitThisLocus[,3]))
        } # end loop over loci
      } # end loop over population pairs

  ####################
  # Fst Calculations #
  ####################

  FstTable <- matrix(0,numPopns,numPopns)
  colnames(FstTable)<- allPopns
  rownames(FstTable)<- allPopns
  
  multilocusFis<-apply(FisMatrix,1,mean)
  multilocusFit<-apply(FitMatrix,1,mean)
  Pairs<-names(multilocusFit)
  for (thisPopnCombo in Pairs) {
      # test thisPopnCombo<-"GA GE"
      Pops<- unlist(strsplit(thisPopnCombo," "))
      Pop1<-Pops[1]
      Pop2<-Pops[2]
      meanFisThisCombo<-as.vector((multilocusFis[Pop1]+multilocusFis[Pop2])/2)
      FitThisCombo<-as.vector(multilocusFit[thisPopnCombo])
      # Fst value for this pairwise comparison
      FstthisCombo<-(FitThisCombo-meanFisThisCombo)/FitThisCombo
      FstTable[Pop2,Pop1]<-round(FstthisCombo,3)
  }
  FstTable[lower.tri(FstTable)==FALSE]<-"*"
  FstTable<- as.data.frame(FstTable)
  return(FstTable)
}




      
