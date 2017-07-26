
mutation<- function(newProgeny,thisLocus,mutationRate,mutationModel,locusBounds,allelicStates) {

    #############
    # Mutations #
    #############
    for (thisAllele in 1:ncol(newProgeny))  {
        #thisAllele<-1
        mutOdds<-1/mutationRate
        mutOccur<-sample(1:mutOdds,1,replace=TRUE)
        if (mutOccur == 1) {
            ########################
            # SMM mutational model #
            ########################
            if (mutationModel == "SMM") {
                odds<-runif(1,0,1)
                # Increase in bp
                if (odds > 0.5 ) {
                    if (as.numeric(newProgeny[thisAllele])+1 < as.numeric(locusBounds[thisLocus,"max"]))  {
                        newProgeny[thisAllele]<-as.numeric(newProgeny[thisAllele])+1
                    }
                    if (as.numeric(newProgeny[thisAllele])+1 > as.numeric(locusBounds[thisLocus,"max"]))  {
                        newProgeny[thisAllele]<-newProgeny[thisAllele]
                    }
                }
                # Decrease in bp
                if (odds < 0.5 ) {
                    if (as.numeric(newProgeny[thisAllele])-1 > as.numeric(locusBounds[thisLocus,"min"]))  {
                        newProgeny[thisAllele]<-as.numeric(newProgeny[thisAllele])-1
                    }
                    if (as.numeric(newProgeny[thisAllele])-1 < as.numeric(locusBounds[thisLocus,"min"]))  {
                        newProgeny[thisAllele]<-newProgeny[thisAllele]
                    }
                }
            } # close SMM
            ########################
            # IAM mutational model #
            ########################
            if (mutationModel == "IAM") {
                minchange<-as.numeric(locusBounds[thisLocus,"min"])
                maxchange<-as.numeric(locusBounds[thisLocus,"max"])
                newProgeny[thisAllele]<-sample(minchange:maxchange,1,replace=TRUE)
            } # close IAM
            
            ########################
            # KAM mutational model #
            ########################
            if (mutationModel == "KAM") {
                possibleChange <- allelicStates[thisLocus,]
                newProgeny[thisAllele] <- sample(possibleChange,1,replace=TRUE)
            } # close KAM
        } # close odds of mutation
   }  # close thisAllele
 return(newProgeny)
}
