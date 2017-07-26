# Import a simulated data set for power testing
importDataFile<- function(filename) {
                     DataSet <- read.csv(filename,header=TRUE,quote="",strip.white = TRUE, na.strings=c("*",""),fill=FALSE)
                     DataSet<- as.data.frame(DataSet)
                     return(DataSet)
                     }
