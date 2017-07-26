######################################
#                                    #
# 8 Functions to obtain the probs    #
#  of gametes from a given phenotype #
#  or genotype, for double reduction #
#   or no double reduction           #
#                                    #
#      Author: David Field           #
#      Date: 1/2/2010                #
#                                    #
######################################

###############################################################
getMotherGameteInfoHexPhenoMin <- function (phenotype) {
##################################################
# 1. Hexaploid, marker="phenotype", DRR="min"    #
##################################################
  switch(length(phenotype),
         {## Monoallele
           gametes.table <- c("aaa")
           probs.table <- c(1.0)
         },
         {## Bialleles
           gametes.table <- c("aaa", "aab", "bbb", "abb")
           probs.table <- c(0.15, 0.35, 0.15, 0.35)
         },
         {## Trialleles
           gametes.table <- c("aaa", "aab", "aac", "bbb", "abb",
                              "bbc", "ccc", "acc", "bcc", "abc")  
           probs.table <- c(0.03, 0.105, 0.105, 0.03, 0.105, 0.105,
                            0.03, 0.105, 0.105, 0.28)
         },
         {## Quadrialleles
           gametes.table <-c("aaa", "aab", "aac", "aad", "bbb", "abb",
                             "bbc", "bbd", "ccc", "acc", "bcc", "ccd",
                             "ddd", "add", "bdd", "cdd", "abc", "abd",
                             "acd", "bcd")  
           probs.table <- c(0.005, 0.035, 0.035, 0.035, 0.005, 0.035,
                            0.035, 0.035, 0.005, 0.035, 0.035, 0.035,
                            0.005, 0.035, 0.035, 0.035, 0.14, 0.14,
                            0.14, 0.14)
         },
         {## Pentalleles
           gametes.table <- c("aab", "aac", "aad", "aae", "abb",
                              "bbc", "bbd", "bbe", "acc", "bcc",
                              "ccd", "cce", "add", "bdd", "cdd",
                              "dde", "aee", "bee", "cee", "dee",
                              "abc", "abd", "abe", "acd", "ace",
                              "ade", "bcd", "bce", "bde", "cde")  
           probs.table <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
                            0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
                            0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.08,
                            0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
                            0.08, 0.08)
         },
         {## Hexalleles
           gametes.table <- c("abc", "abd", "abe", "abf", "acd",
                              "ace", "acf", "ade", "adf", "aef",
                              "bcd", "bce", "bcf", "bde", "bdf",
                              "bef", "cde", "cdf", "cef", "def") 
           probs.table <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                            0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,phenotype)  {
                            y <- match(thisVec,c("a","b","c","d","e","f"))
                            return(phenotype[y])
                          },
                          phenotype)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}
##################################################
# 2. Hexaploid, marker="phenotype", DRR="max"    #
##################################################
getMotherGameteInfoHexPhenoMax <- function (phenotype) {
  switch(length(phenotype),
         {## Monoallele
           gametes.table <- c("aaa")
           probs.table <- c(1.0)},
         {## Bialleles
           gametes.table <- c("aaa", "aab", "bba", "bbb")
           probs.table <- c(0.1818, 0.3182, 0.3182, 0.1818)},
         {## Trialleles
           gametes.table <- c("aaa", "aab", "aac", "abc", "bba",
                              "bbb", "bbc", "cca", "ccb", "ccc") 
           probs.table <- c(0.049, 0.108, 0.108, 0.204, 0.108, 0.049,
                            0.108, 0.108, 0.108, 0.049)
         },  
         {## Quadrialleles
           gametes.table <- c("aaa", "aab", "aac", "aad", "abc",
                              "abd", "acd", "bba", "bbb", "bbc",
                              "bbd", "bcd", "cca", "ccb", "ccc",
                              "ccd", "dda", "ddb", "ddc", "ddd") 
           probs.table <- c(0.015, 0.045, 0.045, 0.045, 0.102, 0.102,
                            0.102, 0.045, 0.015, 0.045, 0.045, 0.102,
                            0.044, 0.045, 0.015, 0.045, 0.045, 0.045,
                            0.045, 0.015)
         },  
         {## Pentalleles
           gametes.table <- c("aaa", "aab", "aac", "aad", "aae",
                              "abc", "abd", "abe", "acd", "ace",
                              "ade", "bba", "bbb", "bbc", "bbd",
                              "bbe", "bcd", "bce", "bde", "cca",
                              "ccb", "ccc", "ccd", "cce", "cde",
                              "dda", "ddb", "ddc", "ddd", "dde",
                              "eea", "eeb", "eec", "eed", "eee") 
           probs.table <- c(0.004, 0.020, 0.020, 0.020, 0.020, 0.058,
                            0.058, 0.058, 0.058, 0.058, 0.058, 0.020,
                            0.004, 0.020, 0.020, 0.020, 0.058, 0.058,
                            0.058, 0.020, 0.020, 0.004, 0.020, 0.020,
                            0.058, 0.020, 0.020, 0.020, 0.004, 0.020,
                            0.020, 0.020, 0.020, 0.020, 0.004) 
         },
         {## Hexalleles
           gametes.table <- c("aab", "aac", "aad", "aae", "aaf",
                              "abc", "abd", "abe", "abf", "acd",
                              "ace", "acf", "ade", "adf", "aef",
                              "bba", "bbc", "bbd", "bbe", "bbf",
                              "bcd", "bce", "bcf", "bde", "bdf",
                              "bef", "cca", "ccb", "ccd", "cce",
                              "ccf", "cde", "cdf", "cef", "dda",
                              "ddb", "ddc", "dde", "ddf", "def",
                              "eea", "eeb", "eec", "eed", "eef",
                              "ffa", "ffb", "ffc", "ffd", "ffe") 
           probs.table <- c(0.009090667, 0.009090667, 0.009090667,
                            0.009090667, 0.009090667, 0.036364000,
                            0.036364000, 0.036364000, 0.036364000,
                            0.036364000, 0.036364000, 0.036364000,
                            0.036364000, 0.036364000, 0.036364000,
                            0.009090667, 0.009090667, 0.009090667,
                            0.009090667, 0.009090667, 0.036364000,
                            0.036364000, 0.036364000, 0.036364000,
                            0.036364000, 0.036364000, 0.009090667,
                            0.009090667, 0.009090667, 0.009090667,
                            0.009090667, 0.036364000, 0.036364000,
                            0.036364000, 0.009090667, 0.009090667,
                            0.009090667, 0.009090667, 0.009090667,
                            0.036364000, 0.009090667, 0.009090667,
                            0.009090667, 0.009090667, 0.009090667,
                            0.009090667, 0.009090667, 0.009090667,
                            0.009090667, 0.009090667)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,phenotype)  {
                            y <- match(thisVec,c("a","b","c","d","e","f"))
                            return(phenotype[y])
                          },
                          phenotype)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}

########################################################
getMotherGameteInfoHexGenoMin <- function (motherGenotype) {
##################################################
# 3. Hexaploid, marker="genotype", DRR="min"     #
##################################################
 
  switch(length(unique(motherGenotype)),
         {## Monoallele
           thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
           gametes.table <- c("aaa")
           probs.table <- c(1.0)
         },
         {## Biallele
             #count table of each allele, sorted so most frequent comes first 
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             if (max(thisTable)==5) { #Biallele type I, if max no. of allele is 5 must be type I
                gametes.table <- c("aaa", "aab")
                probs.table <- c(0.5, 0.5)
             } else if (max(thisTable)==4) { #Biallele type II, if max no. of allele is 4 must be type II
                gametes.table <- c("aaa", "aab", "abb")
                probs.table <- c(0.2, 0.6, 0.2)
             } else if (max(thisTable)==3) { #Biallele type III, if max no. of allele is 3 must be type III
                gametes.table <- c("aaa", "aab", "bbb", "abb")
                probs.table <- c(0.05, 0.45, 0.05, 0.45)
             }   
         }, #end switch
         {## Triallele
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             if (max(thisTable)==4) { #Triallele type I if max no. of allele is 4 must be type I
                gametes.table <- c("aaa", "aab", "aac", "abc")  
                probs.table <- c(0.2, 0.3, 0.3, 0.2)
             } else if (max(thisTable)==3) { #Triallele type II, if max no. of allele is 3 must be type II
                gametes.table <- c("aaa", "aab", "aac", "abb","bbc","abc")  
                probs.table <- c(0.05, 0.3, 0.15, 0.15, 0.05, 0.3)
             } else if (max(thisTable)==2) { #Triallele type III, if max no. of allele is 2 must be type III
                gametes.table <- c("aab", "aac", "abb", "bbc", "acc", "bcc", "abc")  
                probs.table <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4)
             }
         },
         {## Quadriallele
             thisTable <- sort(table(motherGenotype),decreasing=TRUE)  #count table of each allele
             if (max(thisTable)==3) { #Quadriallele type I if max no. of allele is 3 must be type I
                  gametes.table <-c("aaa", "aab", "aac", "aad","abc", "abd", "acd", "bcd")  
                  probs.table <- c(0.05, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.05)
             } else if (max(thisTable)==2) { #Quadriallele type II if max no. of allele is 2 must be type II
                  gametes.table <-c("aab", "aac", "aad", "abb", "bbc", "bbd", "abc", "abd", "acd", "bcd")  
                  probs.table <- c(0.1, 0.05, 0.05, 0.1, 0.05, 0.05, 0.2, 0.2, 0.1, 0.1)
             }
         },
         {## Pentallele
           thisTable <- sort(table(motherGenotype),decreasing=TRUE)  #count table of each allele
           gametes.table <- c("aab", "aac", "aad", "aae", "abc", "abd", "abe", 
                              "acd", "ace", "ade", "bcd", "bce", "bde", "cde")  
           probs.table <- c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1,
                            0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05)
         },
         {## Hexalleles
           thisTable <- sort(table(motherGenotype),decreasing=TRUE)  #count table of each allele
           gametes.table <- c("abc", "abd", "abe", "abf", "acd",
                              "ace", "acf", "ade", "adf", "aef",
                              "bcd", "bce", "bcf", "bde", "bdf",
                              "bef", "cde", "cdf", "cef", "def") 
           probs.table <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                            0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  
  thisGenoAllelesOrder<-names(thisTable)
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,thisGenoAllelesOrder)  {
                          # gametes.table<-strsplit(gametes.table,"")
                          #thisVec<-gametes.table[[2]]
                            y <- match(thisVec,c("a","b","c","d","e","f"))
                            return(thisGenoAllelesOrder[y])
                          },
                          thisGenoAllelesOrder)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}

##########################################################
getMotherGameteInfoHexGenoMax <- function (motherGenotype) {
##################################################
# 4. Hexaploid, marker="genotype", DRR="max"     #
##################################################
  switch(length(unique(motherGenotype)),
         {## Monoallele
           thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
           gametes.table <- c("aaa")
           probs.table <- c(1.0)
         },
         {## Biallele
             #count table of each allele, sorted so most frequent comes first 
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             if (max(thisTable)==5) { #Biallele type I, if max no. of allele is 5 must be type I
                gametes.table <- c("aaa", "aab", "abb")
                probs.table <- c(0.546, 0.409, 0.045)
             } else if (max(thisTable)==4) { #Biallele type II, if max no. of allele is 4 must be type II
                gametes.table <- c("aaa", "aab", "bbb", "abb")
                probs.table <- c(0.255, 0.509, 0.018, 0.218)
             } else if (max(thisTable)==3) { #Biallele type III, if max no. of allele is 3 must be type III
                gametes.table <- c("aaa", "aab", "bbb", "abb")
                probs.table <- c(0.091, 0.409, 0.091, 0.409)
             }   
         }, #end switch
         {## Triallele
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             if (max(thisTable)==4) { #Triallele type I if max no. of allele is 4 must be type I
                gametes.table <- c("aaa", "aab", "aac", "abb", "bbc", "acc", "bcc", "abc")  
                probs.table <- c(0.255, 0.255, 0.255, 0.036, 0.009, 0.036, 0.009, 0.145)
             } else if (max(thisTable)==3) { #Triallele type II, if max no. of allele is 3 must be type II
                gametes.table <- c("aaa", "aab", "aac", "bbb", "abb", "bbc", "acc", "bcc", "abc")    
                probs.table <- c(0.091, 0.273, 0.136, 0.018, 0.164, 0.055, 0.027, 0.018, 0.218)
             } else if (max(thisTable)==2) { #Triallele type III, if max no. of allele is 2 must be type III
                gametes.table <- c("aaa", "aab", "aac", "bbb", "abb", "bbc", "ccc","acc", "bcc", "abc") 
                probs.table <- c(0.018, 0.109, 0.109, 0.018, 0.109, 0.109, 0.018, 0.109, 0.109, 0.292)
             }
         },
         {## Quadriallele
             thisTable <- sort(table(motherGenotype),decreasing=TRUE)  #count table of each allele
             if (max(thisTable)==3) { #Quadriallele type I if max no. of allele is 3 must be type I
                  gametes.table <-c("aaa", "aab", "aac", "aad", "abb", "bbc", "bbd", "acc", "bcc", 
                                    "ccd", "add", "bdd", "cdd", "abc", "abd", "acd", "bcd")  
                  probs.table <- c(0.0909, 0.1364, 0.1364, 0.1364, 0.0273, 0.0091, 0.0091, 0.0273, 0.0091, 
                                    0.0091, 0.0273, 0.0091, 0.0091, 0.1091, 0.1091, 0.1091, 0.0364)
             } else if (max(thisTable)==2) { #Quadriallele type II if max no. of allele is 2 must be type II
                  gametes.table <-c("aaa", "aab", "aac", "aad", "bbb", "abb", "bbc", "bbd", "acc", 
                                    "bcc", "ccd", "add", "bdd", "cdd", "abc", "abd", "acd", "bcd")  
                  probs.table <- c(0.0182, 0.1091, 0.0545, 0.0545, 0.0182, 0.1091, 0.0545, 0.0545, 0.0182, 
                                   0.0182, 0.0091, 0.0182, 0.0182, 0.0091, 0.1455, 0.1455, 0.0727, 0.0727)
             }
         },
         {## Pentallele
           thisTable <- sort(table(motherGenotype),decreasing=TRUE)  #count table of each allele
           gametes.table <- c("aaa", "aab", "aac", "aad", "aae", "abb",
                              "bbc", "bbd", "bbe", "acc", "bcc", "ccd", 
                              "cce", "add", "bdd", "cdd", "dde", "aee", 
                              "bee", "cee", "dee", "abc", "abd", "abe", 
                              "acd", "ace", "ade", "bcd", "bce", "bde", 
                              "cde")  
           probs.table <- c(0.018, 0.055, 0.055, 0.055, 0.055, 0.018, 0.009,
                            0.009, 0.009, 0.018, 0.009, 0.009, 0.009, 0.018,
                            0.009, 0.009, 0.009, 0.018, 0.009, 0.009, 0.009,
                            0.073, 0.073, 0.073, 0.073, 0.073, 0.073, 0.036,
                            0.036, 0.036, 0.036)
         },
         {## Hexalleles
           thisTable <- sort(table(motherGenotype),decreasing=TRUE)  #count table of each allele
           gametes.table <- c("aab", "aac", "aad", "aae", "aaf", "abb", "bbc", 
                              "bbd", "bbe", "bbf", "acc", "bcc", "ccd", "cce", 
                              "ccf", "add", "bdd", "cdd", "dde", "ddf", "aee",
                              "bee", "cee", "dee", "eef", "aff", "bff", "cff",
                              "dff", "eff", "abc", "abd", "abe", "abf", "acd",
                              "ace", "acf", "ade", "adf", "aef", "bcd", "bce", 
                              "bcf", "bde", "bdf", "bef", "cde", "cdf", "cef", 
                              "def") 
           probs.table <- c(0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009,   
                            0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009,
                            0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009,
                            0.036, 0.036, 0.036, 0.036, 0.036,0.036, 0.036, 0.036, 0.036, 0.036,
                            0.036, 0.036, 0.036, 0.036, 0.036,0.036, 0.036, 0.036, 0.036, 0.036)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  
  thisGenoAllelesOrder<-names(thisTable)
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,thisGenoAllelesOrder)  {
                          # gametes.table<-strsplit(gametes.table,"")
                          #thisVec<-gametes.table[[2]]
                            y <- match(thisVec,c("a","b","c","d","e","f"))
                            return(thisGenoAllelesOrder[y])
                          },
                          thisGenoAllelesOrder)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}


##########################################################
getMotherGameteInfoTetrPhenoMin <- function (phenotype) {
##################################################
# 5. Tetraploid, marker="phenotype", DRR="min"   #
##################################################
  switch(length(phenotype),
         {## Monoallele
           gametes.table <- c("aa")
           probs.table <- c(1.0)
         },
         {## Bialleles
           gametes.table <- c("aa","ab","bb")
           probs.table <- c(0.222,0.556,0.222)
         },
         {## Trialleles
           gametes.table <- c("aa", "ab", "ac", "bb", "bc", "cc")
           probs.table <- c(0.056, 0.222, 0.278, 0.111, 0.278, 0.056)
         },
         {## Quadrialleles
           gametes.table <- c("ab", "ac", "ad", "bc", "bd", "cd")
           probs.table <- c(0.167, 0.167, 0.167, 0.167, 0.167, 0.167)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,phenotype)  {
                            y <- match(thisVec,c("a","b","c","d"))
                            return(phenotype[y])
                          },
                          phenotype)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}


############################################################
getMotherGameteInfoTetrPhenoMax <- function (phenotype) {
  ##################################################
  # 6. Tetraploid, marker="phenotype", DRR="max"   #
  ##################################################

  switch(length(phenotype),
         {## Monoallele
           gametes.table <- c("aa")
           probs.table <- c(1.0)
         },
         {## Bialleles
           gametes.table <- c("aa","ab","bb") 
           probs.table <- c(0.2619,0.4762,0.2619)
         },
         {## Trialleles
           gametes.table <- c("aa", "ab", "ac", "bb", "bc", "cc")
           probs.table <- c(0.095, 0.238, 0.238, 0.095, 0.238, 0.095)
         },
         {## Quadrialleles
           gametes.table <- c("aa", "ab", "ac", "ad", "bb", "bc", "bd", "cc", "cd", "dd")
           probs.table <- c(0.036, 0.143, 0.143, 0.143, 0.036, 0.143, 0.143, 0.036, 0.143, 0.036)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,phenotype)  {
                            y <- match(thisVec,c("a","b","c","d"))
                            return(phenotype[y])
                          },
                          phenotype)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}

               #test     motherGenotype<-c("188","206")

##########################################################
getMotherGameteInfoTetrGenoMin <- function (motherGenotype) {
##################################################
# 7. Tetraploid, marker="genotype", DRR="min"    #
##################################################
  switch(length(unique(motherGenotype)),
         {## Monoallele
           thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
           gametes.table <- c("aa")
           probs.table <- c(1.0)
         },
         {## Bialleles
             #count table of each allele, sorted so most frequent comes first 
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             if (max(thisTable)==3) { #Biallele type I (simplex), if max no. of allele is 3
                gametes.table <- c("aa","ab")
                probs.table <- c(0.5,0.5)
             } else if (max(thisTable)==2) { #Biallele type II (duplex), if max no. of allele is 2
                gametes.table <- c("aa","ab","bb")
                probs.table <- c(0.167,0.667,0.167)
             }   
         }, #end switch
         {## Trialleles
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             gametes.table <- c("aa", "ab", "ac", "bc")
             probs.table <- c(0.16666, 0.333333, 0.333333, 0.16666)
         },
         {## Quadrialleles
           thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
           gametes.table <- c("ab", "ac", "ad", "bc", "bd", "cd")
           probs.table <- c(0.167, 0.167, 0.167, 0.167, 0.167, 0.167)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  thisGenoAllelesOrder<-names(thisTable)
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,thisGenoAllelesOrder)  {
                          # gametes.table<-strsplit(gametes.table,"")
                          #thisVec<-gametes.table[[2]]
                            y <- match(thisVec,c("a","b","c","d","e","f"))
                            return(thisGenoAllelesOrder[y])
                          },
                          thisGenoAllelesOrder)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}

##########################################################
getMotherGameteInfoTetrGenoMax <- function (motherGenotype) {
##################################################
# 8. Tetraploid, marker="genotype", DRR="max"    #
##################################################
  switch(length(unique(motherGenotype)),
         {## Monoallele
           thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
           gametes.table <- c("aa")
           probs.table <- c(1.0)
         },
         {## Bialleles
             #count table of each allele, sorted so most frequent comes first 
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             if (max(thisTable)==3) { #Biallele type I (simplex), if max no. of allele is 3
                gametes.table <- c("aa","ab","bb")
                probs.table <- c(0.536,0.429,0.035)
             } else if (max(thisTable)==2) { #Biallele type II (duplex), if max no. of allele is 2
                gametes.table <- c("aa","ab","bb")
                probs.table <- c(0.214,0.572,0.214)
             }   
         }, #end switch
         {## Trialleles
             thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
             gametes.table <- c("aa", "ab", "ac", "bc", "bb", "cc")
             probs.table <- c(0.214286, 0.285714, 0.285714, 0.142857, 0.03571, 0.03571)
         },
         {## Quadrialleles
           thisTable <- sort(table(motherGenotype),decreasing=TRUE) 
           gametes.table <- c("aa", "ab", "ac", "ad", "bc", "bd", "bb", "cc", "cd", "dd")
           probs.table <- c(0.036, 0.143, 0.143, 0.143, 0.143, 0.143, 0.036, 0.036, 0.143, 0.036)
         }
         )  ##End switch()
  ## Convert gametes.table from a vector of abc strings to a list of
  ## allele vectors
  thisGenoAllelesOrder<-names(thisTable)
  gametes.table <- lapply(strsplit(gametes.table,""),
                          function (thisVec,thisGenoAllelesOrder)  {
                          # gametes.table<-strsplit(gametes.table,"")
                          #thisVec<-gametes.table[[2]]
                            y <- match(thisVec,c("a","b","c","d","e","f"))
                            return(thisGenoAllelesOrder[y])
                          },
                          thisGenoAllelesOrder)
  ##Form the gamete names and assign these to gametes.table and
  ##probs.table 
  nms <- sapply(gametes.table,
                function(thisVec) {
                  paste(thisVec,collapse=" ")
                })
  names(gametes.table) <- nms
  names(probs.table) <- nms
  return(list(gametes=gametes.table, prob=probs.table))
##                    stringsAsFactors = FALSE))
}






















##There are more codes to come, that implement the more general cases
##of 0 < DRR < maximum for each ploidy.

##Note also that the above probability tables depend only on the
##number of alleles present in each possible gamete - can we
##abbreviate the above code then?  use some sort of combinatorial
##selection code to derive the set of possible gametes from the
##specified phenotype alleles, then assign the relevant probabilities
##based on the number of unique alleles in each Gamete?  Hmmm....
##The reason I am interested in this approach rather than the above,
##is that I beleive that DF's formulas for the general DRR case will
##depend only on the number of unique alleles in the potential
##Gamete.  Hence t'would be sensible to have code that assigns probs
##on the basis of the number of unique alleles in the gamete...

## the combinat package, function combn()does NOT appear to be the way to 
## generate the possible gametes, unfortuantely....

## require(combinat)
## combn(c("a","b","c"),3)
## combn(c("a","b","c"),2)
## ##Nope - that's not what I want...

## (Combns and Permns do not allow repetition of the elements...)

##Although there is the following:
##    combn(rep(c("a","b","c"),3),3)
##This gets closer, but Combn does not recognise that elements of the
##sample are repeated - hence we do not get true combinations (since,
##for example, aab, aba, and baa are all provided in the output.  To
##get what I want, the result of this command would need to be sorted
##by columns (the result is a matrix), then reduced to the unique
##columns...  We could do this, of course, but there must be an easier
##way to get what I want...

##Expand.grid is another way to get all possible combinations, but I
##have to repeat the sample as many times as the gamete length...
##E.g., for the input Phenotype c("a","b","c"), with gamete length 3,
##we use:
##  expand.grid(c("a","b","c"),c("a","b","c"),c("a","b","c"))
##However, this still has the unwanted repetitions which need to be
##removed. 

## What do I need?  I need to generate the unique, unordered
## selections (with repetition) of m objects from k.  There are
## formulas to determing the NUMBER of such selections, but I want a
## simple algorithm (or better, an R function) to generate all of the
## selections themselves...

## I don't think this is going to happen, actually - better to just
## have the tables of selections ready as needed.  Still, we can
## combine each such table with a vector of "numbers of indices" that
## also is used to index a shorter vector of probabilities.  For
## example, in getMotherGameteInfoHexMin() the trialleles case is:
##            gametes.table <- c("aaa", "aab", "aac", "bbb", "abb",
##                              "bbc", "ccc", "acc", "bcc", "abc")  
##           probs.table <- c(0.03, 0.105, 0.105, 0.03, 0.105, 0.105,
##                            0.03, 0.105, 0.105, 0.28)
##We could replace these lines by:
##            gametes.table <- c("aaa", "aab", "aac", "bbb", "abb",
##                              "bbc", "ccc", "acc", "bcc", "abc")  
##           probs.table <- c(0.03, 0.105, 0.28)
##           probs.ind <- c(1,2,2,1,2,2,1,2,2,3)
## hence:
##          probs.table[probs.ind]
## will give the appropriate probability vector to match the vector of
## gametes... Probs.ind can be specified explicitely as above, but we
## could generate it using some R code to count the number of unique
## letters in each entry of gametes.table...


