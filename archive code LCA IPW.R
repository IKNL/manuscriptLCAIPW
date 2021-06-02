###############################################################################################
# last edit: 02/06/2021

# accompanying code to the manuscript:
# "A new three-step method for using inverse propensity weighting with latent class analysis"

# F.J. Clouth, MSc.
# Department of Methodology and Statistics, Tilburg University
# f.j.clouth@tilburguniversity.edu

library(foreign)
library(rio)
library(dplyr)
library(ggplot2)
library(ggsci)
library(survey)
library(tableone)

setwd("adjust")
LG <- "adjust/LatentGOLD6.0/lg60.exe"
# latent Gold can be downloaded at https://www.statisticalinnovations.com/latent-gold-6-0/
# latent Gold syntax is generate and run from this code

###############################################################################################
# Simulation study
# generate latent Gold syntax

# simulate 1000 data sets per condition
GenData <- function(syntaxName, infile, outfile, N, B, G){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", infile,"'
                                                    
                                                    model
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=", N," tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors profile probmeans=posterior
                                                    loadings bivariateresiduals estimatedvalues=model;
                                                    outfile '", outfile,"' simulation=1000;
                                                    variables
                                                    caseweight n1000;
                                                    dependent Z 2, Y1 2, Y2 2, Y3 2, Y4 2, Y5 2, Y6 2;
                                                    independent C1, C2;
                                                    latent
                                                    Zlat nominal 2, Cluster nominal 3;
                                                    equations
                                                    Zlat <- 1 + C1 + C2;
                                                    Cluster <- 1 + C1 + C2 + Zlat;
                                                    Z <- (w~wei) 1 | Zlat; 
                                                    Y1 <- 1 | Cluster;
                                                    Y2 <- 1 | Cluster;
                                                    Y3 <- 1 | Cluster;
                                                    Y4 <- 1 | Cluster;
                                                    Y5 <- 1 | Cluster;
                                                    Y6 <- 1 | Cluster;
                                                    
                                                    w = {1 0 0 1};
                                                    
                                                    {0 1 ", B," .5 .5 1 1 1 1 1 ", G,"
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# first step (measurement model, only needed in the three-step approach)
Step1Syntax <- function(syntaxName, Step1Infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", Step1Infile,"'
                                                    
                                                    
                                                    model
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors profile bivariateresiduals estimatedvalues=model;
                                                    outfile  '", outfile,"' classification=posterior keep simulation_., C1, C2, Z, IPW;
                                                    variables
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1;
                                                    Y1 - Y6 <- 1 | Cluster;
                                                    {0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# third step using IPW
Step3Syntax <- function(syntaxName, Step3Infile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", Step3Infile,"'
                                                    model
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    step3 proportional ml;  
                                                    output      
                                                    parameters=first betaopts=wl standarderrors=robust profile estimatedvalues=model marginaleffects append='", outfile,"';
                                                    variables
                                                    samplingweight IPW ipw;
                                                    independent Z nominal;
                                                    latent Cluster nominal posterior = ( Cluster#1 Cluster#2 Cluster#3 ) ;
                                                    equations
                                                    Cluster <- 1 + Z;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# one-step analysis using IPW
LanzaSyntax <- function(syntaxName, LanzaInfile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", LanzaInfile,"'
                                                    model
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors=robust profile bivariateresiduals estimatedvalues=model marginaleffects append='", outfile,"';
                                                    variables
                                                    samplingweight IPW rescale;
                                                    independent Z nominal;
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + Z;
                                                    Y1 - Y6 <- 1 | Cluster;
                                                    {0 0 0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# one-step analysis using regression-adjustment
NoWeightSyntax <- function(syntaxName, NoWeightInfile, outfile){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile '", NoWeightInfile,"'
                                                    model  
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors=robust profile probmeans=posterior
                                                    loadings bivariateresiduals estimatedvalues=model marginaleffects append='", outfile,"';
                                                    variables
                                                    dependent Y1, Y2, Y3, Y4, Y5, Y6;
                                                    independent C1, C2, Z nominal;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + C1 + C2 + Z;
                                                    Y1 - Y6 <- 1 | Cluster;
                                                    {0 0 0 0 0 0 0 0
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361    1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    1.386294361   -1.386294361    -1.386294361
                                                    }
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# simulation conditions
B <- c(1, 2, 3)
G <- c(1, 2, 3)
N <- c(500, 1000, 2500)


# run syntax to generate 1000 data sets per condition (stored in one big data set per condition)
for(b in B) {
  for(g in G) {
    for(n in N) {
      GenData(syntaxName = paste0("GenData_B",b, "G",g, "N",n), infile = paste0("example",n, ".sav") , outfile = paste0("simB",b, "G",g, "N",n, ".sav"), N = n, B = b, G = g)
      shell(paste0(LG, " ", "GenData_B",b, "G",g, "N",n, ".lgs", " ", "/b"))
      print(c(b, g, n))
    }
  }
}


# separate the 1000 data sets per condition and estimate IPWs
for(b in B) {
  for(g in G) {
    for(n in N) {
      sim <- read.spss(paste0("simB",b, "G",g, "N",n, ".sav"), to.data.frame = T)
      sim$Z[sim$Z==2] <- 0
      for(r in 1:1000) {
        sim_r <- sim %>%
          filter(simulation_. == r)
        
        
        MyGLM <- glm(Z ~ C1 + C2, family = "binomial", data = sim_r)
        PS_est <- MyGLM$fitted.values
        
        sim_r$PS <- PS_est
        
        sim_r$IPW <- ifelse(sim_r$Z == 1, 1/sim_r$PS, 1/(1-sim_r$PS))
        
        
        export(sim_r, paste0("simB",b, "G",g, "N",n, "_with_IPW_Rep",r, ".sav"))

        print(c(b, g, n, r))
      }
    }
  }
}


# run first step (measurement model)
for(b in B) {
  for(g in G) {
    for(n in N) {
      for(r in 1:1000) {
        Step1Syntax(syntaxName = paste0("Step1Syntax_B",b, "G",g, "N",n, "Rep",r), Step1Infile = paste0("simB",b, "G",g, "N",n, "_with_IPW_Rep",r, ".sav") , outfile = paste0("classificationB",b, "G",g, "N",n, "Rep",r, ".sav"))

        shell(paste0(LG, " ", "Step1Syntax_B",b, "G",g, "N",n, "Rep",r,".lgs", " ", "/b"))
        
        print(c(b, g, n, r))
      }
    }
  }
}


# run the three different methods
for(b in B) {
  for(g in G) {
    for(n in N) {
      for(r in 1:1000) {
        Step3Syntax(syntaxName = paste0("Step3Syntax_B",b, "G",g, "N",n, "Rep",r), Step3Infile = paste0("classificationB",b, "G",g, "N",n, "Rep",r, ".sav") , outfile = paste0("Step3Output_B",b, "G",g, "N",n, "Rep",r, ".csv"))
        LanzaSyntax(syntaxName = paste0("LanzaSyntax_B",b, "G",g, "N",n, "Rep",r), LanzaInfile = paste0("simB",b, "G",g, "N",n, "_with_IPW_Rep",r, ".sav") , outfile = paste0("LanzaOutput_B",b, "G",g, "N",n, "Rep",r, ".csv"))
        NoWeightSyntax(syntaxName = paste0("NoweightSyntax_B",b, "G",g, "N",n, "Rep",r), NoWeightInfile = paste0("simB",b, "G",g, "N",n, "_with_IPW_Rep",r, ".sav") , outfile = paste0("NoWeightOutput_B",b, "G",g, "N",n, "Rep",r, ".csv"))
        
        shell(paste0(LG, " ", "Step3Syntax_B",b, "G",g, "N",n, "Rep",r,".lgs", " ", "/b"))
        shell(paste0(LG, " ", "LanzaSyntax_B",b, "G",g, "N",n, "Rep",r,".lgs", " ", "/b"))
        shell(paste0(LG, " ", "NoweightSyntax_B",b, "G",g, "N",n, "Rep",r,".lgs", " ", "/b"))      
        
        print(c(b, g, n, r))
      }
    }
  }
}



# read results of the simulation study and calculate bias of the ATE, sd of the ATE, and bias of the SE

UnadjustedATE <- array(data = NA, dim = c(3, 3, 3, 1000))
UnadjustedATESE <- array(data = NA, dim = c(3, 3, 3, 1000))
LanzaATE <- array(data = NA, dim = c(3, 3, 3, 1000))
LanzaATESE <- array(data = NA, dim = c(3, 3, 3, 1000))
Step3ATE <- array(data = NA, dim = c(3, 3, 3, 1000))
Step3ATESE <- array(data = NA, dim = c(3, 3, 3, 1000))


SampleN <- 500
for(b in 1:3) {
  for(g in 1:3) {
    for(r in 1:1000) {
      UnadjustedResults <- read.table(paste0("NoWeightOutput_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      UnadjustedEstVal <- as.character(UnadjustedResults[7, ])
      UnadjustedEstVal <- strsplit(UnadjustedEstVal, ",")
      UnadjustedATE[b, g, 1, r] <- matrix(as.numeric(UnadjustedEstVal[[1]][17]), 1, 1, byrow = T)
      UnadjustedEstValSE <- as.character(UnadjustedResults[11, ])
      UnadjustedEstValSE <- strsplit(UnadjustedEstValSE, ",")
      UnadjustedATESE[b, g, 1, r] <- matrix(as.numeric(UnadjustedEstValSE[[1]][17]), 1, 1, byrow = T)
      
      LanzaResults <- read.table(paste0("LanzaOutput_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      LanzaEstVal <- as.character(LanzaResults[7, ])
      LanzaEstVal <- strsplit(LanzaEstVal, ",")
      LanzaATE[b, g, 1, r] <- matrix(as.numeric(LanzaEstVal[[1]][11]), 1, 1, byrow = T)
      LanzaEstValSE <- as.character(LanzaResults[11, ])
      LanzaEstValSE <- strsplit(LanzaEstValSE, ",")
      LanzaATESE[b, g, 1, r] <- matrix(as.numeric(LanzaEstValSE[[1]][11]), 1, 1, byrow = T)
      
      Step3Results <- read.table(paste0("Step3Output_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      Step3EstVal <- as.character(Step3Results[6, ])
      Step3EstVal <- strsplit(Step3EstVal, ",")
      Step3ATE[b, g, 1, r] <- matrix(as.numeric(Step3EstVal[[1]][11]), 1, 1, byrow = T)
      Step3EstValSE <- as.character(Step3Results[10, ])
      Step3EstValSE <- strsplit(Step3EstValSE, ",")
      Step3ATESE[b, g, 1, r] <- matrix(as.numeric(Step3EstValSE[[1]][11]), 1, 1, byrow = T)
      
      print(c(b, g, r))
    }
  }
}

SampleN <- 1000
for(b in 1:3) {
  for(g in 1:3) {
    for(r in 1:1000) {
      UnadjustedResults <- read.table(paste0("NoWeightOutput_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      UnadjustedEstVal <- as.character(UnadjustedResults[7, ])
      UnadjustedEstVal <- strsplit(UnadjustedEstVal, ",")
      UnadjustedATE[b, g, 2, r] <- matrix(as.numeric(UnadjustedEstVal[[1]][17]), 1, 1, byrow = T)
      UnadjustedEstValSE <- as.character(UnadjustedResults[11, ])
      UnadjustedEstValSE <- strsplit(UnadjustedEstValSE, ",")
      UnadjustedATESE[b, g, 2, r] <- matrix(as.numeric(UnadjustedEstValSE[[1]][17]), 1, 1, byrow = T)
      
      LanzaResults <- read.table(paste0("LanzaOutput_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      LanzaEstVal <- as.character(LanzaResults[7, ])
      LanzaEstVal <- strsplit(LanzaEstVal, ",")
      LanzaATE[b, g, 2, r] <- matrix(as.numeric(LanzaEstVal[[1]][11]), 1, 1, byrow = T)
      LanzaEstValSE <- as.character(LanzaResults[11, ])
      LanzaEstValSE <- strsplit(LanzaEstValSE, ",")
      LanzaATESE[b, g, 2, r] <- matrix(as.numeric(LanzaEstValSE[[1]][11]), 1, 1, byrow = T)
      
      Step3Results <- read.table(paste0("Step3Output_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      Step3EstVal <- as.character(Step3Results[6, ])
      Step3EstVal <- strsplit(Step3EstVal, ",")
      Step3ATE[b, g, 2, r] <- matrix(as.numeric(Step3EstVal[[1]][11]), 1, 1, byrow = T)
      Step3EstValSE <- as.character(Step3Results[10, ])
      Step3EstValSE <- strsplit(Step3EstValSE, ",")
      Step3ATESE[b, g, 2, r] <- matrix(as.numeric(Step3EstValSE[[1]][11]), 1, 1, byrow = T)
      
      print(c(b, g, r))
    }
  }
}


SampleN <- 2500
for(b in 1:3) {
  for(g in 1:3) {
    for(r in 1:1000) {
      UnadjustedResults <- read.table(paste0("NoWeightOutput_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      UnadjustedEstVal <- as.character(UnadjustedResults[7, ])
      UnadjustedEstVal <- strsplit(UnadjustedEstVal, ",")
      UnadjustedATE[b, g, 3, r] <- matrix(as.numeric(UnadjustedEstVal[[1]][17]), 1, 1, byrow = T)
      UnadjustedEstValSE <- as.character(UnadjustedResults[11, ])
      UnadjustedEstValSE <- strsplit(UnadjustedEstValSE, ",")
      UnadjustedATESE[b, g, 3, r] <- matrix(as.numeric(UnadjustedEstValSE[[1]][17]), 1, 1, byrow = T)
      
      LanzaResults <- read.table(paste0("LanzaOutput_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      LanzaEstVal <- as.character(LanzaResults[7, ])
      LanzaEstVal <- strsplit(LanzaEstVal, ",")
      LanzaATE[b, g, 3, r] <- matrix(as.numeric(LanzaEstVal[[1]][11]), 1, 1, byrow = T)
      LanzaEstValSE <- as.character(LanzaResults[11, ])
      LanzaEstValSE <- strsplit(LanzaEstValSE, ",")
      LanzaATESE[b, g, 3, r] <- matrix(as.numeric(LanzaEstValSE[[1]][11]), 1, 1, byrow = T)
      
      Step3Results <- read.table(paste0("Step3Output_B",b, "G",g, "N",SampleN, "Rep",r, ".csv"))
      Step3EstVal <- as.character(Step3Results[6, ])
      Step3EstVal <- strsplit(Step3EstVal, ",")
      Step3ATE[b, g, 3, r] <- matrix(as.numeric(Step3EstVal[[1]][11]), 1, 1, byrow = T)
      Step3EstValSE <- as.character(Step3Results[10, ])
      Step3EstValSE <- strsplit(Step3EstValSE, ",")
      Step3ATESE[b, g, 3, r] <- matrix(as.numeric(Step3EstValSE[[1]][11]), 1, 1, byrow = T)
      
      print(c(b, g, r))
    }
  }
}


TrueATE <- c(0.0664870, 0.3019906, 0.4824338)
biasAdj <- array(data = NA, dim = c(3, 3, 3, 1000))
biasAdjSE <- array(data = NA, dim = c(3, 3, 3, 1000))
biasLanza <- array(data = NA, dim = c(3, 3, 3, 1000))
biasLanzaSE <- array(data = NA, dim = c(3, 3, 3, 1000))
biasStep3 <- array(data = NA, dim = c(3, 3, 3, 1000))
biasStep3SE <- array(data = NA, dim = c(3, 3, 3, 1000))
AdjSD <- array(data = NA, dim = c(3, 3, 3))
LanzaSD <- array(data = NA, dim = c(3, 3, 3))
Step3SD <- array(data = NA, dim = c(3, 3, 3))

for(b in 1:3) {
  for(n in 1:3) {
    for(r in 1:1000) {
      biasAdj[b, 1, n, r] <- UnadjustedATE[b, 1, n, r] - TrueATE[1]
      biasAdj[b, 2, n, r] <- UnadjustedATE[b, 2, n, r] - TrueATE[2]
      biasAdj[b, 3, n, r] <- UnadjustedATE[b, 3, n, r] - TrueATE[3]
      biasAdjSE[b, 1, n, r] <- UnadjustedATESE[b, 1, n, r] - sd(UnadjustedATE[b, 1, n, ])
      biasAdjSE[b, 2, n, r] <- UnadjustedATESE[b, 2, n, r] - sd(UnadjustedATE[b, 2, n, ])
      biasAdjSE[b, 3, n, r] <- UnadjustedATESE[b, 3, n, r] - sd(UnadjustedATE[b, 3, n, ])
      
      biasLanza[b, 1, n, r] <- LanzaATE[b, 1, n, r] - TrueATE[1]
      biasLanza[b, 2, n, r] <- LanzaATE[b, 2, n, r] - TrueATE[2]
      biasLanza[b, 3, n, r] <- LanzaATE[b, 3, n, r] - TrueATE[3]
      biasLanzaSE[b, 1, n, r] <- LanzaATESE[b, 1, n, r] - sd(LanzaATE[b, 1, n, ])
      biasLanzaSE[b, 2, n, r] <- LanzaATESE[b, 2, n, r] - sd(LanzaATE[b, 2, n, ])
      biasLanzaSE[b, 3, n, r] <- LanzaATESE[b, 3, n, r] - sd(LanzaATE[b, 3, n, ])
      
      biasStep3[b, 1, n, r] <- Step3ATE[b, 1, n, r] - TrueATE[1]
      biasStep3[b, 2, n, r] <- Step3ATE[b, 2, n, r] - TrueATE[2]
      biasStep3[b, 3, n, r] <- Step3ATE[b, 3, n, r] - TrueATE[3]
      biasStep3SE[b, 1, n, r] <- Step3ATESE[b, 1, n, r] - sd(Step3ATE[b, 1, n, ])
      biasStep3SE[b, 2, n, r] <- Step3ATESE[b, 2, n, r] - sd(Step3ATE[b, 2, n, ])
      biasStep3SE[b, 3, n, r] <- Step3ATESE[b, 3, n, r] - sd(Step3ATE[b, 3, n, ])
    }
    AdjSD[b, 1, n] <- sd(UnadjustedATE[b, 1, n, ])
    AdjSD[b, 2, n] <- sd(UnadjustedATE[b, 2, n, ])
    AdjSD[b, 3, n] <- sd(UnadjustedATE[b, 3, n, ])
    
    LanzaSD[b, 1, n] <- sd(LanzaATE[b, 1, n, ])
    LanzaSD[b, 2, n] <- sd(LanzaATE[b, 2, n, ])
    LanzaSD[b, 3, n] <- sd(LanzaATE[b, 3, n, ])
    
    Step3SD[b, 1, n] <- sd(Step3ATE[b, 1, n, ])
    Step3SD[b, 2, n] <- sd(Step3ATE[b, 2, n, ])
    Step3SD[b, 3, n] <- sd(Step3ATE[b, 3, n, ])
  }
}


meanbiasAdj <- array(data = NA, dim = c(3, 3, 3))
meanbiasLanza <- array(data = NA, dim = c(3, 3, 3))
meanbiasStep3 <- array(data = NA, dim = c(3, 3, 3))

meanbiasAdjSE <- array(data = NA, dim = c(3, 3, 3))
meanbiasLanzaSE <- array(data = NA, dim = c(3, 3, 3))
meanbiasStep3SE <- array(data = NA, dim = c(3, 3, 3))

for(b in 1:3) {
  for(g in 1:3) {
    for(n in 1:3) {
      meanbiasAdj[b, g, n] <- mean(biasAdj[b, g, n, ])
      meanbiasLanza[b, g, n] <- mean(biasLanza[b, g, n, ])
      meanbiasStep3[b, g, n] <- mean(biasStep3[b, g, n, ])
      
      meanbiasAdjSE[b, g, n] <- mean(biasAdjSE[b, g, n, ])
      meanbiasLanzaSE[b, g, n] <- mean(biasLanzaSE[b, g, n, ])
      meanbiasStep3SE[b, g, n] <- mean(biasStep3SE[b, g, n, ])
    }
  }
}

method <- rep(c("one-step", "adjusted", "three-step"), each = 27)
matrix(meanbiasAdj, 27, 1)
betas <- rep(rep(c(1, 2, 3), each = 9), 3)
gammas <- rep(rep(c(1, 2, 3), each = 3), 9)
Nsample <- rep(c(500, 1000, 2500), 27)
bias <- rep(NA, 27)
SE <- rep(NA, 27)
SD <- rep(NA, 27)  

Myfinalresults1000 <- data.frame(method, betas, gammas, Nsample, bias, SE, SD)

Myfinalresults1000$bias[1:3] <- meanbiasLanza[1, 1, ]
Myfinalresults1000$bias[4:6] <- meanbiasLanza[1, 2, ]
Myfinalresults1000$bias[7:9] <- meanbiasLanza[1, 3, ]
Myfinalresults1000$bias[10:12] <- meanbiasLanza[2, 1, ]
Myfinalresults1000$bias[13:15] <- meanbiasLanza[2, 2, ]
Myfinalresults1000$bias[16:18] <- meanbiasLanza[2, 3, ]
Myfinalresults1000$bias[19:21] <- meanbiasLanza[3, 1, ]
Myfinalresults1000$bias[22:24] <- meanbiasLanza[3, 2, ]
Myfinalresults1000$bias[25:27] <- meanbiasLanza[3, 3, ]

Myfinalresults1000$bias[28:30] <- meanbiasAdj[1, 1, ]
Myfinalresults1000$bias[31:33] <- meanbiasAdj[1, 2, ]
Myfinalresults1000$bias[34:36] <- meanbiasAdj[1, 3, ]
Myfinalresults1000$bias[37:39] <- meanbiasAdj[2, 1, ]
Myfinalresults1000$bias[40:42] <- meanbiasAdj[2, 2, ]
Myfinalresults1000$bias[43:45] <- meanbiasAdj[2, 3, ]
Myfinalresults1000$bias[46:48] <- meanbiasAdj[3, 1, ]
Myfinalresults1000$bias[49:51] <- meanbiasAdj[3, 2, ]
Myfinalresults1000$bias[52:54] <- meanbiasAdj[3, 3, ]

Myfinalresults1000$bias[55:57] <- meanbiasStep3[1, 1, ]
Myfinalresults1000$bias[58:60] <- meanbiasStep3[1, 2, ]
Myfinalresults1000$bias[61:63] <- meanbiasStep3[1, 3, ]
Myfinalresults1000$bias[64:66] <- meanbiasStep3[2, 1, ]
Myfinalresults1000$bias[67:69] <- meanbiasStep3[2, 2, ]
Myfinalresults1000$bias[70:72] <- meanbiasStep3[2, 3, ]
Myfinalresults1000$bias[73:75] <- meanbiasStep3[3, 1, ]
Myfinalresults1000$bias[76:78] <- meanbiasStep3[3, 2, ]
Myfinalresults1000$bias[79:81] <- meanbiasStep3[3, 3, ]



Myfinalresults1000$SE[1:3] <- meanbiasLanzaSE[1, 1, ]
Myfinalresults1000$SE[4:6] <- meanbiasLanzaSE[1, 2, ]
Myfinalresults1000$SE[7:9] <- meanbiasLanzaSE[1, 3, ]
Myfinalresults1000$SE[10:12] <- meanbiasLanzaSE[2, 1, ]
Myfinalresults1000$SE[13:15] <- meanbiasLanzaSE[2, 2, ]
Myfinalresults1000$SE[16:18] <- meanbiasLanzaSE[2, 3, ]
Myfinalresults1000$SE[19:21] <- meanbiasLanzaSE[3, 1, ]
Myfinalresults1000$SE[22:24] <- meanbiasLanzaSE[3, 2, ]
Myfinalresults1000$SE[25:27] <- meanbiasLanzaSE[3, 3, ]

Myfinalresults1000$SE[28:30] <- meanbiasAdjSE[1, 1, ]
Myfinalresults1000$SE[31:33] <- meanbiasAdjSE[1, 2, ]
Myfinalresults1000$SE[34:36] <- meanbiasAdjSE[1, 3, ]
Myfinalresults1000$SE[37:39] <- meanbiasAdjSE[2, 1, ]
Myfinalresults1000$SE[40:42] <- meanbiasAdjSE[2, 2, ]
Myfinalresults1000$SE[43:45] <- meanbiasAdjSE[2, 3, ]
Myfinalresults1000$SE[46:48] <- meanbiasAdjSE[3, 1, ]
Myfinalresults1000$SE[49:51] <- meanbiasAdjSE[3, 2, ]
Myfinalresults1000$SE[52:54] <- meanbiasAdjSE[3, 3, ]

Myfinalresults1000$SE[55:57] <- meanbiasStep3SE[1, 1, ]
Myfinalresults1000$SE[58:60] <- meanbiasStep3SE[1, 2, ]
Myfinalresults1000$SE[61:63] <- meanbiasStep3SE[1, 3, ]
Myfinalresults1000$SE[64:66] <- meanbiasStep3SE[2, 1, ]
Myfinalresults1000$SE[67:69] <- meanbiasStep3SE[2, 2, ]
Myfinalresults1000$SE[70:72] <- meanbiasStep3SE[2, 3, ]
Myfinalresults1000$SE[73:75] <- meanbiasStep3SE[3, 1, ]
Myfinalresults1000$SE[76:78] <- meanbiasStep3SE[3, 2, ]
Myfinalresults1000$SE[79:81] <- meanbiasStep3SE[3, 3, ]



Myfinalresults1000$SD[1:3] <- LanzaSD[1, 1, ]
Myfinalresults1000$SD[4:6] <- LanzaSD[1, 2, ]
Myfinalresults1000$SD[7:9] <- LanzaSD[1, 3, ]
Myfinalresults1000$SD[10:12] <- LanzaSD[2, 1, ]
Myfinalresults1000$SD[13:15] <- LanzaSD[2, 2, ]
Myfinalresults1000$SD[16:18] <- LanzaSD[2, 3, ]
Myfinalresults1000$SD[19:21] <- LanzaSD[3, 1, ]
Myfinalresults1000$SD[22:24] <- LanzaSD[3, 2, ]
Myfinalresults1000$SD[25:27] <- LanzaSD[3, 3, ]

Myfinalresults1000$SD[28:30] <- AdjSD[1, 1, ]
Myfinalresults1000$SD[31:33] <- AdjSD[1, 2, ]
Myfinalresults1000$SD[34:36] <- AdjSD[1, 3, ]
Myfinalresults1000$SD[37:39] <- AdjSD[2, 1, ]
Myfinalresults1000$SD[40:42] <- AdjSD[2, 2, ]
Myfinalresults1000$SD[43:45] <- AdjSD[2, 3, ]
Myfinalresults1000$SD[46:48] <- AdjSD[3, 1, ]
Myfinalresults1000$SD[49:51] <- AdjSD[3, 2, ]
Myfinalresults1000$SD[52:54] <- AdjSD[3, 3, ]

Myfinalresults1000$SD[55:57] <- Step3SD[1, 1, ]
Myfinalresults1000$SD[58:60] <- Step3SD[1, 2, ]
Myfinalresults1000$SD[61:63] <- Step3SD[1, 3, ]
Myfinalresults1000$SD[64:66] <- Step3SD[2, 1, ]
Myfinalresults1000$SD[67:69] <- Step3SD[2, 2, ]
Myfinalresults1000$SD[70:72] <- Step3SD[2, 3, ]
Myfinalresults1000$SD[73:75] <- Step3SD[3, 1, ]
Myfinalresults1000$SD[76:78] <- Step3SD[3, 2, ]
Myfinalresults1000$SD[79:81] <- Step3SD[3, 3, ]




Nlabels <- c('500' = "N = 500", '1000' = "N = 1000", '2500' = "N = 2500")
blabels <- c('1' = "confounding = 1", '2' = "confounding = 2", '3' = "confounding = 3")

# Fig 2. Bias of the ATE
png("bias1000.png", width = 3100, height = 2100, res = 300)
ggplot(data = Myfinalresults1000, aes(x = gammas, y = bias, group = method, color = method)) +
  geom_point(size = 2) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  scale_y_continuous(breaks = c(-0.02, -0.01, 0, 0.01, 0.02), limits = c(-0.02, 0.02)) +
  scale_color_jco() +
  xlab("effect size") +
  ylab("bias of the average treatment effect") + 
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Nsample ~ betas, labeller = labeller(Nsample = as_labeller(Nlabels),
                                                  betas = as_labeller(blabels)))
dev.off()


# Fig 3. Standard deviation of the ATE
png("SD1000.png", width = 2900, height = 2100, res = 300)
ggplot(data = Myfinalresults1000, aes(x = gammas, y = SD, group = method, color = method)) +
  geom_point(size = 2) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1), limits = c(0, 0.12)) +
  scale_color_jco() +
  xlab("effect size") +
  ylab("SD of the average treatment effect") + 
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Nsample ~ betas, labeller = labeller(Nsample = as_labeller(Nlabels),
                                                  betas = as_labeller(blabels)))
dev.off()


# Fig 4. Bias of the standard errors
png("SEbias1000.png", width = 3100, height = 2100, res = 300)
ggplot(data = Myfinalresults1000, aes(x = gammas, y = SE, group = method, color = method)) +
  geom_point(size = 2) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  scale_y_continuous(breaks = c(-0.02, -0.01, 0, 0.01, 0.02), limits = c(-0.02, 0.02)) +
  scale_color_jco() +
  xlab("effect size") +
  ylab("bias of the standard errors") + 
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Nsample ~ betas, labeller = labeller(Nsample = as_labeller(Nlabels),
                                                  betas = as_labeller(blabels)))
dev.off()





##############################################################################################
# prostate cancer application
profiles.data <- read.spss("pr11a_EN_1.0.sav", to.data.frame = T)


# data preparation
profiles.data$ql <- profiles.data$pr11a09s01
profiles.data$pf <- profiles.data$pr11a09s02
profiles.data$rf <- profiles.data$pr11a09s03
profiles.data$ef <- profiles.data$pr11a09s04
profiles.data$cf <- profiles.data$pr11a09s05
profiles.data$sf <- profiles.data$pr11a09s06
profiles.data$fa <- profiles.data$pr11a09s07
profiles.data$nv <- profiles.data$pr11a09s08
profiles.data$pa <- profiles.data$pr11a09s09
profiles.data$dy <- profiles.data$pr11a09s10
profiles.data$sl <- profiles.data$pr11a09s11
profiles.data$ap <- profiles.data$pr11a09s12
profiles.data$co <- profiles.data$pr11a09s13
profiles.data$di <- profiles.data$pr11a09s14
profiles.data$fi <- profiles.data$pr11a09s15

profiles.data$fa <- 100 - profiles.data$fa
profiles.data$nv <- 100 - profiles.data$nv
profiles.data$pa <- 100 - profiles.data$pa
profiles.data$dy <- 100 - profiles.data$dy
profiles.data$sl <- 100 - profiles.data$sl
profiles.data$ap <- 100 - profiles.data$ap
profiles.data$co <- 100 - profiles.data$co
profiles.data$di <- 100 - profiles.data$di
profiles.data$fi <- 100 - profiles.data$fi

profiles.data$ID <- profiles.data$pr11a01pat_id
profiles.data$response <- profiles.data$pr11a01response
profiles.data$stage <- profiles.data$pr11a01stage
profiles.data$grade <- profiles.data$pr11a01grade
profiles.data$gleason <- profiles.data$pr11a01gleason
profiles.data$age <- profiles.data$pr11a01ageques
profiles.data$gender <- profiles.data$pr11a01gend
profiles.data$BMI <- profiles.data$pr11a01BMI
profiles.data$diag <- profiles.data$pr11a01yrsdiag
profiles.data$as <- profiles.data$pr11a01ww
profiles.data$smoke <- profiles.data$pr11a03q01
profiles.data$alc <- profiles.data$pr11a03q06

# treat: 1=treatment, 0=active surveillance
profiles.data$treat <- 0
profiles.data$treat[is.na(profiles.data$as)] <- 1

prostate.data <- profiles.data %>%
  filter(response == "respondent") %>%
  filter(stage == "Stage I" | stage == "Stage II") %>%
  filter(diag == "2-4 years" | diag == "4-6 years")

prostate.data <- prostate.data %>%
  select(ID, stage, gleason, grade, age, BMI, smoke, alc, diag, treat, ql, pf, rf, ef, cf, sf, fa, nv, pa, dy, sl, ap, co, di, fi)

prostate.data$age2 <- as.numeric(prostate.data$age)
prostate.data$age2[is.na(prostate.data$age2)] <- 9

prostate.data$gleason2 <- as.numeric(prostate.data$gleason)
prostate.data$gleason2[is.na(prostate.data$gleason2)] <- 9

prostate.data <- prostate.data %>%
  filter(!is.na(stage)) %>%
  filter(!is.na(gleason2)) %>%
  filter(!is.na(age2))

prostate.data$gleason3 <- car::recode(prostate.data$gleason2, "1 = 'gleason 2-6'; 2 = 'gleason 7-10'; 3 = 'gleason 7-10'; 9 = 'MV'")
prostate.data$gleason2 <- car::recode(prostate.data$gleason2, "1 = 'gleason 2-6'; 2 = 'gleason 7'; 3 = 'gleason 8-10'; 9 = 'MV'")

prostate.data$age2 <- case_when(
  prostate.data$age2 == 1 ~ "<=60 years",
  prostate.data$age2 == 2 ~ "<=60 years",
  prostate.data$age2 == 3~ ">60 years and <=65 years",
  prostate.data$age2 == 4 ~ ">65 years and <=70 years",
  prostate.data$age2 == 5 ~ ">70 years and <=75 years",
  prostate.data$age2 == 6 ~ ">75 years and <=80 years",
  prostate.data$age2 == 7 ~ ">80 years",
  prostate.data$age2 == 8 ~ ">80 years",
  prostate.data$age2 == 9 ~ "MV"
)

# estimate propensity score model
MyGLM <- glm(treat ~ stage + gleason2 + age2, family = "binomial", data = prostate.data)
PS_est <- MyGLM$fitted.values

prostate.data$PS <- PS_est

# generate inverse propensity weights
prostate.data$IPW <- ifelse(prostate.data$treat == 1, 1/prostate.data$PS, 1/(1-prostate.data$PS))


# create new dummy variables to assess balancing
prostate.data.w <- prostate.data

prostate.data.w$stage1 <- ifelse(prostate.data.w$stage == "Stage I", 1, 0)
prostate.data.w$stage2 <- ifelse(prostate.data.w$stage == "Stage II", 1, 0)

prostate.data.w$gleason2_6 <- ifelse(prostate.data.w$gleason == "gleason 2-6", 1, 0)
prostate.data.w$gleason7 <- ifelse(prostate.data.w$gleason == "gleason 7", 1, 0)
prostate.data.w$gleason8 <- ifelse(prostate.data.w$gleason == "gleason 8-10", 1, 0)
prostate.data.w$gleason7_10 <- ifelse(prostate.data.w$gleason3 == "gleason 7-10", 1, 0)
prostate.data.w$gleasonMV <- ifelse(prostate.data.w$gleason2 == "MV", 1, 0)

prostate.data.w$age60 <- ifelse(prostate.data.w$age2 == "<=60 years", 1, 0)
prostate.data.w$age65 <- ifelse(prostate.data.w$age2 == ">60 years and <=65 years", 1, 0)
prostate.data.w$age70 <- ifelse(prostate.data.w$age2 == ">65 years and <=70 years", 1, 0)
prostate.data.w$age75 <- ifelse(prostate.data.w$age2 == ">70 years and <=75 years", 1, 0)
prostate.data.w$age80 <- ifelse(prostate.data.w$age2 == ">75 years and <=80 years", 1, 0)
prostate.data.w$ageplus80 <- ifelse(prostate.data.w$age2 == ">80 years", 1, 0)
prostate.data.w$ageMV <- ifelse(prostate.data.w$age2 == "MV", 1, 0)

xvars <- c("stage1", "stage2", "gleason2_6", "gleason7", "gleason8", "gleasonMV", "age60", "age65", "age70", "age75", "age80", "ageplus80", "ageMV")


# trim most extreme weights
prostate.data.w.trim99 <- prostate.data.w %>%
  filter(IPW > quantile(prostate.data.w$IPW, probs = c(.01))) %>%
  filter(IPW < quantile(prostate.data.w$IPW, probs = c(.99)))

prostate.data.w.trim95 <- prostate.data.w %>%
  filter(IPW > quantile(prostate.data.w$IPW, probs = c(.05))) %>%
  filter(IPW < quantile(prostate.data.w$IPW, probs = c(.95)))  

# truncate most extreme weights
prostate.data.w.trunc99 <- prostate.data.w
prostate.data.w.trunc95 <- prostate.data.w

prostate.data.w.trunc99$IPW[prostate.data.w$IPW < quantile(prostate.data.w$IPW, probs = c(.01))] <- quantile(prostate.data.w$IPW, probs = c(.01))
prostate.data.w.trunc99$IPW[prostate.data.w$IPW > quantile(prostate.data.w$IPW, probs = c(.99))] <- quantile(prostate.data.w$IPW, probs = c(.99))

prostate.data.w.trunc95$IPW[prostate.data.w$IPW < quantile(prostate.data.w$IPW, probs = c(.05))] <- quantile(prostate.data.w$IPW, probs = c(.05))
prostate.data.w.trunc95$IPW[prostate.data.w$IPW > quantile(prostate.data.w$IPW, probs = c(.95))] <- quantile(prostate.data.w$IPW, probs = c(.95))

# generate stabilized weights
prop.table(table(prostate.data$treat))
prostate.data$IPWstab <- ifelse(prostate.data$treat == 1, 0.7802419/prostate.data$PS, 0.2197581/(1-prostate.data$PS))
prostate.data.w$IPWstab <- prostate.data$IPWstab

prostate.data.w.stab <- prostate.data.w


# check distribution of weights
plot(sort(prostate.data$IPW))
table(sort(prostate.data$IPW))

# check distribution of stabilized weights
table(sort(prostate.data$IPWstab))
plot(sort(prostate.data$IPWstab))


# create "Table One's" for different weigthing methods
table1 <- CreateTableOne(vars = xvars, strata = "treat", data = prostate.data.w, test = T)
print(table1, smd = T)

prostate.data.full <- svydesign(ids = ~1, data = prostate.data.w, weights = prostate.data.w$IPW)
table1.full <- svyCreateTableOne(vars = xvars, strata = "treat", data = prostate.data.full, test = T)
print(table1.full, smd = T)

prostate.data.trim95 <- svydesign(ids = ~1, data = prostate.data.w.trim95, weights = prostate.data.w.trim95$IPW)
table1.trim95 <- svyCreateTableOne(vars = xvars, strata = "treat", data = prostate.data.trim95, test = T)
print(table1.trim95, smd = T)

prostate.data.trim99 <- svydesign(ids = ~1, data = prostate.data.w.trim99, weights = prostate.data.w.trim99$IPW)
table1.trim99 <- svyCreateTableOne(vars = xvars, strata = "treat", data = prostate.data.trim99, test = T)
print(table1.trim99, smd = T)

prostate.data.trunc95 <- svydesign(ids = ~1, data = prostate.data.w.trunc95, weights = prostate.data.w.trunc95$IPW)
table1.trunc95 <- svyCreateTableOne(vars = xvars, strata = "treat", data = prostate.data.trunc95, test = T)
print(table1.trunc95, smd = T)

prostate.data.trunc99 <- svydesign(ids = ~1, data = prostate.data.w.trunc99, weights = prostate.data.w.trunc99$IPW)
table1.trunc99 <- svyCreateTableOne(vars = xvars, strata = "treat", data = prostate.data.trunc99, test = T)
print(table1.trunc99, smd = T)

prostate.data.stab <- svydesign(ids = ~1, data = prostate.data.w.stab, weights = prostate.data.w.stab$IPWstab)
table1.stab <- svyCreateTableOne(vars = xvars, strata = "treat", data = prostate.data.stab, test = T)
print(table1.stab, smd = T)


# export data set for application
export(prostate.data, "ProstateData.sav")


# check overlap of propensity scores
prostate.data$PSt <- ifelse(prostate.data$treat == 1, prostate.data$PS, 0)
prostate.data$PSnt <- ifelse(prostate.data$treat == 0, prostate.data$PS, 0)

prostate.data$PSt[prostate.data$PSt ==0] <- NA
prostate.data$PSnt[prostate.data$PSnt ==0] <- NA

png("PSoverlap.png", width = 3100, height = 2100, res = 300)
ggplot(prostate.data, aes(x=PS) ) +
  geom_histogram( aes(x = PSt, y = ..density..), fill="#0073C2FF" ) +
  geom_label( aes(x=.2, y=3, label="treatment"), color="#0073C2FF") + 
  geom_histogram( aes(x = PSnt, y = -..density..), fill= "#EFC000FF") +
  geom_label( aes(x=.2, y=-4, label="active surveillance"), color="#EFC000FF") +
  xlab("propensity score") +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15))
dev.off()


# generate latent Gold Syntax to run the different methods on the prostate cancer data
Step1Syntax <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile 'ProstateData.sav'
                                                    model
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  includeall;
                                                    output      
                                                    parameters=effect  betaopts=wl standarderrors profile probmeans=posterior
                                                    loadings bivariateresiduals estimatedvalues=model reorderclasses;
                                                    outfile  'prostate3class.sav'
                                                    classification=posterior      keep ID, stage, gleason, age, diag, treat, PS, IPW, gleason2, age2;
                                                    variables
                                                    dependent ql, pf, rf, ef, cf, sf, fa, nv, pa, dy, sl, ap, co, di,
                                                    fi;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1;
                                                    ql <- 1 + Cluster;
                                                    pf <- 1 + Cluster;
                                                    rf <- 1 + Cluster;
                                                    ef <- 1 + Cluster;
                                                    cf <- 1 + Cluster;
                                                    sf <- 1 + Cluster;
                                                    fa <- 1 + Cluster;
                                                    nv <- 1 + Cluster;
                                                    pa <- 1 + Cluster;
                                                    dy <- 1 + Cluster;
                                                    sl <- 1 + Cluster;
                                                    ap <- 1 + Cluster;
                                                    co <- 1 + Cluster;
                                                    di <- 1 + Cluster;
                                                    fi <- 1 + Cluster;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}



Step3Syntax <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile 'prostate3class.sav'
                                                    model
                                                    options
                                                    maxthreads=all;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    step3 proportional ml;  
                                                    output      
                                                    parameters=first betaopts=wl standarderrors=robust profile estimatedvalues=model append='prostateoutput.csv';
                                                    variables
                                                    samplingweight IPW ipw;
                                                    independent treat;
                                                    latent Cluster nominal posterior = ( Cluster#1 Cluster#2 Cluster#3 ) ;
                                                    equations
                                                    Cluster <- 1 + treat;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}



LanzaSyntax <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile 'ProstateData.sav'
                                                    model
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors profile bivariateresiduals reorderclasses estimatedvalues=model append='prostateLanza.csv';
                                                    variables
                                                    samplingweight IPW;
                                                    independent treat;
                                                    dependent ql, pf, rf, ef, cf, sf, fa, nv, pa, dy, sl, ap, co, di, fi;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + treat;
                                                    ql <- 1 + Cluster;
                                                    pf <- 1 + Cluster;
                                                    rf <- 1 + Cluster;
                                                    ef <- 1 + Cluster;
                                                    cf <- 1 + Cluster;
                                                    sf <- 1 + Cluster;
                                                    fa <- 1 + Cluster;
                                                    nv <- 1 + Cluster;
                                                    pa <- 1 + Cluster;
                                                    dy <- 1 + Cluster;
                                                    sl <- 1 + Cluster;
                                                    ap <- 1 + Cluster;
                                                    co <- 1 + Cluster;
                                                    di <- 1 + Cluster;
                                                    fi <- 1 + Cluster;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


NoWeightSyntax <- function(syntaxName){
  
  newSyntaxToBe <- utils::capture.output(cat(paste0("
                                                    //LG6.0//
                                                    version = 6.0
                                                    infile 'ProstateData.sav'
                                                    model
                                                    options
                                                    maxthreads=8;
                                                    algorithm 
                                                    tolerance=1e-008 emtolerance=0.01 emiterations=250 nriterations=50 ;
                                                    startvalues
                                                    seed=0 sets=16 tolerance=1e-005 iterations=50;
                                                    bayes
                                                    categorical=1 variances=1 latent=1 poisson=1;
                                                    montecarlo
                                                    seed=0 sets=0 replicates=500 tolerance=1e-008;
                                                    quadrature  nodes=10;
                                                    missing  excludeall;
                                                    output      
                                                    parameters=first  betaopts=wl standarderrors profile bivariateresiduals reorderclasses estimatedvalues=model append='prostateNoWeight.csv';
                                                    variables
                                                    independent treat, stage, gleason2, age2;
                                                    dependent ql, pf, rf, ef, cf, sf, fa, nv, pa, dy, sl, ap, co, di, fi;
                                                    latent
                                                    Cluster nominal 3;
                                                    equations
                                                    Cluster <- 1 + treat + stage + gleason2 + age2;
                                                    ql <- 1 + Cluster;
                                                    pf <- 1 + Cluster;
                                                    rf <- 1 + Cluster;
                                                    ef <- 1 + Cluster;
                                                    cf <- 1 + Cluster;
                                                    sf <- 1 + Cluster;
                                                    fa <- 1 + Cluster;
                                                    nv <- 1 + Cluster;
                                                    pa <- 1 + Cluster;
                                                    dy <- 1 + Cluster;
                                                    sl <- 1 + Cluster;
                                                    ap <- 1 + Cluster;
                                                    co <- 1 + Cluster;
                                                    di <- 1 + Cluster;
                                                    fi <- 1 + Cluster;
                                                    end model")))
  utils::write.table(newSyntaxToBe, paste0(syntaxName, ".lgs"), row.names = FALSE, quote = FALSE, col.names = FALSE)
}


# run different methods on prostate cancer data
Step1Syntax(syntaxName = "Step1Syntax")
Step3Syntax(syntaxName = "Step3Syntax")
LanzaSyntax(syntaxName = "LanzaSyntax")
NoWeightSyntax(syntaxName = "NoWeightSyntax")

shell(paste0(LG, " ", "Step1Syntax.lgs", " ", "/b"))
shell(paste0(LG, " ", "Step3Syntax.lgs", " ", "/b"))
shell(paste0(LG, " ", "LanzaSyntax.lgs", " ", "/b"))
shell(paste0(LG, " ", "NoWeightSyntax.lgs", " ", "/b"))


# Estimated values output from latent Gold was combined in this csv
tidyresults <- read.table("final results tidy data.csv", sep = ",", header = T)

tidyresults$class <- as.factor(tidyresults$class)
tidyresults$group <- as.factor(tidyresults$group)
tidyresults$method <- factor(tidyresults$method, levels = c('three-step', 'one-step', 'adjusted'))


# Fig 8. Main results of the real data application
png("resultsapllication.png", width = 3100, height = 2100, res = 300)
ggplot(tidyresults, aes(x = class, y = score, group = group, color = class, shape = group,
                        ymin = CI_low, ymax = CI_high)) +
  geom_point(position = position_dodge(width = .8), size = 1.5) +
  geom_pointrange(position = position_dodge(width = .8), size = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_jco() +
  ylab("class membership porbability") +
  facet_grid(. ~ method) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12))
dev.off()
