#Adapted from Klaus Oberauer; Colin Stoneking; Dominik Wabersich; Hsuan-Yu Lin(2017): "Hierarchical Bayesian measurement models for continuous reproduction of visual features from working memory"


#clear variables and graphs
rm(list=ls(all=TRUE))
graphics.off()

#load libraries and helpers
library(rjags)
library(circular)
library(coda)
load.module("vonmises") #source: Oberauer et al., 2017 (osf.io/wjg7y/) 
load.module("dic")
root<-"C:/Users/ytabi/Documents/Hierarchical Bayesian measurement models/Bayes Model Data/" #where your data is stored

for(EXPT in 1:3) {
for(effects in 1:7){
postpred <- 0
diags <- 0 #0 = none, 1 = Gelman-Rubin Rhat
parameters <- c("seqpos1_effect_prec_mean", "seqpos2_effect_prec_mean", "seqpos3_effect_prec_mean", "seqpos1_effect_prec_sd", "seqpos2_effect_prec_sd", "seqpos3_effect_prec_sd", "seqpos1_effect_ptar_mean", "seqpos2_effect_ptar_mean", "seqpos3_effect_ptar_mean", "seqpos1_effect_ptar_sd", "seqpos2_effect_ptar_sd", "seqpos3_effect_ptar_sd", "factor1_effect_prec_mean", "factor1_effect_prec_sd", "factor1_effect_ptar_mean", "factor1_effect_ptar_sd", "factor1_seqpos1_effect_prec_mean", "factor1_seqpos2_effect_prec_mean", "factor1_seqpos3_effect_prec_mean", "factor1_seqpos1_effect_prec_sd", "factor1_seqpos2_effect_prec_sd", "factor1_seqpos3_effect_prec_sd", "factor1_seqpos1_effect_ptar_mean", "factor1_seqpos2_effect_ptar_mean", "factor1_seqpos3_effect_ptar_mean", "factor1_seqpos1_effect_ptar_sd", "factor1_seqpos2_effect_ptar_sd", "factor1_seqpos3_effect_ptar_sd","pitem_i","prec_i","factor1_eff_prec","factor1_eff_ptar","seqpos1_eff_prec","seqpos1_eff_ptar","seqpos2_eff_prec","seqpos2_eff_ptar","seqpos3_eff_prec","seqpos3_eff_ptar","factor1_seqpos1_eff_prec","factor1_seqpos1_eff_ptar","factor1_seqpos2_eff_prec","factor1_seqpos2_eff_ptar","factor1_seqpos3_eff_prec","factor1_seqpos3_eff_ptar""factor2_effect_prec_mean","factor2_effect_prec_sd","factor2_effect_ptar_mean","factor2_effect_ptar_sd","factor1_factor2_effect_prec_mean","factor1_factor2_effect_prec_sd","factor1_factor2_effect_ptar_mean","factor1_factor2_effect_ptar_sd","factor2_eff_prec","factor2_eff_ptar","factor1_factor2_eff_prec","factor1_factor2_eff_ptar") #model parameter names
parameters<-c(parameters, "loglik")


nadapt <- 5000        #number of adaptations
niterations <- 10000  #number of iterations
nchains <- 4          #number of chains

setwd(root)                                   #use scripts that can be found in this folder

#load the .dat files we have created with RecallCues_Convert.m
m           <- read.table(paste0(root,"EXPT ", toString(EXPT), "/m.dat"))       #memoryset
x           <- read.table(paste0(root,"EXPT ", toString(EXPT), "/x.dat"))       #response
id          <- read.table(paste0(root,"EXPT ", toString(EXPT), "/id.dat"))      #id number
setsize     <- read.table(paste0(root,"EXPT ", toString(EXPT), "/setsize.dat")) #setsize
condition   <- read.table(paste0(root,"EXPT ", toString(EXPT), "/condition.dat")) #condition number
condnum     <- nrow(unique(condition))

n.params.prepred <- length(parameters)                      #number of parameters included
if (postpred == 1) parameters <- c(parameters, "postpred")

  ntrials <- length(t(id))  
  data = list (
    "m" = cbind(m, rep(0,ntrials)),                               #add column of zeros at the end for the Eb component
    "x" = as.numeric(unlist(x)),                                  #response
    "N" = ntrials,                                                #number of trials
    "nsubj" = length(t(unique(id))),                              #number of subjects
    "id" = as.numeric(unlist(id)),                                #id
    "kappafactor" = c(rep(1,length(m)), 0),                       #kappa factor
    "setsize" = as.numeric(unlist(setsize)),                      #setsize

#dummy variables

    "factor1"= NULL,
    "factor2"= NULL,                                              

    "seqpos1"= NULL,
    "seqpos2"= NULL,
    "seqpos3"= NULL,
    "seqpos4"= NULL
)

if (EXPT==1 | EXPT==2){

  if (effects==1){ #factor1 and factor2 have an influence on precision and target response
    modelname   <-"mixture3_2x2.bug" #use the 2x2-design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3)*1)) #probe type
    data$factor2<-as.numeric(unlist((condition == 3 | condition == 4)*1)) #setsize

  } else if (effects==2){ #only factor1 has an influence on precision and target response
    modelname   <-"mixture3_2.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3)*1)) 

  } else if (effects==3){ #only factor2 has an influence on precision and target response
    modelname   <-"mixture3_2.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 3 | condition == 4)*1))

  } else if (effects==4){ #only factor2 has an influence on precision and target response
    modelname   <-"mixture3.bug" #use the one factor, 2 levels design

  } else if (effects==5){ #only factor2 has an influence on precision and target response
    modelname   <-"mixture3_2x2_nointeraction.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3)*1)) #probe type
    data$factor2<-as.numeric(unlist((condition == 3 | condition == 4)*1)) #setsize

  } else if (effects==6){ #only factor2 has an influence on precision and target response
    modelname   <-"mixture3_2_precisiononly.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3)*1)) #probe type

  } else if (effects==7){ #only factor2 has an influence on precision and target response
    modelname   <-"mixture3_2_ptargetonly.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3)*1)) #probe type


}
 



} else {

  if (effects==1){ #factor1 and sequence have an influence precision and target response
    modelname   <-"mixture3_2x4.bug" #use the 2x4 design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3 | condition == 5 | condition == 7)*1)); #probe type
    data$seqpos1<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 3 | condition == 4)*1))  #sequence positions
    data$seqpos2<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 5 | condition == 6)*1))
    data$seqpos3<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 7 | condition == 8)*1))

  } else if (effects==3){ #only sequence has an influence precision and target response
    modelname   <-"mixture3_4.bug" #use the 2x4 design
    data$seqpos1<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 3 | condition == 4)*1))
    data$seqpos2<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 5 | condition == 6)*1))
    data$seqpos3<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 7 | condition == 8)*1))
 
 } else if (effects==2){ #only factor1 has an influence precision and target response
    modelname   <-"mixture3_2.bug" #use the 2x4 design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3 | condition == 5 | condition == 7)*1));


  } else if (effects==4){ #only factor1 has an influence precision and target response
    modelname   <-"mixture3.bug" #use the one factor, 2 levels design
 

  } else if (effects==5){ #only factor1 has an influence precision and target response
    modelname   <-"mixture3_2x4_nointeraction.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3 | condition == 5 | condition == 7)*1)); #probe type
    data$seqpos1<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 3 | condition == 4)*1))  #sequence positions
    data$seqpos2<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 5 | condition == 6)*1))
    data$seqpos3<-as.numeric(unlist((condition == 1 | condition == 2 | condition == 7 | condition == 8)*1))
 

  } else if (effects==6){ #only factor1 has an influence precision and target response
    modelname   <-"mixture3_2_precisiononly.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3 | condition == 5 | condition == 7)*1));
 

  } else if (effects==7){ #only factor1 has an influence precision and target response
    modelname   <-"mixture3_2_ptargetonly.bug" #use the one factor, 2 levels design
    data$factor1<-as.numeric(unlist((condition == 1 | condition == 3 | condition == 5 | condition == 7)*1));
  }

}

  try (
    if (1 == 1) {
      # fitting with rjags
      jagsmodel<-jags.model(modelname, data=data, n.chains=nchains, n.adapt=nadapt)
      codasamples <- coda.samples(jagsmodel, variable.names=parameters,
                                  n.iter=niterations, thin=1)  # turn object of class "jags" into "mcmc.list"
      

     
      samples <- as.matrix(codasamples)

      # calculate DIC
      Devchain <- samples[, "deviance"]
      meanDev <- mean(Devchain)
      pD <- 0.5*var(Devchain)
      DIC <- meanDev + pD
      dicfile <- paste0(root, "Recall_Cues_Interfere.DIC")
      dicthere = file.access(dicfile, mode=0)
      if (dicthere == 0) load(dicfile) else DICX <- matrix(0, 3, 7)
      colnames(DICX) <- c("both factors(*)", "probe type only", "additional factor only", "no factor", "both factors", "presiononly", "ptargetonly")
      DICX[EXPT,effects] <- DIC
      save(DICX, file=dicfile)
      
      # calculate WAIC
      lppd <- 0  #log pointwise predicted density (Equation 5 in Gelman et al, 2014)
      pwaic <- 0 # WAIC penalty
      for (i in 1:ntrials) {
        LLidx <- paste("loglik[", as.character(i), "]", sep="")
        lppd <- lppd + log(mean(exp(samples[, LLidx]))) # adds the mean lppd for each data point
        pwaic <- pwaic + var(samples[, LLidx]) # adds the penalty pWAIC2 (Equation 12 in Gelman et al., 2014)
      }
      WAIC <- -2*(lppd-pwaic) # compute WAIC, multiply with -2 to translate into deviance scale
      waicfile <- paste0(root, "Recall_Cues_Interfere.WAIC")
      waicthere = file.access(waicfile, mode=0)
      if (waicthere == 0) load(waicfile) else WAICX <- matrix(0, 3, 7)
      colnames(WAICX) <- c("both factors(*)", "probe type only", "additional factor only","no factor", "both factors", "presiononly", "ptargetonly")
      WAICX[EXPT,effects] <- WAIC
      save(WAICX, file=waicfile)
      

      save(file=paste0(root, "EXPT_", toString(EXPT), "_effects_",  toString(effects),".RData"))
    } 
  )  #try
  } #effects
} # Experiment

