# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")
#sessionInfo()

mainDir <- "P:/FDR dependence/Model Selection/logoffset_in_simulation_logoffset_in_model_negative_beta_positive_beta_50_50/result"

#mainDir <- "/home/ntyet/research/modelselection2" # linux
# mainDir <- "P:/research/modelselection2" # Windows
#dir.create(file.path(mainDir, "OldData_Offset_QLfit"), showWarnings = FALSE)
newDir <- file.path(mainDir,"EdgeR_glmQLFTest")

pbeta1 <- "pbeta_1"
dir.create(file.path(newDir, pbeta1), showWarnings = FALSE)
oldDir.pbeta1 <- file.path(mainDir,pbeta1)

newDir.pbeta1 <- file.path(newDir, pbeta1)

library("xtable")
library("plyr")
library("edgeR")
library("qvalue")
library("locfit")
#library("QuasiSeq")
#detach("QuasiSeq")
library("maps")
library(clinfun)
library("fields")
library("fdrtool")
source('P:/qvalue_1.34.0/qvalue/R/qvalue_m0.R')
source("http://www.public.iastate.edu/~dnett/microarray/multtest.txt")
source("P:/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/hybrid-poisson-test-vc.R")
source("P:/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/Hyprid_poisson_test_vc_modified0.R")
source("P:/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/Hyprid_poisson_test_vc_modified1.R")
source("P:/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/fdrtool_1.2.10/fdrtool/R/ecdf.pval.R")
source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\QL.fit2.R")
source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\NBDev.R")
source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\PoisDev.R")
source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\QL.results.R")

# # RFI output data
# dat <- read.table("U:/R/RA/Data/RFI_uniq_comb_count_corrected.txt")
# meta.data <- read.csv("U:/R/RA/Data/RIN values RNAseq Obj3a Martine.csv")
# meta.data2 <- read.table("U:/R/RA/Data/Meta.data2.txt")
# 
# name <- paste("X",meta.data$Sample.Name, sep = "")
# dat2 <- dat[, name]
# #colnames(dat2)
# Lane.RNAseq <- as.factor(meta.data$Lane.RNAseq)
# Diet <- meta.data$diet
# Line <- meta.data$line
# #RFI.value <- meta.data2$RFI.value
# RFI.value <- meta.data2$RFI.value
# RIN.before.GD <- meta.data$RIN.before.GD
# RIN.after.GD  <- meta.data$RIN.after.GD
# Conc.after.GD <- meta.data$Conc.after.GD.ng.ul
# date.GD <- meta.data$date.GD
# date.RNA.extraction <- meta.data$date.RNA.extraction
# counts <- as.matrix(dat2[rowSums(dat2>0)>1&
#                            rowMeans(dat2)>1,])
# dim(counts)
# 
# summary(log(counts+1))
# mean.1 <- apply(log(counts[,Line ==1]+1), 1, mean)
# mean.2 <- apply(log(counts[,Line ==2]+1), 1, mean)
# summary(mean.1 - mean.2)
# 
# wd <- "P:\\FDR dependence\\Model Selection\\logoffset_in_simulation_logoffset_in_model_negative_beta_positive_beta_50_50\\result\\pbeta_1"
# 
# sim.dat1 <- load(file = paste(wd, "\\p.beta_1i.beta_1.5e.beta_2S_2L_0.1U_0.5m_29.RData", sep = ""))
# str(sim.dat1)
# counts <- sim1$counts
# summary(log(counts+1))
# mean.1 <- apply(log(counts[,1:20]+1), 1, mean)
# mean.2 <- apply(log(counts[,21:40]+1), 1, mean)
# summary(mean.1 - mean.2)
# str(sim1)


# source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\QL.fit2.R")
# source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\NBDev.R")
# source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\PoisDev.R")
# source("P:\\stevescode\\QuasiSeq_1.0-2\\QuasiSeq\\R\\QL.results.R")

# Load RFI count.mean and NBdisp from model0RFI.line

# 
# load(file = "U:/R/RA/Data/Additional Plot/Model0.line.rfi.RData")
# load(file = "U:/R/RA/Data/Additional Plot/Model0.result.line.rfi.RData")
# #counts <- as.matrix(dat2[rowSums(dat2>0)>1&
# #                           rowMeans(dat2)>1,])
# #################################################
# mean.count <- fit.line.rfi$mn.cnt
# #counts
# nb.disp <- fit.line.rfi$NB.disp



I <- 2
J <- 1000
K <- 20
DE <- round(J*.2)
EE <- J - DE

#######################################################
#######################################################
#SAT.LIKE2 to calculate likelihood 

# function to compute the constant term in pdf of NegBin with 
# paramters mean \mu, observation y and dispersion parameter 1/disp
log.gamma <- function(counts, disp){
  log.g <- NULL
  n <- length(counts)
  for(i in 1:n){
    log.g[i] <- sum(log(0:(counts[i]-1)+disp)) - sum(log(1:counts[i]) )
  }
  return(log.g)  
}

SAT.LIKE2<-function(data,disp){
  means<-data
  like<-disp*log(disp/(disp+means))
  like[data!=0]<-like[data!=0]+data[data!=0]*log(means[data!=0]/(disp+means[data!=0]))+
    log.gamma(data[data!=0],disp )
  sum(like)
}

## Function to calculate AIC of the QL.fit model 

# AIC.QL <- function(counts,QL.fit.object){
#   n <- dim(counts)[2]
#   m <- dim(counts)[1]
#   disp <- 1/QL.fit.object$NB.disp
#   den.df <- QL.fit.object$den.df
#   phi.hat.dev <- QL.fit.object$phi.hat.dev
#   p <- n - den.df
#   dev <- phi.hat.dev*den.df
#   L0 <- NULL
#   for (i in 1:m){
#     L0[i] <- SAT.LIKE2(counts[i,],disp[i])
#   }
#   
#   return(dev-2*L0+2*p)
# }


##############


#######AIC of the glmfit 
#glmfitobject <- fit

AIC.glmfit <- function(counts,glmfitobject){
  n <- dim(counts)[2]
  m <- dim(counts)[1]
  disp <- 1/glmfitobject$dispersion
  p <- dim(glmfitobject$design)[2]
  dev <- glmfitobject$deviance
  L0 <- NULL
  for (i in 1:m){
    L0[i] <- SAT.LIKE2(counts[i,],disp[i])
  }
  
  return(dev-2*L0+2*p)
}
#i <- 1
###

# the following function is the same as function pval.hist
# in the paper of Pounds et al. 2012, EBT, with the estimation of 
# density obtained from the paper by Korbinian Strimmer 2008
# 
pval.hist.grenander <- function(p.value){
  grenander.out <- grenander(ecdf(p.value))
  p.brks <- c(0, grenander.out$x.knots)
  b.edf <- c(0, grenander.out$F.knots)
  p.diffs <- diff(p.brks)
  h.cdf <- approx(p.brks, b.edf, xout = p.value)$y  # get the histogram-based CDF estimates from each p-value
  p.hist <- exp(log(diff(b.edf))-log(diff(p.brks))) # get the hight for each histogram bar
  pi0.hat <- min(p.hist)                            # get the pi0 estimate from histogram bar
  h.ebp <- approx(p.brks, pi0.hat/c(p.hist, p.hist[length(p.hist)]), xout = p.value)$y # get the conservative EBP interpolation 
  h.fdr <- exp(log(pi0.hat) + log(p.value) - log(h.cdf))                                     # Get the histogram based FDR estimate
  h.ebp[p.value==0] <- 0
  h.fdr[p.value==0] <- 0
  return(list( p.value = p.value,          # input p-value,
               h.cdf = h.cdf,              # the histogram Grenander based cdf estimate
               h.fdr = h.fdr,              # the histogram Grenander based FDR
               h.ebp = h.ebp,              # the histogram Grenander based EBP
               p.brks = p.brks,            # the p-value break-points of the histogram
               p.hist = p.hist,            # the heights of each histogram bar
               edf.brks = b.edf,           # the breaks points in the EDF of the histogram estimator
               pi0.hat = pi0.hat))         # the histogram Grenander based estimate of the proportion of tests with a true null hypothesis
}

# function to claculate auc based on test value

auc <- function(test){
  marker <- 1 - test
  D <- c(rep(1,DE), rep(0,EE))
  out <- roc.curve(marker = marker, status = D, method = "empirical")
  return(as.numeric(strsplit(capture.output(print(out)), " ")[[1]][8]))
}



## FDP 
# function to calculate FDP at threshold qvalue fixed controled by Storey FDR, estimate.m0 by Nettleton

fdp <- function(pvalue, qvalue.threshold){
  n <- length(qvalue.threshold)
  out <- NULL
  q.value <- qvalue.m0(pvalue)
  for(i in 1:n){
    if (sum(q.value<=qvalue.threshold[i])==0){
      out[i] <- 0
    } else{
      out[i] <- sum((q.value<=qvalue.threshold[i])&(1:J>DE))/sum(q.value<=qvalue.threshold[i])
    }
  }
  return(out)
}



##Sim counts data
# 
# sim.counts <- function(p.beta, i.beta, e.beta, S, L, U){
#   omega <- NULL
#   tau <- matrix(0, nrow = I, ncol = J) # treatment value
#   x <- matrix(0, nrow = I, ncol = K) # covariate value
#   beta <- matrix(0, nrow = I, ncol = J) # coefficients of covariate
#   phi <- matrix(0, nrow = I, ncol = K) # offset value
#   b <- NULL
#   b1 <- NULL
#   b2 <- NULL
#   
#   # generate phi and covariate x
#   
#   for(i in 1:I){
#     for(k in 1:K){
#       x[i,k] <- rnorm(1,0,1)
#       #phi[i,k] <- rnorm(1,0,.125^2)
#       phi[i,k] <- rnorm(1,0,0)
#     }
#   }
#   
#   # genetating beta which is the coefficients of covariates
#   
#   beta.ind <- NULL
#   
#   for(j in 1:J){
#     beta.ind[j] <- sample(c(0,1,-1), size = 1,
#                           prob = c(1-p.beta, p.beta/2,p.beta/2))
#     beta[1,j] <- beta[2,j] <- beta.ind[j]*runif(1,i.beta, e.beta)
#   }
#   
#   
#   mu <- array(0, dim = c(I,J,K))
#   y <- array(0, dim = c(I,J,K))
#   mean.sim <- NULL
#   
#   for(j in 1:J){
#     repeat{
#       if (1<=j &j <= DE){
#         j.p <- sample(1:length(mean.count), 1 )
#         mean.sim[j] <- mean.count[j.p]
#         b1[j] <- 1/rgamma(1, shape = S*log(mean.count[j.p])^(1/8), scale = 1)
#         b2[j] <- runif(1, L,U)
#         b[j] <- b1[j] +b2[j]
#         trt.ind <- rbinom(1,1,.5)
#         if (trt.ind ==1){
#           tau[1,j] <- log(mean.count[j.p]/sqrt(b[j]))
#           tau[2,j] <- log(mean.count[j.p]*sqrt(b[j]))
#         } else{
#           tau[1,j] <- log(mean.count[j.p]/sqrt(b[j]))
#           tau[2,j] <- log(mean.count[j.p]*sqrt(b[j]))
#         }
#         omega[j] <- nb.disp[j.p]
#       }
#       if ((DE +1) <= j){
#         j.p <- sample(1:length(mean.count), 1 )
#         mean.sim[j] <- mean.count[j.p]
#         tau[1,j] <- log(mean.count[j.p])
#         tau[2,j] <- log(mean.count[j.p])
#         omega[j] <- nb.disp[j.p] 
#       }
#       for(i in 1:I){
#         for(k in 1:K){
#           mu[i,j,k] <- exp(phi[i,k]+tau[i,j]+beta[i,j]*x[i,k])
#           y[i,j,k]  <- rnbinom(n=1, size=1/omega[j], mu=mu[i,j,k])
#         }
#       }
#       if (mean(y[,j,])>1& sum(y[,j,]>0)>1) break
#     }
#   }
#   
#   counts <- cbind(y[1,,],y[2,,])
#   return(list(counts= counts, x = x, 
#               beta.ind = beta.ind, 
#               trt.ind = trt.ind, 
#               beta = beta, mean.sim = mean.sim, 
#               tau = tau, 
#               mu =cbind(mu[1,,],mu[2,,]),
#               b = b))
# }
# 
# 
# 
# 
# # Initial parameter for check functions
# # S <- 1.25
# # L <- 1
# # U <- 1.5
# # 
# # p.beta <- 1
# # i.beta <- .1
# # e.beta <- .5
# 

#############################################3
# ?glmLRT
# attributes(fit)
# str(fit$table)
# 
# mean.count <- log(apply(sim1$counts, 1, mean))
# mean.count1 <- log(apply(sim1$counts[,1:K],1,mean)+1)
# mean.count2 <- log(apply(sim1$counts[,K+1:K],1,mean)+1)
# log.fold <- -mean.count1 + mean.count2
# plot(mean.count,log.fold)
# #sim1$trt.ind
# length(mean.count)
# length(log.fold)

###############################################
sim.edger <- function(sim1){
  counts <- sim1$counts
  beta.ind <- sim1$beta.ind
  x <- sim1$x
  design <- model.matrix(~rep(1:2, each = K))
  counts <- DGEList(counts, group = design[,2])
  counts <- calcNormFactors(counts)
  counts <- estimateGLMCommonDisp(counts, design, verbose = TRUE)
  counts <- estimateGLMTrendedDisp(counts, design)
  counts <- estimateGLMTagwiseDisp(counts, design)
  fit <- glmFit(counts, design)
  lrt <- glmLRT(fit, coef = 2)
  
  #   counts <- estimateDisp(counts, design, robust = TRUE)
  #   fit <- glmQLFTest(counts, design, robust=TRUE, plot = F, coef = 2)
  #   
  #   pvalue.trt.nocov <- lrt$table$PValue
  #   hist(pvalue.trt.nocov, nclass = 100)
  #   #AICAIC.glmfit(sim1$counts, fit)
  #   
  
  #############################################################
  
  # Design list with covariate in full model 
  design <- model.matrix(~rep(1:2, each = K)+x)
  counts <- sim1$counts
  counts <- DGEList(counts, group = design[,2])
  counts <- calcNormFactors(counts)
  #   counts <- estimateDisp(counts, design, robust = T)
  #   fit2 <- glmQLFTest(counts, design, robust=TRUE, plot = F, coef = 2)
  #   fit3 <- glmQLFTest(counts, design, robust=TRUE, plot = F, coef = 3)
  #   fit2$dispersion
  #   fit3$dispersion
  
  #   str(fit2$table$PValue)
  counts <- estimateGLMCommonDisp(counts, design, verbose = TRUE)
  counts <- estimateGLMTrendedDisp(counts, design)
  counts <- estimateGLMTagwiseDisp(counts, design)
  fit2 <- glmFit(counts, design)
  lrt2 <- glmLRT(fit2, coef = 2) # trt
  lrt3 <- glmLRT(fit2, coef = 3) # cov
  #   
  
  #dev.off()
  #pvalue.trt.cov <- fit2$table$PValue
  pvalue.trt.cov <- lrt2$table$PValue
  #hist(pvalue.trt.cov, nclass = 100)
  
  #pvalue.cov <- fit3$table$PValue
  pvalue.cov <- lrt3$table$PValue
  #hist(pvalue.cov, nclass = 100)
  pvalue.trt.nocov <- lrt$table$PValue
  #pvalue.trt.cov <- lrt2$table$PValue
  
  ebp.func <- pval.hist.modified0(pvalue.cov)
  g.ebp.func <- pval.hist.grenander(pvalue.cov)
  #fit$dispersion
  aic.nocov <- AIC.glmfit(sim1$counts, fit)
  
  aic.cov <- AIC.glmfit(sim1$counts, fit2)
  
  pvalue.trt.ebp <- NULL
  pvalue.trt.g.ebp <- NULL
  pvalue.trt.aic <- NULL
  pvalue.trt.true <- NULL
  #aic.nocov
  
  for (j in 1:J){
    if (ebp.func$h.ebp[j] < .5) pvalue.trt.ebp[j] <- pvalue.trt.cov[j]
    if (ebp.func$h.ebp[j] >= .5) pvalue.trt.ebp[j] <- pvalue.trt.nocov[j]
    
    if (g.ebp.func$h.ebp[j] < .5) pvalue.trt.g.ebp[j] <- pvalue.trt.cov[j]
    if (g.ebp.func$h.ebp[j] >= .5) pvalue.trt.g.ebp[j] <- pvalue.trt.nocov[j]
    
    if (aic.nocov[j] <= aic.cov[j]) pvalue.trt.aic[j] <- pvalue.trt.nocov[j]
    if (aic.nocov[j] >aic.cov[j]) pvalue.trt.aic[j] <- pvalue.trt.cov[j]
    
    if (beta.ind[j]==0) pvalue.trt.true[j] <- pvalue.trt.nocov[j]
    if (beta.ind[j]!=0) pvalue.trt.true[j] <- pvalue.trt.cov[j]
  }
  
  #########################AAA methods with AHE and Grenander Estimate###
  
  ebp.func <- pval.hist.modified0(pvalue.cov)
  ebp.trt.nocov <- pval.hist.modified0(pvalue.trt.nocov)
  ebp.trt.cov <- pval.hist.modified0(pvalue.trt.cov)
  aebp.trt <- ebp.trt.nocov$h.ebp*ebp.func$h.ebp +ebp.trt.cov$h.ebp*(1-ebp.func$h.ebp)
  
  
  g.ebp.func <- pval.hist.grenander(pvalue.cov)
  g.ebp.trt.nocov <- pval.hist.grenander(pvalue.trt.nocov)
  g.ebp.trt.cov <- pval.hist.grenander(pvalue.trt.cov)
  g.aebp.trt <- g.ebp.trt.nocov$h.ebp*g.ebp.func$h.ebp +g.ebp.trt.cov$h.ebp*(1-g.ebp.func$h.ebp)
  
  
  # ROC curves for 7 methods
  auc.nocov <- auc(pvalue.trt.nocov)
  auc.cov <- auc(pvalue.trt.cov)
  auc.ebp <-  auc(pvalue.trt.ebp)
  auc.g.ebp<-  auc(pvalue.trt.g.ebp)
  auc.aic <- auc(pvalue.trt.aic)
  auc.aebp <- auc(aebp.trt)
  auc.g.aebp <- auc(g.aebp.trt)
  auc.true <- auc(pvalue.trt.true)
  
  #qvalue, FDP for first five methods using Storey FDR with estimate.m0 by Nettleton et. al. 2006
  
  
  t <- c(0.01, 0.05, 0.10, .15)
  fdp.nocov <- fdp(pvalue.trt.nocov, t)
  fdp.cov <- fdp(pvalue.trt.cov, t)
  fdp.ebp <- fdp(pvalue.trt.ebp, t)
  fdp.g.ebp <- fdp(pvalue.trt.g.ebp, t)
  fdp.aic <- fdp(pvalue.trt.aic, t)
  fdp.true <- fdp(pvalue.trt.true, t)
  res <- list(auc.nocov = auc.nocov,
              auc.cov = auc.cov ,
              auc.ebp = auc.ebp ,
              auc.g.ebp = auc.g.ebp,
              auc.aic = auc.aic ,
              auc.aebp = auc.aebp,
              auc.g.aebp = auc.g.aebp,
              auc.true = auc.true,
              fdp.nocov = fdp.nocov,
              fdp.cov = fdp.cov ,
              fdp.ebp = fdp.ebp ,
              fdp.g.ebp = fdp.g.ebp,
              fdp.aic = fdp.aic ,
              fdp.true = fdp.true)          
  return(res)
}

##############################
# 
# sim.QLfit2 <- function(sim1){
#   counts <- sim1$counts
#   beta.ind <- sim1$beta.ind
#   x <- sim1$x
#   design.list <- vector("list",2)
#   design.list[[1]] <- rep(1:2, each = K)
#   design.list[[2]] <- rep(1, ncol(counts))
#   #size <- apply(counts, 2, quantile, .75)
#   size <- apply(counts, 2, sum)
#   fit <- QL.fit(counts, design.list, 
#                 #log.offset = log(size),
#                 Model = "NegBin",
#                 print.progress=FALSE)
#   #traceback()
#   result.fit <- QL.results(fit, Plot=FALSE)
#   
#   #############################################################
#   
#   # Design list with covariate in full model 
#   design.list <- vector("list", 3)
#   x1 <- as.factor(rep(1:2, each = K))
#   design.list[[1]] <- model.matrix(~ x1 + x)
#   design.list[[2]] <- model.matrix(~ x1) # test for covariate
#   design.list[[3]] <- model.matrix(~ x) # test for treatment effect
#   test.mat <- rbind(1:2, c(1,3))
#   row.names(test.mat) <- c("Covariate", "Treatment")
#   fit.2 <- QL.fit(counts, design.list, 
#                   #log.offset = log(size),
#                   test.mat,
#                   Model="NegBin", 
#                   print.progress=FALSE)
#   
#   result.fit2 <- QL.results(fit.2,Plot= FALSE)
#   
#   pvalue.cov <- result.fit2$P.values[[3]][,1]
#   
#   pvalue.trt.nocov <- result.fit$P.values[[3]][,1]
#   pvalue.trt.cov <- result.fit2$P.values[[3]][,2]
#   hist(pvalue.trt.nocov, nclass = 100)
#   hist(pvalue.trt.cov, nclass = 100)
#   ebp.func <- pval.hist.modified0(pvalue.cov)
#   g.ebp.func <- pval.hist.grenander(pvalue.cov)
#   
#   aic.nocov <- AIC.QL(counts, fit)
#   aic.cov <- AIC.QL(counts, fit.2)
#   
#   pvalue.trt.ebp <- NULL
#   pvalue.trt.g.ebp <- NULL
#   pvalue.trt.aic <- NULL
#   pvalue.trt.true <- NULL
#   for (j in 1:J){
#     if (ebp.func$h.ebp[j] < .5) pvalue.trt.ebp[j] <- pvalue.trt.cov[j]
#     if (ebp.func$h.ebp[j] >= .5) pvalue.trt.ebp[j] <- pvalue.trt.nocov[j]
#     
#     if (g.ebp.func$h.ebp[j] < .5) pvalue.trt.g.ebp[j] <- pvalue.trt.cov[j]
#     if (g.ebp.func$h.ebp[j] >= .5) pvalue.trt.g.ebp[j] <- pvalue.trt.nocov[j]
#     
#     if (aic.nocov[j] <= aic.cov[j]) pvalue.trt.aic[j] <- pvalue.trt.nocov[j]
#     if (aic.nocov[j] >aic.cov[j]) pvalue.trt.aic[j] <- pvalue.trt.cov[j]
#     
#     if (beta.ind[j]==0) pvalue.trt.true[j] <- pvalue.trt.nocov[j]
#     if (beta.ind[j]!=0) pvalue.trt.true[j] <- pvalue.trt.cov[j]
#   }
#   
#   #########################AAA methods with AHE and Grenander Estimate###
#   
#   ebp.func <- pval.hist.modified0(pvalue.cov)
#   ebp.trt.nocov <- pval.hist.modified0(pvalue.trt.nocov)
#   ebp.trt.cov <- pval.hist.modified0(pvalue.trt.cov)
#   aebp.trt <- ebp.trt.nocov$h.ebp*ebp.func$h.ebp +ebp.trt.cov$h.ebp*(1-ebp.func$h.ebp)
#   
#   
#   g.ebp.func <- pval.hist.grenander(pvalue.cov)
#   g.ebp.trt.nocov <- pval.hist.grenander(pvalue.trt.nocov)
#   g.ebp.trt.cov <- pval.hist.grenander(pvalue.trt.cov)
#   g.aebp.trt <- g.ebp.trt.nocov$h.ebp*g.ebp.func$h.ebp +g.ebp.trt.cov$h.ebp*(1-g.ebp.func$h.ebp)
#   
#   
#   # ROC curves for 7 methods
#   auc.nocov <- auc(pvalue.trt.nocov)
#   auc.cov <- auc(pvalue.trt.cov)
#   auc.ebp <-  auc(pvalue.trt.ebp)
#   auc.g.ebp<-  auc(pvalue.trt.g.ebp)
#   auc.aic <- auc(pvalue.trt.aic)
#   auc.aebp <- auc(aebp.trt)
#   auc.g.aebp <- auc(g.aebp.trt)
#   auc.true <- auc(pvalue.trt.true)
#   
#   #qvalue, FDP for first five methods using Storey FDR with estimate.m0 by Nettleton et. al. 2006
#   
#   
#   t <- c(0.01, 0.05, 0.10, .15)
#   fdp.nocov <- fdp(pvalue.trt.nocov, t)
#   fdp.cov <- fdp(pvalue.trt.cov, t)
#   fdp.ebp <- fdp(pvalue.trt.ebp, t)
#   fdp.g.ebp <- fdp(pvalue.trt.g.ebp, t)
#   fdp.aic <- fdp(pvalue.trt.aic, t)
#   fdp.true <- fdp(pvalue.trt.true, t)
#   res <- list(auc.nocov = auc.nocov,
#               auc.cov = auc.cov ,
#               auc.ebp = auc.ebp ,
#               auc.g.ebp = auc.g.ebp,
#               auc.aic = auc.aic ,
#               auc.aebp = auc.aebp,
#               auc.g.aebp = auc.g.aebp,
#               auc.true = auc.true,
#               fdp.nocov = fdp.nocov,
#               fdp.cov = fdp.cov ,
#               fdp.ebp = fdp.ebp ,
#               fdp.g.ebp = fdp.g.ebp,
#               fdp.aic = fdp.aic ,
#               fdp.true = fdp.true)          
#   return(res)
# }


p.beta <- c(0,  0.2, 0.5 , 0.8, 1)
i.beta <- c(0.1, 1.5)
e.beta <- c(.5,  2)
S <- c(1.25, 2)
L <- c(0.1 ,1.5)
U <- c(0.5, 2)

n.sim <- 100
auc.nocov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.cov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.g.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.aic <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.true <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.aebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))
auc.g.aebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim))

fdp.nocov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim, 4))
fdp.cov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim, 4))
fdp.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim, 4))
fdp.g.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim, 4))
fdp.aic <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim, 4))
fdp.true <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L), n.sim, 4))

res.auc.nocov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.cov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.g.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.aic <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.true <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.aebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
res.auc.g.aebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))

res.fdp.nocov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
res.fdp.cov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
res.fdp.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
res.fdp.g.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
res.fdp.aic <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
res.fdp.true <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))



sd.res.auc.nocov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.cov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.g.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.aic <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.true <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.aebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))
sd.res.auc.g.aebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L)))

sd.res.fdp.nocov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
sd.res.fdp.cov <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
sd.res.fdp.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
sd.res.fdp.g.ebp <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
sd.res.fdp.aic <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))
sd.res.fdp.true <- array(0, dim=c(length(p.beta), length(i.beta), length(S), length(L),4))


#

# set.seed(1)
# sim <- sim.counts(p.beta[5], 3, 4, S[1], 3, 4)
# sim <- sim.counts(p.beta[5], .1, .5, S[1], 3, 4)
# sim <- sim.counts(p.beta[5], 3, 4, S[1], .1, .5)
# counts <- sim$counts
# mu <- sim$mu
# tau <- sim$tau
# dim(mu)
# x <- sim$x
# dim(x)
# beta <- sim$beta
# mean.sim <- sim$mean.sim
# dev.off()
# 
# j <- i <- 1
# log(mean.sim[j])
# plot(c(x[1,],x[2,]), log(mu[j,]))
# plot(c(x[1,],x[2,]), log(counts[j,]+1))
# tau[,j]
# b <- sim$b
# log(b[j])
# log(mean.sim[j])
#     # counts[23,1:K]
#     # counts[23,K+1:K]
#     # 
#     # # i <- 7,13,28,31,34,36
#     # e <- 0
#     # #i <- 2
#     #dev.off()
#     i <- j <- 18
#     e <- 1
#     par(mfrow = c(2,2))
#     for (i in 1:4+ 4*e)
#     {
#       log(counts[i,1:K]+1)-log(counts[i,K+1:K]+1)
#       
#       plot(c(x[1,],x[2,]), log(counts[i,]+1),main =i)
#       points(x[1,], log(counts[i,1:K]+1), pch = 16, col = "red")
#       mod1 <- lm(log(counts[i,]+1)~ rep(c(1,2), each = K))
#       mod2 <- lm(log(counts[i,]+1)~ rep(c(1,2), each = K) + c(x[1,],x[2,]))
#       #str(anova(mod1))
#       anova(mod1)[[5]][1]
#       anova(mod2)[[5]][1]
#       }
#     

################
##############
# i <- 4; j <- 1;k <- 1; l <- 1; m <- 4

for(i in 5){
  for(j in 1:length(i.beta)){
    for(k in 1:length(S)){
      for(l in 1:length(U)){
        for(m in 1:n.sim){
          pathsave <- paste(oldDir.pbeta1, 
                            "/p.beta_",
                            p.beta[i], 
                            "i.beta_",
                            i.beta[j],
                            "e.beta_",
                            e.beta[j],
                            "S_",
                            S[k],
                            "L_",
                            L[l],
                            "U_",
                            U[l],
                            "m_",
                            m,
                            ".RData",sep = "")
          load(file = pathsave)
          out <- sim.edger(sim1)
          traceback()
          auc.nocov[i,j,k,l,m] <- out$auc.nocov
          auc.cov[i,j,k,l,m] <- out$auc.cov
          auc.ebp[i,j,k,l,m] <- out$auc.ebp
          auc.g.ebp[i,j,k,l,m] <- out$auc.g.ebp
          auc.aic[i,j,k,l,m] <- out$auc.aic
          auc.true[i,j,k,l,m] <- out$auc.true
          auc.aebp[i,j,k,l,m] <- out$auc.aebp
          auc.g.aebp[i,j,k,l,m] <- out$auc.g.aebp
          
          fdp.nocov[i,j,k,l,m, ] <- out$fdp.nocov
          fdp.cov[i,j,k,l,m, ] <- out$fdp.cov
          fdp.ebp[i,j,k,l,m, ] <- out$fdp.ebp
          fdp.g.ebp[i,j,k,l,m, ] <- out$fdp.g.ebp
          fdp.aic[i,j,k,l,m, ] <- out$fdp.aic
          fdp.true[i,j,k,l,m, ] <- out$fdp.true
          print(paste("i = ", i, ",j = ", j, ",k = ", k, ",l = ", l, ",m = ", m))
        }
        res.auc.nocov[i,j,k,l] <- mean(auc.nocov[i,j,k,l,])
        res.auc.cov[i,j,k,l] <- mean(auc.cov[i,j,k,l,])
        res.auc.ebp[i,j,k,l] <- mean(auc.ebp[i,j,k,l,])
        res.auc.g.ebp[i,j,k,l] <- mean(auc.g.ebp[i,j,k,l,])
        res.auc.aic[i,j,k,l] <- mean(auc.aic[i,j,k,l,])
        res.auc.true[i,j,k,l] <- mean(auc.true[i,j,k,l,])
        res.auc.aebp[i,j,k,l] <- mean(auc.aebp[i,j,k,l,])
        res.auc.g.aebp[i,j,k,l] <- mean(auc.g.aebp[i,j,k,l,])
        
        sd.res.auc.nocov[i,j,k,l] <- sd(auc.nocov[i,j,k,l,])
        sd.res.auc.cov[i,j,k,l] <- sd(auc.cov[i,j,k,l,])
        sd.res.auc.ebp[i,j,k,l] <- sd(auc.ebp[i,j,k,l,])
        sd.res.auc.g.ebp[i,j,k,l] <- sd(auc.g.ebp[i,j,k,l,])
        sd.res.auc.aic[i,j,k,l] <- sd(auc.aic[i,j,k,l,])
        sd.res.auc.true[i,j,k,l] <- sd(auc.true[i,j,k,l,])
        sd.res.auc.aebp[i,j,k,l] <- sd(auc.aebp[i,j,k,l,])
        sd.res.auc.g.aebp[i,j,k,l] <- sd(auc.g.aebp[i,j,k,l,])
        
        
        
        res.fdp.nocov[i,j,k,l,] <- apply(fdp.nocov[i,j,k,l,,], 2, mean)
        res.fdp.cov[i,j,k,l,] <- apply(fdp.cov[i,j,k,l,,], 2, mean)
        res.fdp.ebp[i,j,k,l,] <- apply(fdp.ebp[i,j,k,l,,], 2, mean)
        res.fdp.g.ebp[i,j,k,l,] <- apply(fdp.g.ebp[i,j,k,l,,], 2, mean)
        res.fdp.aic[i,j,k,l,] <- apply(fdp.aic[i,j,k,l,,], 2, mean)
        res.fdp.true[i,j,k,l,] <- apply(fdp.true[i,j,k,l,,], 2, mean)
        
        
        sd.res.fdp.nocov[i,j,k,l,] <- apply(fdp.nocov[i,j,k,l,,], 2, sd)
        sd.res.fdp.cov[i,j,k,l,] <- apply(fdp.cov[i,j,k,l,,], 2, sd)
        sd.res.fdp.ebp[i,j,k,l,] <- apply(fdp.ebp[i,j,k,l,,], 2, sd)
        sd.res.fdp.g.ebp[i,j,k,l,] <- apply(fdp.g.ebp[i,j,k,l,,], 2, sd)
        sd.res.fdp.aic[i,j,k,l,] <- apply(fdp.aic[i,j,k,l,,], 2, sd)
        sd.res.fdp.true[i,j,k,l,] <- apply(fdp.true[i,j,k,l,,], 2, sd)
        
      }
    }
  }
  
}
traceback()
save(res.auc.nocov, file = paste(newDir.pbeta1, "/res.auc.nocov.RData", sep = ""))
save(res.auc.cov, file = paste(newDir.pbeta1, "/res.auc.cov.RData", sep = ""))
save(res.auc.ebp, file = paste(newDir.pbeta1, "/res.auc.ebp.RData", sep = ""))
save(res.auc.g.ebp, file = paste(newDir.pbeta1, "/res.auc.g.ebp.RData", sep = ""))
save(res.auc.aic, file = paste(newDir.pbeta1, "/res.auc.aic.RData", sep = ""))
save(res.auc.true, file = paste(newDir.pbeta1, "/res.auc.true.RData", sep = ""))
save(res.auc.aebp, file = paste(newDir.pbeta1, "/res.auc.aebp.RData", sep = ""))
save(res.auc.g.aebp, file = paste(newDir.pbeta1, "/res.auc.g.aebp.RData", sep = ""))

save(res.fdp.nocov, file = paste(newDir.pbeta1, "/res.fdp.nocov.RData", sep = ""))
save(res.fdp.cov, file = paste(newDir.pbeta1, "/res.fdp.cov.RData", sep = ""))
save(res.fdp.ebp, file = paste(newDir.pbeta1, "/res.fdp.ebp.RData", sep = ""))
save(res.fdp.g.ebp, file = paste(newDir.pbeta1, "/res.fdp.g.ebp.RData", sep = ""))
save(res.fdp.aic, file = paste(newDir.pbeta1, "/res.fdp.aic.RData", sep = ""))
save(res.fdp.true, file = paste(newDir.pbeta1, "/res.fdp.true.RData", sep = ""))

save(sd.res.auc.nocov, file = paste(newDir.pbeta1, "/sd.res.auc.nocov.RData", sep = ""))
save(sd.res.auc.cov, file = paste(newDir.pbeta1, "/sd.res.auc.cov.RData", sep = ""))
save(sd.res.auc.ebp, file = paste(newDir.pbeta1, "/sd.res.auc.ebp.RData", sep = ""))
save(sd.res.auc.g.ebp, file = paste(newDir.pbeta1, "/sd.res.auc.g.ebp.RData", sep = ""))
save(sd.res.auc.aic, file = paste(newDir.pbeta1, "/sd.res.auc.aic.RData", sep = ""))
save(sd.res.auc.true, file = paste(newDir.pbeta1, "/sd.res.auc.true.RData", sep = ""))
save(sd.res.auc.aebp, file = paste(newDir.pbeta1, "/sd.res.auc.aebp.RData", sep = ""))
save(sd.res.auc.g.aebp, file = paste(newDir.pbeta1, "/sd.res.auc.g.aebp.RData", sep = ""))



save(sd.res.fdp.nocov, file = paste(newDir.pbeta1, "/sd.res.fdp.nocov.RData", sep = ""))
save(sd.res.fdp.cov, file = paste(newDir.pbeta1, "/sd.res.fdp.cov.RData", sep = ""))
save(sd.res.fdp.ebp, file = paste(newDir.pbeta1, "/sd.res.fdp.ebp.RData", sep = ""))
save(sd.res.fdp.g.ebp, file = paste(newDir.pbeta1, "/sd.res.fdp.g.ebp.RData", sep = ""))
save(sd.res.fdp.aic, file = paste(newDir.pbeta1, "/sd.res.fdp.aic.RData", sep = ""))
save(sd.res.fdp.true, file = paste(newDir.pbeta1, "/sd.res.fdp.true.RData", sep = ""))


# auc of those methods: 
i <- 5


ind <- expand.grid(j = c(1,2), k = c(1,2), l = c(1,2))

auc.all <-function(x) {c(     paste(round(res.auc.nocov[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.nocov[i,ind[x,1],ind[x,2],ind[x,3]], 4),")", sep = ""),
                              
                              paste(round(res.auc.cov[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.cov[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep =""),
                              paste(round(res.auc.ebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.ebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep = ""),
                              paste(round(res.auc.g.ebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.g.ebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep=""),
                              paste(round(res.auc.aic[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.aic[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep=""),
                              paste(round(res.auc.true[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.true[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep = ""),
                              paste(round(res.auc.aebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.aebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep=""),
                              paste(round(res.auc.g.aebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),"(",
                                    round(sd.res.auc.g.aebp[i,ind[x,1],ind[x,2],ind[x,3]], 4),")",sep="")
)}



auc.out <- laply(1:dim(ind)[1], function(x)auc.all(x))

dimnames(auc.out) = 
  list(c("(1,1,1)", 
         "(2,1,1)",
         "(1,2,1)",
         "(2,2,1)",
         "(1,1,2)",
         "(2,1,2)",
         "(1,2,2)",
         "(2,2,2)")
       , 
       c("nocov", 
         "cov", 
         "ebp", 
         "g.ebp", 
         "aic", 
         "true", 
         "aepb", 
         "g.aebp"))

auc.out

xtable(auc.out)


# FDR = .01



fdp.all <-function(x) {c(     paste(round(res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],1], 4),"(",
                                    round(sd.res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],1], 4),")", sep = ""),
                              
                              paste(round(res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],1], 4),"(",
                                    round(sd.res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],1], 4),")",sep =""),
                              paste(round(res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],1], 4),"(",
                                    round(sd.res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],1], 4),")",sep = ""),
                              paste(round(res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],1], 4),"(",
                                    round(sd.res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],1], 4),")",sep=""),
                              paste(round(res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],1], 4),"(",
                                    round(sd.res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],1], 4),")",sep=""),
                              paste(round(res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],1], 4),"(",
                                    round(sd.res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],1], 4),")",sep = "")
)}


#library("plyr")
fdp.out <- laply(1:dim(ind)[1], function(x)fdp.all(x))
fdp.out
dimnames(fdp.out) = 
  list(c("(1,1,1)", 
         "(2,1,1)",
         "(1,2,1)",
         "(2,2,1)",
         "(1,1,2)",
         "(2,1,2)",
         "(1,2,2)",
         "(2,2,2)")
       , 
       c("nocov", 
         "cov", 
         "ebp", 
         "g.ebp", 
         "aic", 
         "true"))

fdp.out
#install.packages("xtable")
#library("xtable")
xtable(fdp.out)



#FDR = .05

fdp.all <-function(x) {c(     paste(round(res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],2], 4),"(",
                                    round(sd.res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],2], 4),")", sep = ""),
                              
                              paste(round(res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],2], 4),"(",
                                    round(sd.res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],2], 4),")",sep =""),
                              paste(round(res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],2], 4),"(",
                                    round(sd.res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],2], 4),")",sep = ""),
                              paste(round(res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],2], 4),"(",
                                    round(sd.res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],2], 4),")",sep=""),
                              paste(round(res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],2], 4),"(",
                                    round(sd.res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],2], 4),")",sep=""),
                              paste(round(res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],2], 4),"(",
                                    round(sd.res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],2], 4),")",sep = "")
)}

#fdp
#library("plyr")
fdp.out <- laply(1:dim(ind)[1], function(x)fdp.all(x))
fdp.out
dimnames(fdp.out) = 
  list(c("(1,1,1)", 
         "(2,1,1)",
         "(1,2,1)",
         "(2,2,1)",
         "(1,1,2)",
         "(2,1,2)",
         "(1,2,2)",
         "(2,2,2)")
       , 
       c("nocov", 
         "cov", 
         "ebp", 
         "g.ebp", 
         "aic", 
         "true"))

fdp.out
#install.packages("xtable")
library("xtable")
xtable(fdp.out)


####FDP = .1


fdp.all <-function(x) {c(     paste(round(res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],3], 4),"(",
                                    round(sd.res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],3], 4),")", sep = ""),
                              
                              paste(round(res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],3], 4),"(",
                                    round(sd.res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],3], 4),")",sep =""),
                              paste(round(res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],3], 4),"(",
                                    round(sd.res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],3], 4),")",sep = ""),
                              paste(round(res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],3], 4),"(",
                                    round(sd.res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],3], 4),")",sep=""),
                              paste(round(res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],3], 4),"(",
                                    round(sd.res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],3], 4),")",sep=""),
                              paste(round(res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],3], 4),"(",
                                    round(sd.res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],3], 4),")",sep = "")
)}


#library("plyr")
fdp.out <- laply(1:dim(ind)[1], function(x)fdp.all(x))
fdp.out
dimnames(fdp.out) = 
  list(c("(1,1,1)", 
         "(2,1,1)",
         "(1,2,1)",
         "(2,2,1)",
         "(1,1,2)",
         "(2,1,2)",
         "(1,2,2)",
         "(2,2,2)")
       , 
       c("nocov", 
         "cov", 
         "ebp", 
         "g.ebp", 
         "aic", 
         "true"))

fdp.out
#install.packages("xtable")
#library("xtable")
xtable(fdp.out)

# FDR = .15


fdp.all <-function(x) {c(     paste(round(res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],4], 4),"(",
                                    round(sd.res.fdp.nocov[i,ind[x,1],ind[x,2],ind[x,3],4], 4),")", sep = ""),
                              
                              paste(round(res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],4], 4),"(",
                                    round(sd.res.fdp.cov[i,ind[x,1],ind[x,2],ind[x,3],4], 4),")",sep =""),
                              paste(round(res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],4], 4),"(",
                                    round(sd.res.fdp.ebp[i,ind[x,1],ind[x,2],ind[x,3],4], 4),")",sep = ""),
                              paste(round(res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],4], 4),"(",
                                    round(sd.res.fdp.g.ebp[i,ind[x,1],ind[x,2],ind[x,3],4], 4),")",sep=""),
                              paste(round(res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],4], 4),"(",
                                    round(sd.res.fdp.aic[i,ind[x,1],ind[x,2],ind[x,3],4], 4),")",sep=""),
                              paste(round(res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],4], 4),"(",
                                    round(sd.res.fdp.true[i,ind[x,1],ind[x,2],ind[x,3],4], 4),")",sep = "")
)}


#library("plyr")
fdp.out <- laply(1:dim(ind)[1], function(x)fdp.all(x))
fdp.out
dimnames(fdp.out) = 
  list(c("(1,1,1)", 
         "(2,1,1)",
         "(1,2,1)",
         "(2,2,1)",
         "(1,1,2)",
         "(2,1,2)",
         "(1,2,2)",
         "(2,2,2)")
       , 
       c("nocov", 
         "cov", 
         "ebp", 
         "g.ebp", 
         "aic", 
         "true"))

fdp.out
#install.packages("xtable")
#library("xtable")
xtable(fdp.out)
