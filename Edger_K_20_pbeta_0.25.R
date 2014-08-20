library("edgeR");library("plyr");library("fdrtool");library("AUC"); library("maps") ;library("fields")
I <- 2; J <- 1000
K <- 20
DE <- round(J*.2)
EE <- J - DE
S <- 1.25
L <- 0.1
U <- 0.5
p.beta <- 0.25
i.beta <- c(0.1, 1)
e.beta <- c(0.5, 1.5)
n.sim <- 100
mainDir <- "/home/ntyet/research/nuisancecovariate" # linux server
mainDir1 <- paste("/home/ntyet/research/nuisancecovariate/edgeR/K_", K, sep = "")  # linux server

# mainDir <- "/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet/research/nuisancecovariate" # Linux laptop
# mainDir1 <- paste("/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet/research/nuisancecovariate/edgeR/K_", K,sep = "")  # linux server
# mainDir <- "P:/research/nuisancecovariate/edgeR" # Linux laptop
# mainDir1 <- paste("P:/research/nuisancecovariate/edgeR/K_", K,sep = "")  # linux server

dir.create(mainDir1, showWarnings = FALSE)
pbeta1 <- paste("pbeta_", p.beta, sep = "")
dir.create(file.path(mainDir1, pbeta1), showWarnings = FALSE)
sources <- "sources"
#dir.create(file.path(mainDir1, sources), showWarnings = FALSE)
dir.source <- file.path(mainDir, sources)
dir.pbeta1 <- file.path(mainDir1, pbeta1)

# for linux 
source(paste(dir.source, '/qvalue_1.34.0/qvalue/R/qvalue_m0.R', sep = ""))
source("http://www.public.iastate.edu/~dnett/microarray/multtest.txt")
source(paste(dir.source, "/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/hybrid-poisson-test-vc.R", sep = ""))

source(paste(dir.source, "/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/Hyprid_poisson_test_vc_modified0.R",sep =""))
source(paste(dir.source, "/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/Hyprid_poisson_test_vc_modified1.R",sep =""))


load(file = paste(dir.source, "/Model0.line.rfi.RData",sep = ""))
load(file = paste(dir.source,"/Model0.result.line.rfi.RData",sep = ""))
str(result.line.rfi)
# hist(result.line.rfi$P.values[[3]][,"RFI.value"])
# str(result.line.rfi$P.values[[3]])
#counts <- as.matrix(dat2[rowSums(dat2>0)>1&
#                           rowMeans(dat2)>1,])

mean.count <- fit.line.rfi$mn.cnt
#counts
nb.disp <- fit.line.rfi$NB.disp
#SAT_LIKE2 to calculate likelihood 

# function to compute the constant term in pdf of NegBin with 
# paramters mean \mu, observation y and dispersion parameter 1/disp
log_gamma <- function(counts, disp){
  log.g <- laply(1:length(counts), function(i) sum(log(0:(counts[i]-1)+disp)) - sum(log(1:counts[i]) ))
  return(log.g)  
}

SAT_LIKE2<-function(count,disp){
  means<-count
  like<-disp*log(disp/(disp+means))
  like[count!=0]<-like[count!=0]+count[count!=0]*log(means[count!=0]/(disp+means[count!=0]))+
    log_gamma(count[count!=0],disp )
  sum(like)
}

# AIC_edger
AIC_glmfit <- function(counts,glmfitobject){
  n <- dim(counts)[2]
  m <- dim(counts)[1]
  disp <- 1/glmfitobject$dispersion
  p <- dim(glmfitobject$design)[2]
  dev <- glmfitobject$deviance
  L0 <- laply(1:m, function(i)SAT_LIKE2(counts[i,],disp[i]))
  return(dev-2*L0+2*p)
}

# the following function is the same as function pval.hist
# in the paper of Pounds et al. 2012, EBT, with the estimation of 
# density obtained from the paper by Korbinian Strimmer 2008
# 


pval_hist_grenander <- function(p.value){
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


# auc

lab <- as.factor(c(rep(1,DE), rep(0,EE)))
auc_out <- function(test.vector){
  roc.out <- roc(1-test.vector, lab)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  auc1 <- auc(roc.out)
  auc2 <- auc(roc.out, min =roc.min)
  return(c(auc1,auc2))
}
## FDP 
# function to calculate FDP at threshold qvalue fixed controled by Storey FDR, estimate.m0 by Nettleton

fdp <- function(pvalue){
  q.value <- jabes.q(pvalue, B = 20)
    ifelse(sum(q.value<=0.05)==0, 0, 
           sum((q.value<=0.05)&(1:J>DE))/sum(q.value <= 0.05))
}
#ebp.temp <- c(0.1, 0.5, 0.4, .3, .02, 1)

#qvalue.threshold <- .05

fdr_ebp <- function(ebp.temp){
  ord <- order(ebp.temp)
  ebp.temp <- ebp.temp[ord]
  fdr.temp <- cumsum(ebp.temp)/1:length(ebp.temp)
  fdr.temp[ord] <- fdr.temp    
  ifelse(sum(fdr.temp<=0.05)==0, 0, 
         sum((fdr.temp<=0.05)&(1:J>DE))/sum(fdr.temp<=0.05))
}

### simcount 

sim_counts <- function(p.beta, i.beta, e.beta, S, L, U){
  omega <- NULL
  tau <- matrix(0, nrow = I, ncol = J) # treatment value
  x <- matrix(0, nrow = I, ncol = K) # covariate value
  beta <- matrix(0, nrow = I, ncol = J) # coefficients of covariate
  phi <- matrix(0, nrow = I, ncol = K) # offset value
  b <- NULL
  b1 <- NULL
  b2 <- NULL
  
  #  generate phi and covariate x
  for(k in 1:K){
    x[1,k] <-x[2,k] <-  rnorm(1,0,1)
    phi[1,k] <- phi[2,k] <- rnorm(1,0,0)
  }
  
  
  # genetating beta which is the coefficients of covariates
  
  beta.ind <- NULL
  
  for(j in 1:J){
    beta.ind[j] <- sample(c(0,1,-1), size = 1,
                          prob = c(1-p.beta, p.beta/2,p.beta/2))
    beta[1,j] <- beta[2,j] <- beta.ind[j]*runif(1,i.beta, e.beta)
  }
  
  
  mu <- array(0, dim = c(I,J,K))
  y <- array(0, dim = c(I,J,K))
  mean.sim <- NULL
  trt.pos <- NULL
  for(j in 1:J){
    repeat{
      if (1<=j &j <= DE){
        j.p <- sample(1:length(mean.count), 1 )
        mean.sim[j] <- mean.count[j.p]
        trt.pos[j] <- j.p
        b1[j] <- 1/rgamma(1, shape = S*log(mean.count[j.p])^(1/8), scale = 1)
        b2[j] <- runif(1, L,U)
        b[j] <- b1[j] +b2[j]
        trt.ind <- rbinom(1,1,.5)
        if (trt.ind ==1){
          tau[1,j] <- log(mean.count[j.p]/sqrt(b[j]))
          tau[2,j] <- log(mean.count[j.p]*sqrt(b[j]))
        } else{
          tau[2,j] <- log(mean.count[j.p]/sqrt(b[j]))
          tau[1,j] <- log(mean.count[j.p]*sqrt(b[j]))
        }
        omega[j] <- nb.disp[j.p]
      }
      if ((DE +1) <= j){
        j.p <- sample(1:length(mean.count), 1 )
        mean.sim[j] <- mean.count[j.p]
        tau[1,j] <- log(mean.count[j.p])
        tau[2,j] <- log(mean.count[j.p])
        omega[j] <- nb.disp[j.p] 
      }
      for(i in 1:I){
        for(k in 1:K){
          mu[i,j,k] <- exp(phi[i,k]+tau[i,j]+beta[i,j]*x[i,k])
          y[i,j,k]  <- rnbinom(n=1, size=1/omega[j], mu=mu[i,j,k])
        }
      }
      if (mean(y[,j,])>1& sum(y[,j,]>0)>1) break
    }
  }
  
  counts <- cbind(y[1,,],y[2,,])
  return(list(counts= counts, x = x, 
              beta.ind = beta.ind, 
              trt.pos = trt.pos, 
              beta = beta, mean.sim = mean.sim, 
              tau = tau, 
              mu =cbind(mu[1,,],mu[2,,]),
              b = b))
}


##edger_fit


edger_fit <- function(p.beta, i.beta, e.beta, S, L, U){
  sim.data <- sim_counts(p.beta, i.beta, e.beta, S, L, U) # i.beta <- 0.1, e.beta <- 0.5
  counts <- sim.data$counts
  beta.ind <- sim.data$beta.ind
  x <- sim.data$x
  design <- model.matrix(~rep(1:2, each = K))
  counts <- DGEList(counts, group = design[,2])
  counts <- calcNormFactors(counts)
  counts <- estimateGLMCommonDisp(counts, design, verbose = TRUE)
  counts <- estimateGLMTrendedDisp(counts, design)
  counts <- estimateGLMTagwiseDisp(counts, design)
  fit <- glmFit(counts, design)
  lrt <- glmLRT(fit, coef = 2)
  
  # Design list with covariate in full model 
  design <- model.matrix(~rep(1:2, each = K)+c(x[1,],x[2,]))
  counts <- sim.data$counts
  counts <- DGEList(counts, group = design[,2])
  counts <- calcNormFactors(counts)
  counts <- estimateGLMCommonDisp(counts, design, verbose = TRUE)
  counts <- estimateGLMTrendedDisp(counts, design)
  counts <- estimateGLMTagwiseDisp(counts, design)
  fit2 <- glmFit(counts, design)
  lrt2 <- glmLRT(fit2, coef = 2) # trt
  lrt3 <- glmLRT(fit2, coef = 3) # cov
  
  pvalue.trt.cov <- lrt2$table$PValue
  
  
  pvalue.cov <- lrt3$table$PValue
  
  pvalue.trt.nocov <- lrt$table$PValue
  
  aic.nocov <- AIC_glmfit(sim.data$counts, fit)
  
  aic.cov <- AIC_glmfit(sim.data$counts, fit2)
  
  ebp.cov <- pval.hist.modified0(pvalue.cov)
  g.ebp.cov <- pval_hist_grenander(pvalue.cov)
  
  ebp_cov <- laply(1:J, function(j)ifelse(ebp.cov$h.ebp[j] < .5,  1, 0))
  #   sum(abs(beta.ind)*ebp_cov)
  #   p.beta
  
  if ((sum(ebp_cov ==0)!= 0) & (sum(ebp_cov ==1)!=0))
  {
    ## code to fit QL.fit for each set of genes (cov and nocov)
    # Design list with covariate in full model 
    design <- model.matrix(~rep(1:2, each = K)+c(x[1,],x[2,]))
    counts_ebp <- counts[ebp_cov ==1, ]
    counts_ebp <- DGEList(counts_ebp, group = design[,2])
    counts_ebp <- calcNormFactors(counts_ebp)
    counts_ebp <- estimateGLMCommonDisp(counts_ebp, design, verbose = TRUE)
    counts_ebp <- estimateGLMTrendedDisp(counts_ebp, design)
    counts_ebp <- estimateGLMTagwiseDisp(counts_ebp, design)
    fit2_ebp <- glmFit(counts_ebp, design)
    lrt2_ebp <- glmLRT(fit2_ebp, coef = 2) # trt
    lrt3_ebp <- glmLRT(fit2_ebp, coef = 3) # cov
    
    pvalue.trt.cov_ebp <- lrt2_ebp$table$PValue
    
    design <- model.matrix(~rep(1:2, each = K))
    counts_ebp0 <- DGEList(counts[ebp_cov==0,], group = design[,2])
    counts_ebp0 <- calcNormFactors(counts_ebp0)
    counts_ebp0 <- estimateGLMCommonDisp(counts_ebp0, design, verbose = TRUE)
    counts_ebp0 <- estimateGLMTrendedDisp(counts_ebp0, design)
    counts_ebp0 <- estimateGLMTagwiseDisp(counts_ebp0, design)
    fit_ebp0 <- glmFit(counts_ebp0, design)
    lrt_ebp0 <- glmLRT(fit_ebp0, coef = 2)
    
    pvalue.trt.nocov_ebp <- lrt_ebp0$table$PValue
    
    pvalue.trt.ebp <- rep(0, J)
    pvalue.trt.ebp[which(ebp_cov==1)] <- pvalue.trt.cov_ebp
    pvalue.trt.ebp[which(ebp_cov==0)] <- pvalue.trt.nocov_ebp
    
  }
  
  if ((sum(ebp_cov ==0)== 0)) pvalue.trt.ebp <- pvalue.trt.cov
  
  
  if ((sum(ebp_cov ==1)== 0)) pvalue.trt.ebp <- pvalue.trt.nocov
  # gene classification when using grenander estimator
  
  
  gebp_cov <- laply(1:J, function(j)
    ifelse(g.ebp.cov$h.ebp[j] < .5,  1, 0))
  ## code to fit QL.fit for each set of genes (cov and nocov)
  if ((sum(gebp_cov ==0)!= 0) & (sum(gebp_cov ==1)!=0))
  {
    ## code to fit QL.fit for each set of genes (cov and nocov)
    # Design list with covariate in full model 
    design <- model.matrix(~rep(1:2, each = K)+c(x[1,],x[2,]))
    counts_gebp <- counts[gebp_cov ==1, ]
    counts_gebp <- DGEList(counts_gebp, group = design[,2])
    counts_gebp <- calcNormFactors(counts_gebp)
    counts_gebp <- estimateGLMCommonDisp(counts_gebp, design, verbose = TRUE)
    counts_gebp <- estimateGLMTrendedDisp(counts_gebp, design)
    counts_gebp <- estimateGLMTagwiseDisp(counts_gebp, design)
    fit2_gebp <- glmFit(counts_gebp, design)
    lrt2_gebp <- glmLRT(fit2_gebp, coef = 2) # trt
    lrt3_gebp <- glmLRT(fit2_gebp, coef = 3) # cov
    
    pvalue.trt.cov_gebp <- lrt2_gebp$table$PValue
    
    design <- model.matrix(~rep(1:2, each = K))
    counts_gebp0 <- DGEList(counts[gebp_cov==0,], group = design[,2])
    counts_gebp0 <- calcNormFactors(counts_gebp0)
    counts_gebp0 <- estimateGLMCommonDisp(counts_gebp0, design, verbose = TRUE)
    counts_gebp0 <- estimateGLMTrendedDisp(counts_gebp0, design)
    counts_gebp0 <- estimateGLMTagwiseDisp(counts_gebp0, design)
    fit_gebp0 <- glmFit(counts_gebp0, design)
    lrt_gebp0 <- glmLRT(fit_gebp0, coef = 2)
    
    pvalue.trt.nocov_gebp <- lrt_gebp0$table$PValue
    
    pvalue.trt.g.ebp <- rep(0, J)
    pvalue.trt.g.ebp[which(gebp_cov==1)] <- pvalue.trt.cov_gebp
    pvalue.trt.g.ebp[which(gebp_cov==0)] <- pvalue.trt.nocov_gebp
    
  }
  
  if ((sum(gebp_cov ==0)== 0)) pvalue.trt.g.ebp <- pvalue.trt.cov
  
  
  if ((sum(gebp_cov ==1)== 0)) pvalue.trt.g.ebp <- pvalue.trt.nocov
  # gene classification when using aic
  
  
  aic_cov <- laply(1:J, function(j)
    ifelse(aic.nocov[j] > aic.cov[j],  1, 0))
  ## code to fit QL.fit for each set of genes (cov and nocov)
  if ((sum(aic_cov ==0)!= 0) & (sum(aic_cov ==1)!=0))
  { ## code to fit QL.fit for each set of genes (cov and nocov)
    # Design list with covariate in full model 
    design <- model.matrix(~rep(1:2, each = K)+c(x[1,],x[2,]))
    counts_aic <- counts[aic_cov ==1, ]
    counts_aic <- DGEList(counts_aic, group = design[,2])
    counts_aic <- calcNormFactors(counts_aic)
    counts_aic <- estimateGLMCommonDisp(counts_aic, design, verbose = TRUE)
    counts_aic <- estimateGLMTrendedDisp(counts_aic, design)
    counts_aic <- estimateGLMTagwiseDisp(counts_aic, design)
    fit2_aic <- glmFit(counts_aic, design)
    lrt2_aic <- glmLRT(fit2_aic, coef = 2) # trt
    lrt3_aic <- glmLRT(fit2_aic, coef = 3) # cov
    
    pvalue.trt.cov_aic <- lrt2_aic$table$PValue
    
    design <- model.matrix(~rep(1:2, each = K))
    counts_aic0 <- DGEList(counts[aic_cov==0,], group = design[,2])
    counts_aic0 <- calcNormFactors(counts_aic0)
    counts_aic0 <- estimateGLMCommonDisp(counts_aic0, design, verbose = TRUE)
    counts_aic0 <- estimateGLMTrendedDisp(counts_aic0, design)
    counts_aic0 <- estimateGLMTagwiseDisp(counts_aic0, design)
    fit_aic0 <- glmFit(counts_aic0, design)
    lrt_aic0 <- glmLRT(fit_aic0, coef = 2)
    
    pvalue.trt.nocov_aic <- lrt_aic0$table$PValue
    pvalue.trt.aic <- rep(0, J)
    pvalue.trt.aic[which(aic_cov==1)] <- pvalue.trt.cov_aic
    pvalue.trt.aic[which(aic_cov==0)] <- pvalue.trt.nocov_aic
  }
  if ((sum(aic_cov ==0)== 0)) pvalue.trt.aic <- pvalue.trt.cov
  
  
  if ((sum(aic_cov ==1)== 0)) pvalue.trt.aic <- pvalue.trt.nocov
  
  ## oracle
  

  oracle_cov  <- laply(1:J, function(j)
    ifelse(beta.ind[j]==0,  0, 1))
  ## code to fit QL.fit for each set of genes (cov and nocov)
  if ((sum(oracle_cov ==0)!= 0) & (sum(oracle_cov ==1)!=0))
  {
    design <- model.matrix(~rep(1:2, each = K)+c(x[1,],x[2,]))
    counts_oracle <- counts[oracle_cov ==1, ]
    counts_oracle <- DGEList(counts_oracle, group = design[,2])
    counts_oracle <- calcNormFactors(counts_oracle)
    counts_oracle <- estimateGLMCommonDisp(counts_oracle, design, verbose = TRUE)
    counts_oracle <- estimateGLMTrendedDisp(counts_oracle, design)
    counts_oracle <- estimateGLMTagwiseDisp(counts_oracle, design)
    fit2_oracle <- glmFit(counts_oracle, design)
    lrt2_oracle <- glmLRT(fit2_oracle, coef = 2) # trt
    lrt3_oracle <- glmLRT(fit2_oracle, coef = 3) # cov
    
    pvalue.trt.cov_oracle <- lrt2_oracle$table$PValue
    
    design <- model.matrix(~rep(1:2, each = K))
    counts_oracle0 <- DGEList(counts[oracle_cov==0,], group = design[,2])
    counts_oracle0 <- calcNormFactors(counts_oracle0)
    counts_oracle0 <- estimateGLMCommonDisp(counts_oracle0, design, verbose = TRUE)
    counts_oracle0 <- estimateGLMTrendedDisp(counts_oracle0, design)
    counts_oracle0 <- estimateGLMTagwiseDisp(counts_oracle0, design)
    fit_oracle0 <- glmFit(counts_oracle0, design)
    lrt_oracle0 <- glmLRT(fit_oracle0, coef = 2)
    
    pvalue.trt.nocov_oracle <- lrt_oracle0$table$PValue
    pvalue.trt.oracle <- rep(0, J)
    pvalue.trt.oracle[which(abs(beta.ind)==1)] <- pvalue.trt.cov_oracle
    pvalue.trt.oracle[which(beta.ind==0)] <- pvalue.trt.nocov_oracle
  }
  if ((sum(oracle_cov ==1)== 0)) pvalue.trt.oracle <- pvalue.trt.cov
  
  
  if ((sum(oracle_cov ==1)== 0)) pvalue.trt.oracle <- pvalue.trt.nocov
  
  

  #########################AAA methods with AHE and Grenander Estimate###
  
  ebp.cov <- pval.hist.modified0(pvalue.cov)
  ebp.trt.nocov <- pval.hist.modified0(pvalue.trt.nocov)
  ebp.trt.cov <- pval.hist.modified0(pvalue.trt.cov)
  aaa.trt <- ebp.trt.nocov$h.ebp*ebp.cov$h.ebp +ebp.trt.cov$h.ebp*(1-ebp.cov$h.ebp)
  
  
  g.ebp.cov <- pval_hist_grenander(pvalue.cov)
  g.ebp.trt.nocov <- pval_hist_grenander(pvalue.trt.nocov)
  g.ebp.trt.cov <- pval_hist_grenander(pvalue.trt.cov)
  g.aaa.trt <- g.ebp.trt.nocov$h.ebp*g.ebp.cov$h.ebp +g.ebp.trt.cov$h.ebp*(1-g.ebp.cov$h.ebp)
  
  
  # ROC curves for 7 methods
  
  
  auc.nocov <- auc_out(pvalue.trt.nocov)[1]
  auc.nocov2 <- auc_out(pvalue.trt.nocov)[2]
  
  auc.cov <- auc_out(pvalue.trt.cov)[1]
  auc.cov2 <- auc_out(pvalue.trt.cov)[2]
  
  auc.ebp <- auc_out(pvalue.trt.ebp)[1]
  auc.ebp2 <- auc_out(pvalue.trt.ebp)[2]
  
  auc.g.ebp <- auc_out(pvalue.trt.g.ebp)[1]
  auc.g.ebp2 <- auc_out(pvalue.trt.g.ebp)[2]
  
  auc.aic <- auc_out(pvalue.trt.aic)[1]
  auc.aic2 <- auc_out(pvalue.trt.aic)[2]
  
  
  auc.aaa <- auc_out(aaa.trt)[1]
  auc.aaa2 <- auc_out(aaa.trt)[2]
  
  auc.g.aaa <- auc_out(g.aaa.trt)[1]
  auc.g.aaa2 <- auc_out(g.aaa.trt)[2]
  
  auc.oracle <- auc_out(pvalue.trt.oracle)[1]
  auc.oracle2 <- auc_out(pvalue.trt.oracle)[2]
  
  #qvalue, FDP for first five methods using Storey FDR with estimate.m0 by Nettleton et. al. 2006
  

  fdp.nocov <- fdp(pvalue.trt.nocov)
  fdp.cov <- fdp(pvalue.trt.cov)
  fdp.ebp <- fdp(pvalue.trt.ebp)
  fdp.g.ebp <- fdp(pvalue.trt.g.ebp)
  fdp.aic <- fdp(pvalue.trt.aic)
  fdp.aaa <- fdr_ebp(aaa.trt)
  fdp.g.aaa <- fdr_ebp(g.aaa.trt)
  fdp.oracle <- fdp(pvalue.trt.oracle)
  
  res <- list(pvalue.trt.nocov=pvalue.trt.nocov, 
              pvalue.trt.cov = pvalue.trt.cov, 
              pvalue.trt.ebp = pvalue.trt.ebp, 
              pvalue.trt.g.ebp = pvalue.trt.g.ebp, 
              pvalue.trt.aic = pvalue.trt.aic, 
              aaa.trt = aaa.trt, 
              g.aaa.trt = g.aaa.trt, 
              pvalue.trt.oracle = pvalue.trt.oracle,
            
              auc.nocov = auc.nocov,
              auc.cov = auc.cov,
              auc.ebp = auc.ebp,
              auc.g.ebp = auc.g.ebp,
              auc.aic = auc.aic,
              auc.aaa = auc.aaa,
              auc.g.aaa = auc.g.aaa,
              auc.oracle = auc.oracle,
              
              auc.nocov2 = auc.nocov2,
              auc.cov2 = auc.cov2 ,
              auc.ebp2 = auc.ebp2 ,
              auc.g.ebp2 = auc.g.ebp2,
              auc.aic2 = auc.aic2 ,
              auc.aaa2 = auc.aaa2,
              auc.g.aaa2 = auc.g.aaa2,
              auc.oracle2 = auc.oracle2,
              
              fdp.nocov = fdp.nocov,
              fdp.cov = fdp.cov ,
              fdp.ebp = fdp.ebp ,
              fdp.g.ebp = fdp.g.ebp,
              fdp.aic = fdp.aic ,
              fdp.aaa = fdp.aaa,
              fdp.g.aaa =fdp.g.aaa,
              fdp.oracle = fdp.oracle)          
  return(res)
}

# j <- 2;i <- 63 

edger_out_20_25 <- llply(1:length(i.beta), function(j){
  out1 <- laply(1:n.sim, function(i){
    edger_fit_out <- edger_fit(p.beta, i.beta[j], e.beta[j], S, L, U)
    pathsave <- paste(dir.pbeta1, 
                      "/p.beta_",
                      p.beta, 
                      "i.beta_",
                      i.beta[j],
                      "e.beta_",
                      e.beta[j],
                      "m_",
                      i,
                      ".RData",sep = "")
    save(edger_fit_out, file = pathsave)
    print(paste("pbeta = ", p.beta,", i.beta = ", i.beta[j],  ", m = ", i))
    return(c(auc.nocov = edger_fit_out$auc.nocov, auc.cov = edger_fit_out$auc.cov, 
             auc.ebp = edger_fit_out$auc.ebp, auc.g.ebp = edger_fit_out$auc.g.ebp, 
             auc.aic = edger_fit_out$auc.aic, auc.aaa = edger_fit_out$auc.aaa, 
             auc.g.aaa = edger_fit_out$auc.g.aaa, auc.oracle = edger_fit_out$auc.oracle,
             
             auc.nocov2 = edger_fit_out$auc.nocov2, auc.cov2 = edger_fit_out$auc.cov2, 
             auc.ebp2 = edger_fit_out$auc.ebp2, auc.g.ebp2 = edger_fit_out$auc.g.ebp2, 
             auc.aic2 = edger_fit_out$auc.aic2, auc.aaa2 = edger_fit_out$auc.aaa2, 
             auc.g.aaa2 = edger_fit_out$auc.g.aaa2, auc.oracle2 = edger_fit_out$auc.oracle2, 
             
             fdp.nocov = edger_fit_out$fdp.nocov, fdp.cov = edger_fit_out$fdp.cov, 
             fdp.ebp = edger_fit_out$fdp.ebp, fdp.g.ebp = edger_fit_out$fdp.g.ebp, 
             fdp.aic = edger_fit_out$fdp.aic, fdp.aaa = edger_fit_out$fdp.aaa, 
             fdp.g.aaa = edger_fit_out$fdp.g.aaa, fdp.oracle = edger_fit_out$fdp.oracle
    ))
  })
  
  out1
} )

head(edger_out_20_25[[1]])
head(edger_out_20_25[[2]])
save(edger_out_20_25, file = paste(dir.pbeta1, "/edger_out_20_25.RData", sep = ""))