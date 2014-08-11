library("edgeR");library("plyr");library("fdrtool");library("AUC"); library("maps") ;library("fields")
I <- 2; J <- 1000
K <- 20
DE <- round(J*.2)
EE <- J - DE
S <- 1.25
L <- 0.1
U <- 0.5
p.beta <- 0
i.beta <- c(0.1, 1)
e.beta <- c(0.5, 1.5)
n.sim <- 100
mainDir <- "/home/ntyet/research/modelselection" # linux server
mainDir1 <- paste("/home/ntyet/research/K_", K, sep = "")  # linux server
dir.create(mainDir1, showWarnings = FALSE)
#mainDir <- "/run/user/1000/gvfs/smb-share:server=smb.stat.iastate.edu,share=ntyet/research/modelselection" # Linux laptop
pbeta1 <- "pbeta_0"
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
#source(paste(dir.source, "/QuasiSeq_Method_CompareFDR_BH_EBP_AHB_m0/fdrtool_1.2.10/fdrtool/R/ecdf.pval.R",sep =""))
source(paste(dir.source, "/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/QL.fit.R",sep=""))
source(paste(dir.source, "/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/NBDev.R",sep =""))
source(paste(dir.source, "/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/PoisDev.R",sep =""))
source(paste(dir.source, "/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/QL.results.R",sep =""))

# Load RFI count.mean and NBdisp from model0RFI.line

# 
# load(file = "U:/R/RA/Data/Additional Plot/Model0.line.rfi.RData")
# load(file = "U:/R/RA/Data/Additional Plot/Model0.result.line.rfi.RData")

#load(file = "U:/R/RA/Data/Additional Plot/Model0.line.rfi.RData")
#load(file = "U:/R/RA/Data/Additional Plot/Model0.result.line.rfi.RData")
# 
# load(file ="/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet//R/RA/Data/Additional Plot/Model0.line.rfi.RData")
# load(file ="/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet//R/RA/Data/Additional Plot/Model0.result.line.rfi.RData")
load(file = paste(dir.source, "/Model0.line.rfi.RData",sep = ""))
load(file = paste(dir.source,"/Model0.result.line.rfi.RData",sep = ""))

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

## Function to calculate AIC of the QL.fit model 

AIC_QL <- function(counts,QL.fit.object){
  n <- dim(counts)[2]
  m <- dim(counts)[1]
  disp <- 1/QL.fit.object$NB.disp
  den.df <- QL.fit.object$den.df
  phi.hat.dev <- QL.fit.object$phi.hat.dev
  p <- n - den.df
  dev <- phi.hat.dev*den.df
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
    if (sum(fdr.temp<=0.05)==0){
      out <- 0
    } else{
      out <- sum((fdr.temp<=0.05)&(1:J>DE))/sum(fdr.temp<=0.05)
    } 
  out
}

##Sim counts data

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


sim_QLfit <- function(p.beta, i.beta, e.beta, S, L, U){
  
  sim.data <- sim_counts(p.beta, i.beta, e.beta, S, L, U)
  counts <- sim.data$counts
  beta.ind <- sim.data$beta.ind
  trt.pos <- sim.data$trt.pos
  x <- sim.data$x
  design.list <- vector("list",2)
  design.list[[1]] <- rep(1:2, each = K)
  design.list[[2]] <- rep(1, ncol(counts))
  fit <- QL.fit(counts, design.list, 
                Model = "NegBin",
                print.progress=FALSE)
  
  result.fit <- QL.results(fit, Plot=FALSE)
  
  #############################################################
  
  # Design list with covariate in full model 
  design.list <- vector("list", 3)
  x1 <- as.factor(rep(1:2, each = K))
  design.list[[1]] <- model.matrix(~ x1 + c(x[1,],x[2,]))
  design.list[[2]] <- model.matrix(~ x1) # test for covariate
  design.list[[3]] <- model.matrix(~ c(x[1,],x[2,])) # test for treatment effect
  test.mat <- rbind(1:2, c(1,3))
  row.names(test.mat) <- c("Covariate", "Treatment")
  fit.2 <- QL.fit(counts, design.list, 
                  test.mat,
                  Model="NegBin", 
                  print.progress=FALSE)
  
  result.fit2 <- QL.results(fit.2,Plot= FALSE)
  
  pvalue.cov <- result.fit2$P.values[[3]][,1]
  pvalue.trt.nocov <- result.fit$P.values[[3]][,1]
  pvalue.trt.cov <- result.fit2$P.values[[3]][,2]
  
  ebp.cov <- pval.hist.modified0(pvalue.cov)
  g.ebp.cov <- pval_hist_grenander(pvalue.cov)
  
  aic.nocov <- AIC_QL(counts, fit)
  aic.cov <- AIC_QL(counts, fit.2)
  
  pvalue.trt.ebp <- laply(1:J, function(j)
    ifelse(ebp.cov$h.ebp[j] < .5,  pvalue.trt.cov[j], pvalue.trt.nocov[j]))

  pvalue.trt.g.ebp <- laply(1:J, function(j)
    ifelse(g.ebp.cov$h.ebp[j] < .5,  pvalue.trt.cov[j], pvalue.trt.nocov[j]))
  
  pvalue.trt.aic <- laply(1:J, function(j)
    ifelse(aic.nocov[j] > aic.cov[j],  pvalue.trt.cov[j], pvalue.trt.nocov[j]))
    
    pvalue.trt.oracle <- laply(1:J, function(j)
      ifelse(beta.ind[j]!=0,  pvalue.trt.cov[j], pvalue.trt.nocov[j]))
  
  
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
  res <- list(beta.ind = beta.ind, 
              trt.pos = trt.pos,
              counts=counts,
              x = c(x[1,], x[2,]),
              S = S, 
              U= U, 
              L = L, 
              
              p.beta = p.beta, 
              i.beta = i.beta, 
              e.beta = e.beta,
              
              pvalue.trt.nocov=pvalue.trt.nocov, 
              pvalue.trt.cov = pvalue.trt.cov, 
              pvalue.trt.ebp = pvalue.trt.ebp, 
              pvalue.trt.g.ebp = pvalue.trt.g.ebp, 
              pvalue.trt.aic = pvalue.trt.aic, 
              aaa.trt = aaa.trt, 
              g.aaa.trt = g.aaa.trt, 
              pvalue.trt.oracle,
            
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


out_20_0 <- llply(1:length(i.beta), function(j){
  out1 <- laply(1:n.sim, function(i){
    sim1 <- sim_QLfit(p.beta, i.beta[j], e.beta[j], S, L, U)
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
    save(sim1, file = pathsave)
    print(paste("pbeta = ", p.beta,", i.beta = ", i.beta[j],  ", m = ", i))
    return(c(auc.nocov = sim1$auc.nocov, auc.cov = sim1$auc.cov, 
             auc.ebp = sim1$auc.ebp, auc.g.ebp = sim1$auc.g.ebp, 
             auc.aic = sim1$auc.aic, auc.aaa = sim1$auc.aaa, 
             auc.g.aaa = sim1$auc.g.aaa, auc.oracle = sim1$auc.oracle,
             
             auc.nocov2 = sim1$auc.nocov2, auc.cov2 = sim1$auc.cov2, 
             auc.ebp2 = sim1$auc.ebp2, auc.g.ebp2 = sim1$auc.g.ebp2, 
             auc.aic2 = sim1$auc.aic2, auc.aaa2 = sim1$auc.aaa2, 
             auc.g.aaa2 = sim1$auc.g.aaa2, auc.oracle2 = sim1$auc.oracle2, 
             
             fdp.nocov = sim1$fdp.nocov, fdp.cov = sim1$fdp.cov, 
             fdp.ebp = sim1$fdp.ebp, fdp.g.ebp = sim1$fdp.g.ebp, 
             fdp.aic = sim1$fdp.aic, fdp.aaa = sim1$fdp.aaa, 
             fdp.g.aaa = sim1$fdp.g.aaa, fdp.oracle = sim1$fdp.oracle
    ))
  })
  
  out1
} )

head(out_20_0[[1]])
head(out_20_0[[2]])
save(out_20_0, file = paste(dir.pbeta1, "/out_20_0.RData", sep = ""))

