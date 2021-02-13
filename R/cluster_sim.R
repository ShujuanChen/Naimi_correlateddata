packages <- c("data.table","tidyverse","skimr","here","lmtest","sandwich",
              "geepack","lme4","RColorBrewer","numDeriv","wgeesel","MASS")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

remotes::install_github("ikosmidis/enrichwith")
library(enrichwith)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

## simulated BMI example to demonstrate MLE with clustered data
set.seed(123)
n=100000
id=rep(1:n,each=6)
rho=0
phi=1
s=3.25
exposure <- rep(1,length(id)) #rbinom(length(id),1,.5)
x=cbind(1,exposure)
beta=c(22,3.25)

cont_fun <- function(beta,rho,phi,sigma,x,id,corstr="independence"){ #N is the sample size, nsubj is cluster size;
  response=c()
  N=length(unique(id))
  i=1
  for(i in 1:N){
    mu_i <- x[which(id==i),]%*%beta 
    if(corstr=="independence"){
      R_i=diag(1,length(mu_i))
    }
    else if (corstr=="exchangeable"){
      R_i=diag(1,length(mu_i))+matrix(rho,ncol=length(mu_i),nrow=length(mu_i))-diag(rho,length(mu_i))
    }
    else {
      H <- abs(outer(1:length(mu_i),1:length(mu_i), "-")) 
      R_i=rho^H
    }
    
    S_i=rep(sigma,length(mu_i))*diag(length(mu_i))
    m=phi*(sqrt(S_i)%*%R_i%*%sqrt(S_i))
    response <- c(response,mvrnorm(n = 1, mu=mu_i, Sigma=m))
  }
  data=cbind(id,x[,2:length(beta)],response)
  data=data.frame(data)
  
  return(data)
}

sim_dat <- cont_fun(beta,rho,phi,s,x,id,corstr="independence")

names(sim_dat)[2] <- "exposure"

head(sim_dat)
mean(sim_dat$exposure)
mean(sim_dat$response)
var(sim_dat$response)

sim_summary <- summary(aov(response ~ as.factor(id),data=sim_dat)) # type II SS
icc <- sim_summary[[1]][1,2]/sum(sim_summary[[1]][,2])
icc