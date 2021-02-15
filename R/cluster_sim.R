userLib <- "~/R/R_LIBS_USER"
.libPaths(userLib)
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

# install.packages("simstudy")
# library(simstudy)

## simulated BMI example to demonstrate MLE with clustered data
set.seed(123)

sim_res <- function(rho){
n=400
id=rep(1:n,each=5)
#rho=.75
phi=1
#s=runif(n,1.25,5)
s=rep(3.25,n)
exposure <- rbinom(length(id),1,.5) #rep(1,length(id)) #
x=cbind(1,exposure)
beta=c(22,3.25)

cont_fun <- function(beta,rho,phi,sigma,x,id,corstr){ #N is the sample size, nsubj is cluster size;
  response=c()
  N=length(unique(id))
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
    
    S_i=rep(sigma[i],length(mu_i))*diag(length(mu_i))
    m=phi*(sqrt(S_i)%*%R_i%*%sqrt(S_i))
    response <- c(response,mvrnorm(n = 1, mu=mu_i, Sigma=m))
  }
  data=cbind(id,x[,2:length(beta)],response)
  data=data.frame(data)
  
  return(data)
}

sim_dat <- cont_fun(beta,rho,phi,s,x,id,corstr="exchangeable")

names(sim_dat)[2] <- "exposure"
sim_dat <- sim_dat %>% 
  group_by(id) %>% 
  mutate(count=1,time=cumsum(count))
# 
# ggplot(sim_dat) + geom_point(aes(x=id,y=response))
# 
# ggplot(sim_dat) + 
#   geom_line(aes(x=time,
#                 y=response,
#                 group=as.factor(id),
#                 color=as.factor(exposure)),
#           size=.5,alpha=.5) +
#   geom_point(aes(x=time,
#                  y=response,
#                  group=as.factor(id), 
#                  color=as.factor(exposure)),
#              size=1.5,alpha=.5) +
#   xlab("Time Since T0") +
#   ylab("Outcome") +
#   labs(color="Exposure") +
#   scale_colour_grey(start = 0, end = .8) +
#   scale_y_continuous(expand=c(0,0)) +
#   scale_x_continuous(expand=c(.01,0))
# 
# head(sim_dat,12)
# mean(sim_dat$exposure)
# mean(sim_dat$response)
# var(sim_dat$response)

sim_summary <- summary(aov(response ~ as.factor(id),data=sim_dat)) # type II SS
icc <- sim_summary[[1]][1,2]/sum(sim_summary[[1]][,2])
icc

mod1 <- glm(response ~ exposure,data=sim_dat,family=gaussian(link="identity"))

res <- c(rho,icc,summary(mod1)$coefficients[2,1:2])

return(res)

}

vec <- rep(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),10)
res <- lapply(vec, function(x) sim_res(rho=x))
res <- do.call(rbind,res)

write.csv(res,"./data/sim_res.csv")




