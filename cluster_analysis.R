packages <- c("data.table","tidyverse","skimr","here","lmtest","sandwich",
              "geepack","lme4","RColorBrewer")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)


## read in data
cluster_trial <- read_csv("./cluster_trial_data.csv")

cluster_trial %>% print(n=5)

cluster_trial %>% 
  group_by(practice,treatment) %>% 
  summarize(meanBMI=mean(BMI),
            numTx=length(BMI))

cluster_trial %>% 
  group_by(treatment) %>% 
  summarize(meanBMI=mean(BMI),
            sdBMI=sd(BMI),
            numTx=length(BMI))

set.seed(123)
cluster_trial %>% ggplot(.) +
  geom_point(aes(x=jitter(treatment,factor=.25),
                 y=BMI,
                 color=as.factor(practice)),
             size=5) +
  xlab("Treatment") +
  theme(legend.position = "none")


## compute intracluster correlation coefficient
bmi_summary <- summary(aov(BMI ~ as.factor(practice),data=cluster_trial)) # type II SS

icc <- bmi_summary[[1]][1,2]/sum(bmi_summary[[1]][,2])

icc

## fit model for BMI ignoring clustering
mod1 <- glm(BMI ~ treatment, data=cluster_trial,family=gaussian(link = "identity"))

summary(mod1)$coefficients

coeftest(mod1)[2,]

coefci(mod1, level = 0.95)[2,]

mod1 <- glm(BMI ~ treatment, data=cluster_trial,family=gaussian(link = "identity"))

summary(mod1)$coefficients


## use robust variance for clustering
coeftest(mod1,vcov=vcovHC(mod1))[2,2]
coeftest(mod1,vcov=vcovHC(mod1,type="HC3",cluster=practice))[2,2]
coefci(mod1,vcov=vcovHC(mod1,type="HC3",cluster=practice), level = 0.95)[2,]

# gives SE identical to stata
coeftest(mod1,vcov=vcovCL)[2,2]
coefci(mod1,vcov=vcovCL, level = 0.95)[2,]



## use clustered bootstrap
set.seed(123)
boot_func <- function(boot_num){
  clusters <- as.numeric(names(table(cluster_trial$practice)))
  index <- sample(1:length(clusters), length(clusters), replace=TRUE)
  bb <- table(clusters[index])
  boot <- NULL
  for(zzz in 1:max(bb)){
    cc <- cluster_trial[cluster_trial$practice %in% names(bb[bb %in% c(zzz:max(bb))]),]
    cc$b_practice<-paste0(cc$practice,zzz)
    boot <- rbind(boot, cc)
  }
  
  mod1 <- glm(BMI ~ treatment, data=boot,family=gaussian(link = "identity"))
  res <- cbind(boot_num,coef(mod1)[2])
  return(res)
}

boot_res <- lapply(1:200, function(x) boot_func(x))
boot_res <- do.call(rbind,boot_res)

sd(boot_res[,2]) ## standard error of the treatment estimate

LCL <- mod1$coefficient[2] - 1.96*sd(boot_res[,2])
UCL <- mod1$coefficient[2] + 1.96*sd(boot_res[,2])

mod1$coefficient[2]
LCL
UCL

## use GEE
mod1_ind <- geeglm(BMI~treatment, id=practice, data=cluster_trial, corstr="independence")
summary(mod1_ind)

mod1_exch <- geeglm(BMI~treatment, id=practice, data=cluster_trial, corstr="exchangeable")
summary(mod1_exch)

mod1_unstr <- geeglm(BMI~treatment, id=practice, data=cluster_trial, corstr="unstructured")
summary(mod1_unstr)

## use LMM
mod1_lmm1 <- lmer(BMI~treatment + (1 | practice), data=cluster_trial)
summ_mod1_lmm1 <- summary(mod1_lmm1)
summ_mod1_lmm1

mean(coef(mod1_lmm1)$practice[,1])





