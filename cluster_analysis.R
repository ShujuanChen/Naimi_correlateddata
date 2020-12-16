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
lead_trial <- read_csv("./longitudinal_lead_data.csv")

lead_trial <- gather(lead_trial, week, lead_value, L0:L6, factor_key = TRUE) %>%
  mutate(week = as.numeric(gsub("L", "", week))) %>% 
  arrange(ID, week) %>% 
  rename(treatment=Treatment)

lead_trial %>% print(n=10)

lead_trial %>% 
  group_by(week,treatment) %>% 
  summarize(mean_lead=mean(lead_value),
            numTx=length(lead_value))

lead_trial %>% 
  group_by(treatment) %>% 
  summarize(mean_lead=mean(lead_value),
            sd_lead=sd(lead_value),
            numTx=length(lead_value))

lead_trial %>% ggplot(.) +
  geom_point(aes(x=week,
                 y=lead_value,
                 color=as.factor(treatment)),
             size=2,
             alpha=.5) +
  geom_line(aes(x=week,
                y=lead_value,
                group=ID,
                color=as.factor(treatment)),
            alpha=.5) +
  xlab("week")

lead_trial <- lead_trial %>% mutate(treatment=as.numeric(treatment=="A"))

## compute intracluster correlation coefficient
lead_summary <- summary(aov(lead_value ~ as.factor(ID),data=lead_trial)) # type II SS

icc <- lead_summary[[1]][1,2]/sum(lead_summary[[1]][,2])

icc

## fit model for BMI ignoring clustering
mod1 <- glm(lead_value ~ treatment, data=lead_trial,family=gaussian(link = "identity"))

summary(mod1)$coefficients

coeftest(mod1)[2,]

coefci(mod1, level = 0.95)[2,]


## use robust variance for clustering
coeftest(mod1,vcov=vcovHC(mod1))[2,2]
coeftest(mod1,vcov=vcovHC(mod1,type="HC3",cluster=ID))[2,2]
coefci(mod1,vcov=vcovHC(mod1,type="HC3",cluster=ID), level = 0.95)[2,]

# gives SE identical to stata
coeftest(mod1,vcov=vcovCL)[2,2]
coefci(mod1,vcov=vcovCL, level = 0.95)[2,]

## use clustered bootstrap
set.seed(123)
boot_func <- function(boot_num){
  clusters <- as.numeric(names(table(lead_trial$ID)))
  index <- sample(1:length(clusters), length(clusters), replace=TRUE)
  bb <- table(clusters[index])
  boot <- NULL
  for(zzz in 1:max(bb)){
    cc <- lead_trial[lead_trial$ID %in% names(bb[bb %in% c(zzz:max(bb))]),]
    cc$b_practice<-paste0(cc$practice,zzz)
    boot <- rbind(boot, cc)
  }
  
  mod1 <- glm(lead_value ~ treatment, data=boot,family=gaussian(link = "identity"))
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
mod1_ind <- geeglm(lead_value~treatment, id=ID, data=lead_trial, corstr="independence")
summary(mod1_ind)

mod1_exch <- geeglm(lead_value~treatment, id=ID, data=lead_trial, corstr="exchangeable")
summary(mod1_exch)

mod1_ar1 <- geeglm(lead_value~treatment, id=ID, data=lead_trial, corstr="ar1")
summary(mod1_ar1)

mod1_unstr <- geeglm(lead_value~treatment, id=ID, data=lead_trial, corstr="unstructured")
summary(mod1_unstr)

## use LMM
mod1_lmm1 <- lmer(lead_value~treatment + (1 | ID), data=lead_trial)
summ_mod1_lmm1 <- summary(mod1_lmm1)
summ_mod1_lmm1

mod1_lmm2 <- lmer(lead_value~treatment + (1 + treatment | ID), data=lead_trial)
summ_mod1_lmm2 <- summary(mod1_lmm2)
summ_mod1_lmm2

mean(coef(mod1_lmm2)$ID[,1])





