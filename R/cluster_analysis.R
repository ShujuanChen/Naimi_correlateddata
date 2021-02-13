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

## empirical lead example to demonstrate robust variance and bootstrap
## read in data
#### data obtained here: https://content.sph.harvard.edu/fitzmaur/ala/tlc.txt
# lead_data <- read_csv("./longitudinal_lead_data.csv")
# Each row of the data set contains the following 6 variables:
# ID, Treatment Group, Lead Level Week 0, Lead Level Week 1, Lead Level Week 4, Lead Level Week 6. 

lead_data <- read_table("https://content.sph.harvard.edu/fitzmaur/ala/tlc.txt", col_names = F, skip=29)
names(lead_data) <- c("ID", "treatment", "L0", "L1", "L4", "L6") 

lead_data <- gather(lead_data, week, lead_value, L0:L6, factor_key=TRUE) %>% 
  mutate(week=as.numeric(gsub("L", "", week)),
         treatment=as.numeric(treatment=="A")) %>% 
  arrange(ID,week)

lead_data %>% print(n=8)

lead_data %>% 
  group_by(treatment) %>% 
  summarize(meanLead=mean(lead_value),
            sdBMI=sd(lead_value),
            numTx=length(lead_value))

lead_data %>% ggplot(.) +
  geom_line(aes(x=week,
                y=lead_value,
                group=as.factor(ID), 
                color=as.factor(treatment)),
            size=.5,alpha=.5) +
  geom_point(aes(x=week,
                 y=lead_value,
                 group=as.factor(ID), 
                 color=as.factor(treatment)),
             size=1.5,alpha=.5) +
  xlab("Week Since Randomization") +
  ylab("Blood Lead Concentration (ug/dL)") +
  labs(color="Treatment") +
  scale_colour_grey(start = 0, end = .8) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(.01,0))


## compute intracluster correlation coefficient
lead_summary <- summary(aov(lead_value ~ as.factor(ID),data=lead_data)) # type II SS

icc <- lead_summary[[1]][1,2]/sum(lead_summary[[1]][,2])

icc

## fit model for lead ignoring clustering
mod1 <- glm(lead_value ~ treatment, data=lead_data,family=gaussian(link = "identity"))
summary(mod1)$coefficients

coefci(mod1, level = 0.95)[2,]


## use robust variance for clustering
coeftest(mod1,vcov=vcovCL(mod1,type="HC3",cluster=~ID))[2,2]
coefci(mod1,vcov=vcovCL(mod1,type="HC3",cluster=~ID), level = 0.95)[2,]

# gives SE identical to stata
coeftest(mod1,vcov=vcovCL(mod1,type="HC0"))[2,2]
coefci(mod1,vcov=vcovCL(mod1,type="HC0"), level = 0.95)[2,]

# gives SE identical to GEE below
coeftest(mod1,vcov=vcovCL(mod1, type="HC0", cluster=~practice, cadjust=F))[2,2]
coefci(mod1,vcov=vcovCL(mod1, type="HC0", cluster=~practice, cadjust=F), level = 0.95)[2,]


## use clustered bootstrap
set.seed(123)
boot_func <- function(boot_num){
  clusters <- as.numeric(names(table(lead_data$ID)))
  index <- sample(1:length(clusters), length(clusters), replace=TRUE)
  bb <- table(clusters[index])
  boot <- NULL
  for(zzz in 1:max(bb)){
    cc <- lead_data[lead_data$ID %in% names(bb[bb %in% c(zzz:max(bb))]),]
    cc$b_ID<-paste0(cc$ID,zzz)
    boot <- rbind(boot, cc)
  }
  
  mod1 <- glm(lead_value ~ treatment, data=boot,family=gaussian(link = "identity"))
  res <- cbind(boot_num,coef(mod1)[2])
  return(res)
}

boot_res <- lapply(1:200, function(x) boot_func(x))
boot_res <- as_tibble(do.call(rbind,boot_res))
names(boot_res)[2] <- "est"

ggplot(boot_res) + 
  geom_histogram(aes(x=est,y=..density..),
                 binwidth = .25, fill = "grey", color = "black") +
  geom_density(aes(est)) +
  xlab("Bootstrap Estimates") +
  ylab("Density") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0))

sd(boot_res[,2]) ## bootstrap standard error of the treatment estimate

LCL <- mod1$coefficient[2] - 1.96*sd(boot_res[,2])
UCL <- mod1$coefficient[2] + 1.96*sd(boot_res[,2])

mod1$coefficient[2]
LCL
UCL

## use GEE
mod1_ind <- geeglm(BMI~treatment, id=practice, data=lead_data, corstr="independence")
summary(mod1_ind)

mod1_exch <- geeglm(BMI~treatment, id=practice, data=lead_data, corstr="exchangeable")
summary(mod1_exch)

mod1_unstr <- geeglm(BMI~treatment, id=practice, data=lead_data, corstr="unstructured")
summary(mod1_unstr)

## use LMM
mod1_lmm1 <- lmer(BMI~treatment + (1 | practice), data=lead_data)
summ_mod1_lmm1 <- summary(mod1_lmm1)
summ_mod1_lmm1

mean(coef(mod1_lmm1)$practice[,1])


## final table

#ignoring
summary(mod1)$coefficients[2,1:2]

#robust var
coeftest(mod1,vcov=vcovHC(mod1,type="const",cluster=practice))[2,2]
coefci(mod1,vcov=vcovHC(mod1,type="HC3",cluster=practice), level = 0.95)[2,]

# robust var: gives SE identical to stata
coeftest(mod1,vcov=vcovCL)[2,2]
coefci(mod1,vcov=vcovCL, level = 0.95)[2,]

################### playing around #####################

vcov(mod1)
mod1$residuals%*%mod1$residuals

set.seed(1)
x <- c(1:4, 7)
y <- c(5 + rnorm(4,sd = 1.2), 35)
plot(x, y)

m <- lm(y ~ x)
summary(m)

mm <- glm(y~x)
summary(mm)


s2 <- deviance(mm)/(length(x)-2) ## recall: deviance generalizes residual sum 
## of squares for GLM
X <- model.matrix(mm)

# recall: if A is a square matrix, solve(A) returns inverse of A
V_hat <- solve(t(X) %*% X)
meat <- (t(X) %*% (s2*diag(length(x))) %*% X)

vce <- V_hat %*% meat %*% V_hat
sqrt(diag(vce))

# boot
sd(boot_res[,2]) ## standard error of the treatment estimate

LCL <- mod1$coefficient[2] - 1.96*sd(boot_res[,2])
UCL <- mod1$coefficient[2] + 1.96*sd(boot_res[,2])

mod1$coefficient[2]
LCL
UCL

# GEE
summary(mod1_ind)
summary(mod1_exch)
summary(mod1_unstr)

# lmem
summ_mod1_lmm1
mean(coef(mod1_lmm1)$practice[,1])