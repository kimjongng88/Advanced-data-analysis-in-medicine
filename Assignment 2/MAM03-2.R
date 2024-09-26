setwd("C:/Users/joey.moonen/Documents/temp/mi/MAM03-2")

library(haven)
library(table1)
library(MASS)

renalt <-  read_sav("renaltx.sav")
renalt <- na.omit(renalt)


########TABLES#######
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

renalt$dgf = as.factor(as.character(renalt$dgf))
renalt$aantalre = as.factor(as.character(renalt$aantalre))
renalt$uprotein = as.factor(as.character(renalt$uprotein))


table1(~dgf + acclft + aantalre + creat + predias + prac + uprotein + cregsh + gsurv|gstatus,data = renalt, overall = F, extra.col = list('P-value' = pvalue))

acclft <- summary(glm(gstatus ~ acclft, family = "binomial", data = renalt, na.action = na.omit))
acclft
dgf <- summary(glm(gstatus ~ dgf, family = "binomial", data = renalt, na.action = na.omit))
dgf
aantalre <- summary(glm(gstatus ~ aantalre, family = "binomial", data = renalt, na.action = na.omit))
aantalre
creat <- summary(glm(gstatus ~ creat, family = "binomial", data = renalt, na.action = na.omit))
creat
predias <- summary(glm(gstatus ~ predias, family = "binomial", data = renalt, na.action = na.omit))
predias
prac <- summary(glm(gstatus ~ prac, family = "binomial", data = renalt, na.action = na.omit))
prac
uprotein <- summary(glm(gstatus ~ uprotein, family = "binomial", data = renalt, na.action = na.omit))
uprotein
cregsh <- summary(glm(gstatus ~ cregsh, family = "binomial", data = renalt, na.action = na.omit))
cregsh

table1(~dgf + acclft + aantalre + creat + predias + prac + uprotein + cregsh + gstatus|gsurv,data = renalt, overall = F, extra.col = list('P-value' = pvalue))

acclft_surv <- summary(lm(gsurv ~ acclft, data = renalt, na.action = na.omit))
acclft_surv
dgf_surv <- summary(lm(gsurv ~ dgf, data = renalt, na.action = na.omit))
dgf_surv
aantalre_surv <- summary(lm(gsurv ~ aantalre, data = renalt, na.action = na.omit))
aantalre_surv
creat_surv <- summary(lm(gsurv ~ creat,data = renalt, na.action = na.omit))
creat_surv
predias_surv <- summary(lm(gsurv ~ predias, data = renalt, na.action = na.omit))
predias_surv
prac_surv <- summary(lm(gsurv ~ prac,  data = renalt, na.action = na.omit))
prac_surv
uprotein_surv <- summary(lm(gsurv ~ uprotein, data = renalt, na.action = na.omit))
uprotein_surv
cregsh_surv <- summary(lm(gsurv ~ cregsh, data = renalt, na.action = na.omit))
cregsh_surv

############
multiStatus <- glm(gstatus ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
summary(multiStatus) # show results

multiSurv <- glm(gsurv ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
summary(multiSurv) # show results

library(survival)
Surv(renalt$gsurv,renalt$gstatus)
plot(survfit(Surv(renalt$gsurv,renalt$gstatus)~1))

# Alleen dgf lineair -> oftewel, hoeven we niet voor te corrigeren 
plot(lm(gsurv ~ acclft, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ dgf, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ aantalre, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ creat,data = renalt, na.action = na.omit))
plot(lm(gsurv ~ predias, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ prac,  data = renalt, na.action = na.omit))
plot(lm(gsurv ~ uprotein,  data = renalt, na.action = na.omit))
plot(lm(gsurv ~ cregsh, data = renalt, na.action = na.omit))

##################
# TRANSFORMATION #
##################
sapply(lapply(renalt, unique), length)
lapply(renalt[c('acclft', 'dgf', 'aantalre', 'creat', 'predias', 'prac', 'uprotein', 'cregsh', 'gsurv', 'gstatus')], unique)

y = renalt$gsurv
geomean = exp(mean(log(y)))
lambda = seq(-2, 2, 0.01)
loglikelihood=c()

for (j in 1:length(lambda)) {
  
  ystar = (y^lambda[j] -1) / (lambda[j]*geomean^(lambda[j]-1))
  #model = lm(ystar ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
  #loglikelihood[j] = logLik(model)
}

sapply(lapply(ystar, unique), length)
lapply(ystar, unique)

modelTrans = lm(ystar ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
modelTrans
loglikelihood[j] = logLik(modelTrans)

plot(lambda, loglikelihood)
# maximale waarde van de correlatie coefficient (even checken waar hij precies zit < Joey + voor status maken)
boxcox(y ~ acclft + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt, plotit=TRUE)

##################
#   POLYNOMIAL   #
##################

## Status ##
modelPolyStatus = lm(gstatus ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh,
            data = renalt, na.action = na.exclude)
modelPolyStatus

polyStatAcclft = lm(gstatus ~ poly(acclft,10), na.action = na.exclude, data = renalt)
polyStatAcclft
polyStatCreat = lm(gstatus ~ poly(creat,10), na.action = na.exclude, data = renalt)
polyStatCreat
polyStatPredias = lm(gstatus ~ poly(predias,10), na.action = na.exclude, data = renalt)
polyStatPredias
polyStatPrac = lm(gstatus ~ poly(prac,10), na.action = na.exclude, data = renalt)
polyStatPrac


# Per blokje runnen. Grafiek laat ziet wat de optimale orde is van de polynoom. De twee laatste regels van het blok berekenen de
# optimale AIC voor je uit. Dus de waarde die je in je console krijgt is de optimale AIC
aic_acclftStat=c()
for (i in 1:15) {
  aic_acclftStat[i] = AIC(lm(gstatus ~ poly(acclft,i),data=renalt))
}
plot(1:15,aic_acclftStat,xlab="order of the polynome - acclft",ylab="AIC", type="b")
AICacclftStat <- AIC(lm(gstatus ~ poly(acclft,2),data=renalt))
AICacclftStat

aic_creatStat=c()
for (i in 1:15) {
  aic_creatStat[i] = AIC(lm(gstatus ~ poly(creat,i),data=renalt))
}
plot(1:15,aic_creatStat,xlab="order of the polynome - creat",ylab="AIC", type="b")
AICcreatStat <- AIC(lm(gstatus ~ poly(creat,14),data=renalt))
AICcreatStat

aic_prediasStat=c()
for (i in 1:15) {
  aic_prediasStat[i] = AIC(lm(gstatus ~ poly(predias,i),data=renalt))
}
plot(1:15,aic_prediasStat,xlab="order of the polynome - predias",ylab="AIC", type="b")
AICprediasstat <- AIC(lm(gstatus ~ poly(predias,1),data=renalt))
AICprediasstat

aic_pracStat=c()
for (i in 1:15) {
  aic_pracStat[i] = AIC(lm(gstatus ~ poly(prac,i),data=renalt))
}
plot(1:15,aic_pracStat,xlab="order of the polynome - prac",ylab="AIC", type="b")
AICpracstat <- AIC(lm(gstatus ~ poly(prac,1),data=renalt))
AICpracstat

## Survival ##
modelPolySurv = lm(gsurv ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh,
                   data = renalt, na.action = na.exclude)

polySurvAcclft = lm(gsurv ~ poly(acclft,10), na.action = na.exclude, data = renalt)
polySurvCreat = lm(gsurv ~ poly(creat,10), na.action = na.exclude, data = renalt)
polySurvPredias = lm(gsurv ~ poly(predias,10), na.action = na.exclude, data = renalt)
polySurvPrac = lm(gsurv ~ poly(prac,10), na.action = na.exclude, data = renalt)

aic_acclftSurv=c()
for (i in 1:15) {
  aic_acclftSurv[i] = AIC(lm(gsurv ~ poly(acclft,i),data=renalt))
}
plot(1:15,aic_acclftSurv,xlab="order of the polynome - acclft",ylab="AIC", type="b")
AICacclftsturv <- AIC(lm(gsurv ~ poly(acclft,2),data=renalt))
AICacclftsturv

aic_creatSurv=c()
for (i in 1:15) {
  aic_creatSurv[i] = AIC(lm(gsurv ~ poly(creat,i),data=renalt))
}
plot(1:15,aic_creatSurv,xlab="order of the polynome - creat",ylab="AIC", type="b")
AICcreatsurv <- AIC(lm(gsurv ~ poly(creat,2),data=renalt))
AICcreatsurv

aic_prediasSurv=c()
for (i in 1:15) {
  aic_prediasSurv[i] = AIC(lm(gsurv ~ poly(predias,i),data=renalt))
}
plot(1:15,aic_prediasSurv,xlab="order of the polynome - predias",ylab="AIC", type="b")
AICprediassurv <- AIC(lm(gsurv ~ poly(predias,1),data=renalt))
AICprediassurv

aic_pracSurv=c()
for (i in 1:15) {
  aic_pracSurv[i] = AIC(lm(gsurv ~ poly(prac,i),data=renalt))
}
plot(1:15,aic_pracSurv,xlab="order of the polynome - prac",ylab="AIC", type="b")
AICpracsurv <- AIC(lm(gsurv ~ poly(prac,1),data=renalt))
AICpracsurv


###########
# SPLINES #
###########

library(splines)

# gekozen voor 5 splines bij surv -> leidt tot 3 verschillende lineairen lijnen
# gekozen voor 10 splines bij stat -> leidt tot 4 verschillende lineairen lijnen
# survival
aic_creg = c()
for (i in 1:5) {
  aic_creg[i] = AIC(lm(gsurv ~ ns(cregsh,i),data=renalt))
}

plot(1:5,aic_creg,xlab="order of the polynome",ylab="AIC", type="b", ylim = c(7240,7305))
lines(1:5,aic_creg,col="orange")

# status  
aic_creg = c()
for (i in 1:10) {
  aic_creg[i] = AIC(lm(gstatus ~ ns(cregsh,i),data=renalt))
}

plot(1:10,aic_creg,xlab="order of the polynome",ylab="AIC", type="b", ylim = c(400,520))
lines(1:10,aic_creg,col="orange")


## Dus: gecorrigeerd voor alles behalve dgf.
## Getransformeerd voor alles
## Polynomiaal voor acclft; creat; predias en prac
## Splines voor cregsh


library('dplyr')
library('funModeling')
renalt %>%
  dplyr::summarise(meanLE=mean(gsurv,na.rm=TRUE),
                   medLE=median(gsurv,na.rm=TRUE),
                   sd=sd(gsurv,na.rm=TRUE),
                   iqr=IQR(gsurv,na.rm=TRUE),
                   Q1=quantile(gsurv,probs=0.25,na.rm=TRUE),
                   Q3=quantile(gsurv,probs=0.75),
                   n=n())

funModeling::profiling_num(renalt)

funModeling::plot_num(renalt)

library("skimr")
varlist <- c("n_missing","complete_rate")
renalt %>% 
  select(gsurv) %>% 
  skimr::skim_without_charts() %>%
  skimr::yank("numeric") %>%
  dplyr::select(-one_of(varlist))





