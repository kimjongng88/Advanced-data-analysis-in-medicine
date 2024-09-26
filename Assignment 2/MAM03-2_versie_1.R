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
renalt$uprotein = as.factor(as.character(renalt$uprotein))

renalt$uprotein = as.numeric(renalt$uprotein)

table1(~dgf + acclft + aantalre + creat + predias + prac + uprotein + cregsh + gsurv|gstatus,data = renalt, overall = F, extra.col = list('P-value' = pvalue))

acclft <- summary(glm(gstatus ~ acclft, family = "binomial", data = renalt, na.action = na.omit))
dgf <- summary(glm(gstatus ~ dgf, family = "binomial", data = renalt, na.action = na.omit))
aantalre <- summary(glm(gstatus ~ aantalre, family = "binomial", data = renalt, na.action = na.omit))
creat <- summary(glm(gstatus ~ creat, family = "binomial", data = renalt, na.action = na.omit))
predias <- summary(glm(gstatus ~ predias, family = "binomial", data = renalt, na.action = na.omit))
prac <- summary(glm(gstatus ~ prac, family = "binomial", data = renalt, na.action = na.omit))
uprotein <- summary(glm(gstatus ~ uprotein, family = "binomial", data = renalt, na.action = na.omit))
cregsh <- summary(glm(gstatus ~ cregsh, family = "binomial", data = renalt, na.action = na.omit))



acclft_surv <- summary(lm(gsurv ~ acclft, data = renalt, na.action = na.omit))
dgf_surv <- summary(lm(gsurv ~ dgf, data = renalt, na.action = na.omit))
aantalre_surv <- summary(lm(gsurv ~ aantalre, data = renalt, na.action = na.omit))
creat_surv <- summary(lm(gsurv ~ creat,data = renalt, na.action = na.omit))
predias_surv <- summary(lm(gsurv ~ predias, data = renalt, na.action = na.omit))
prac_surv <- summary(lm(gsurv ~ prac,  data = renalt, na.action = na.omit))
uprotein_surv <- summary(lm(gsurv ~ uprotein, data = renalt, na.action = na.omit))
cregsh_surv <- summary(lm(gsurv ~ cregsh, data = renalt, na.action = na.omit))

library(survival)
Surv(renalt$gsurv,renalt$gstatus)
plot(survfit(Surv(renalt$gsurv,renalt$gstatus)~1))

plot(lm(gsurv ~ acclft, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ dgf, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ aantalre, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ creat,data = renalt, na.action = na.omit))
plot(lm(gsurv ~ predias, data = renalt, na.action = na.omit))
plot(lm(gsurv ~ prac,  data = renalt, na.action = na.omit))
plot(lm(gsurv ~ cregsh, data = renalt, na.action = na.omit))

############
multiStatus <- glm(gstatus ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
summary(multiStatus) # show results

multiSurv <- glm(gsurv ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
summary(multiSurv) # show results



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
  #model = lm(gstatus ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
  #loglikelihood[j] = logLik(model)
}

sapply(lapply(ystar, unique), length)
lapply(ystar, unique)

modelTransY = lm(ystar ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt)
loglikelihood[j] = logLik(modelTrans)

plot(lambda, loglikelihood)
boxcox(y ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=renalt, plotit=TRUE)

##################
#   POLYNOMIAL   #
##################

## Status ##
modelPolyStatus = lm(gstatus ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh,
            data = renalt, na.action = na.exclude)

polyStatAcclft = lm(gstatus ~ poly(acclft,10), na.action = na.exclude, data = renalt)
polyStatDgf = lm(gstatus ~ poly(dgf,10), na.action = na.exclude, data = renalt)
polyStatAantalre = lm(gstatus ~ poly(aantalre,10), na.action = na.exclude, data = renalt)
polyStatCreat = lm(gstatus ~ poly(creat,10), na.action = na.exclude, data = renalt)
polyStatPredias = lm(gstatus ~ poly(predias,10), na.action = na.exclude, data = renalt)
polyStatPrac = lm(gstatus ~ poly(prac,10), na.action = na.exclude, data = renalt)
polyStatUprotein = lm(gstatus ~ poly(uprotein,10), na.action = na.exclude, data = renalt)
polyStatCregsh = lm(gstatus ~ poly(cregsh,10), na.action = na.exclude, data = renalt)

aic=c()
for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(acclft,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - acclft",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(dgf,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - dgf",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(aantalre,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - aantalre",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(creat,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - creat",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(predias,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - predias",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(prac,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - prac",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(uprotein,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - uprotein",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gstatus ~ poly(cregsh,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - cregsh",ylab="AIC", type="b")

## Survival ##
modelPolySurv = lm(gsurv ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh,
                   data = renalt, na.action = na.exclude)

polySurvAcclft = lm(gsurv ~ poly(acclft,10), na.action = na.exclude, data = renalt)
polySurvDgf = lm(gsurv ~ poly(dgf,10), na.action = na.exclude, data = renalt)
polySurvAantalre = lm(gsurv ~ poly(aantalre,10), na.action = na.exclude, data = renalt)
polySurvCreat = lm(gsurv ~ poly(creat,10), na.action = na.exclude, data = renalt)
polySurvPredias = lm(gsurv ~ poly(predias,10), na.action = na.exclude, data = renalt)
polySurvPrac = lm(gsurv ~ poly(prac,10), na.action = na.exclude, data = renalt)
polySurvUprotein = lm(gsurv ~ poly(uprotein,10), na.action = na.exclude, data = renalt)
polySurvCregsh = lm(gsurv ~ poly(cregsh,10), na.action = na.exclude, data = renalt)

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(acclft,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - acclft",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(dgf,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - dgf",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(aantalre,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - aantalre",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(creat,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - creat",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(predias,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - predias",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(prac,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - prac",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(uprotein,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - uprotein",ylab="AIC", type="b")

for (i in 1:15) {
  aic[i] = AIC(lm(gsurv ~ poly(cregsh,i),data=renalt))
}
plot(1:15,aic,xlab="order of the polynome - cregsh",ylab="AIC", type="b")


#############
##SMOOTHING##
#############

library(splines)
aic=c()
for ( i in 1 : 30) {
  aic[i]=AIC(lm(gsurv ~ ns(acclft,df=i), data=renalt))
}

aic_re=c()
for (i in 1:30) {
  aic_re[i] = AIC(lm(gsurv ~ ns(aantalre,df = i),data=renalt))
}

aic_cr=c()
for (i in 1:30) {
  aic_cr[i] = AIC(lm(gsurv ~ ns(creat,df = i),data=renalt))
}

aic_pre=c()
for (i in 1:30) {
  aic_pre[i] = AIC(lm(gsurv ~ ns(predias,i),data=renalt))
}

aic_prac = c()
for (i in 1:30) {
  aic_prac[i] = AIC(lm(gsurv ~ ns(prac,i),data=renalt))
}

aic_up=c()
for (i in 1:30) {
  aic_up[i] = AIC(lm(gsurv ~ ns(uprotein,i),data=renalt))
}

aic_creg = c()
for (i in 1:30) {
  aic_creg[i] = AIC(lm(gsurv ~ ns(cregsh,i),data=renalt))
}


plot(1:30,aic,xlab="order of the polynome",ylab="AIC", type="b", ylim = c(7240,7305))
lines(1:7,aic_re,col="red")
lines(1:30,aic_cr,col="blue")
lines(1:30,aic_pre,col="yellow")
lines(1:30,aic_prac,col="purple")
lines(1:1,aic_up,col="pink")
lines(1:30,aic_creg,col="orange")



aic=c()
for ( i in 1 : 50) {
  aic[i]=AIC(lm(gstatus ~ ns(acclft,df=i), data=renalt))
}

aic_re=c()
for (i in 1:7) {
  aic_re[i] = AIC(lm(gstatus ~ ns(aantalre,df = i),data=renalt))
}

plot(1:7,aic_re,xlab="order of the polynome",ylab="AIC", type="b")


aic_re


aic_cr=c()
for (i in 1:50) {
  aic_cr[i] = AIC(lm(gstatus ~ ns(creat,df = i),data=renalt))
}

aic_pre=c()
for (i in 1:50) {
  aic_pre[i] = AIC(lm(gstatus ~ ns(predias,i),data=renalt))
}

aic_prac = c()
for (i in 1:50) {
  aic_prac[i] = AIC(lm(gstatus ~ ns(prac,i),data=renalt))
}

aic_up=c()
for (i in 1:50) {
  aic_up[i] = AIC(lm(gstatus ~ ns(uprotein,i),data=renalt))
}

aic_creg = c()
for (i in 1:50) {
  aic_creg[i] = AIC(lm(gstatus ~ ns(cregsh,i),data=renalt))
}

aic_ns[2]


plot(1:50,aic,xlab="order of the polynome",ylab="AIC", type="b", ylim = c(400,520))
lines(1:7,aic_re,col="red")
lines(1:30,aic_cr,col="blue")
lines(1:30,aic_pre,col="yellow")
lines(1:30,aic_prac,col="purple")
lines(1:1,aic_up,col="pink")
lines(1:30,aic_creg,col="orange")


model2=lm(renalt$gsurv~ns(renalt$aantalre,df=3))
AIC(model2)
plot(renalt$gsurv,renalt$aantalre)
lines(renalt$gsurv,predict(model2),col=2,lwd=2)


model1 = glm(renalt$gstatus~renalt$uprotein,family=binomial)
model2 = glm(renalt$gstatus~ns(renalt$uprotein,df=8),family=binomial)
plot(renalt$uprotein,renalt$gstatus,ylab="event risk",xlab="prto")
lines(renalt$uprotein,predict(model1,type="resp"),col=2,lwd=2)
lines(renalt$uprotein,predict(model2,type="resp"),col=3,lwd=2)

model1 = glm(renalt$gstatus~renalt$acclft,family=binomial)
model2 = glm(renalt$gstatus~ns(renalt$acclft,df=2),family=binomial)
plot(renalt$acclft,renalt$gstatus,ylab="event risk",xlab="prto")
lines(renalt$acclft,predict(model1,type="resp"),col=2,lwd=2)
lines(renalt$acclft,predict(model2,type="resp"),col=3,lwd=2)

model1 = glm(renalt$gstatus~renalt$creat,family=binomial)
model2 = glm(renalt$gstatus~ns(renalt$creat,df=2),family=binomial)
plot(renalt$creat,renalt$gstatus,ylab="event risk",xlab="prto")
lines(renalt$creat,predict(model1,type="resp"),col=2,lwd=2)
lines(renalt$creat,predict(model2,type="resp"),col=3,lwd=2)


dgf, panel reactive, diastolic





