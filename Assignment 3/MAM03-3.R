
#install.packages('nlme')
library(readr)
library(nlme)
library(scales)
library(ggplot2)

marfan <- read_csv('marfan.csv')
marfan <- na.omit(marfan)

marfan$patnr <- as.factor(marfan$patnr)


##########
## LMEM ##
##########

# ---> ANOVA methode gebruikt om te bepalen welk type model het beste is. Hieruit blijkt:
# - Age: Random intercept and slope (optie 3)
# - Sex: Repeated measurements  (optie 1)
# mochten modellen beide hetzelfde moeten zijn, dan random intercept and slope (optie 3), want AIC van 1 en 3 verschilt maar 1 bij sex

### AGE ###
lmeAge1 = lme(diameter~age, random = (~1| patnr), data=marfan, , method="REML")
lmeAge2 = lme(diameter~age, random = (~0+age|patnr), data=marfan, , method="REML")
lmeAge3 = lme(diameter~age, random = (~1+age|patnr), data=marfan, , method="REML")
anova(lmeAge1, lmeAge2, lmeAge3)

summary(lmeAge3)

par(mfrow=c(2,3))
plot(marfan$age, resid(lmeAge3))
abline(h=0, lty=2)
lines(unique(marfan$age), tapply(resid(lmeAge3), marfan$age, mean), col=2,lwd=3)
tapply(resid(lmeAge3), marfan$age, function(x) lines(0:9, x, lty=2))
qqnorm(resid(lmeAge3))
qqline(resid(lmeAge3))
hist(lmeAge3$coefficients$random$patnr[,1],xlab="patient specific intercept",main="")
qqnorm(lmeAge3$coefficients$random$patnr[,1])
qqline(lmeAge3$coefficients$random$patnr[,1])
hist(lmeAge3$coefficients$random$patnr[,1],xlab="patient specific slope",main="")
qqnorm(lmeAge3$coefficients$random$patnr[,2])
qqline(lmeAge3$coefficients$random$patnr[,2])


### SEX ###
lmeSex1 = lme(diameter~sexe, random= (~1| patnr), data=marfan, , method="REML")
lmeSex2 = lme(diameter~sexe, random = (~0+sexe|patnr), data=marfan, , method="REML")
lmeSex3 = lme(diameter~sexe, random = (~1+sexe|patnr), data=marfan, , method="REML")
anova(lmeSex1, lmeSex2, lmeSex3)

summary(lmeSex)

par(mfrow=c(2,2))
plot(marfan$sexe, resid(lmeSex1))
abline(h=0, lty=2)
lines(unique(marfan$sexe), tapply(resid(lmeSex1), marfan$sexe, mean), col=2,lwd=3)
tapply(resid(lmeSex1), marfan$sexe, function(x) lines(0:9, x, lty=2))  # ----> x and y lengths differ
qqnorm(resid(lmeSex1))
qqline(resid(lmeSex1))
hist(lmeSex1$coefficients$random$patnr,xlab="patient specific intercept",main="")
qqnorm(lmeSex1$coefficients$random$patnr)
qqline(lmeSex1$coefficients$random$patnr)


####################
# PREDICTION MODEL #
####################

# ----> Prediction model voor leeftijd tussen 20 en 40, gebaseerd op code uit ppt
# lichte groei van diameter is te zien naarmate leeftijd stijgt

newdata <- data.frame(age = 20:40, patnr = 1:21)
designmatix <- model.matrix(~age, newdata)
predvar <- diag(designmatix %*% vcov(lmeAge3) %*% t(designmatix))
newdata$SE <- sqrt(predvar)
newdata$mean <- predict(lmeAge3, newdata = newdata,
                        level = 0)

plot(20:40, newdata$mean, xlab="Age", ylab="Diameter", type="l", col="red",
     ylim=c(0, 80), xlim=c(20,40))
lines(20:40, newdata$mean +1.96*newdata$SE, lty = 2, col="red")
lines(20:40, newdata$mean -1.96*newdata$SE, lty = 2, col="red")
for(i in 1:length(levels(marfan$patnr))){
  subject_id = levels(marfan$patnr)[i]
  lines(marfan$age[marfan$patnr == subject_id], 
        marfan$diameter[marfan$patnr == subject_id],
        col=alpha(i, 0.2))
}

####################
# Figure #
####################
a <- data.frame(..1 = character(), patnr = character(), metingnr = character(), age = character(), sexe = character(), diameter = character())
for(i in unique(marfan$patnr)){
 a[nrow(a) + 1,] = subset(subset(marfan, patnr == i), 
                   metingnr == max(metingnr, na.rm = TRUE))
}
pred_one <- data.frame(patnr = a$patnr, age = a$age, sexe = a$sexe, diameter = a$diameter)
#pred_one$metingnr <- as.numeric(pred_one$metingnr) + 1
pred_one$age <- as.numeric(pred_one$age) + 1
pred_one$patnr <- as.numeric(pred_one$patnr)


lmeAge2_1 = lme(diameter~age, random = (~1| patnr), data=marfan, method="REML")
lmeAge2_2 = lme(diameter~age, random = (~0+age|patnr), data=marfan , method="REML")
lmeAge2_3 = lme(diameter~age, random = (~1+age|patnr), data=marfan , method="REML")



designmatixfig <- model.matrix(~age, pred_one)
predvarfig <- diag(designmatixfig %*% vcov(lmeAge3) %*% t(designmatixfig))
pred_one$SE <- sqrt(predvarfig)
pred_one$mean <- predict(lmeAge3, newdata = pred_one,
                        level = 0)

plot(pred_one$age, pred_one$mean, xlab="Age", ylab="Diameter", type="l", col="red",
     ylim=c(30, 65), xlim=c(20,40))
lines(pred_one$age, pred_one$mean +1.96*pred_one$SE, lty = 2, col="red")
lines(pred_one$age, pred_one$mean -1.96*pred_one$SE, lty = 2, col="red")

for(i in 1:length(levels(marfan$patnr))){
  subject_id_fig = levels(marfan$patnr)[i]
  lines(marfan$age[marfan$patnr == subject_id_fig], 
        marfan$diameter[marfan$patnr == subject_id_fig],
        col=alpha(i, 0.2))
}



for(i in unique(marfan$patnr)) {
  subsetMarfan <- subset(marfan, patnr == i)
  metingen = length(subsetMarfan$metingnr)
  print(metingen)
  tapply(resid(lmeSex1),marfan$sexe, function(x) lines(0:metingen, 0:metingen, lty=2))  # ----> x and y lengths differ
}

#######NEW#########
library("lme4")
library('merTools')
library("scales")
library("ggplot2")
library("tidyverse")
######## Deze regels zijn miss niet nodig, maar in feite heb ik dit gedaan zodat alleen  de laatste meting in pred_one komt, maar niet nodig
a <- data.frame(..1 = character(), patnr = character(), metingnr = character(), age = character(), sexe = character(), diameter = character())
for(i in unique(marfan$patnr)){
  a[nrow(a) + 1,] = subset(subset(marfan, patnr == i), 
                           metingnr == max(metingnr, na.rm = TRUE))
}
pred_one <- data.frame(patnr = a$patnr, age = a$age, sexe = a$sexe, diameter = a$diameter)
#pred_one$metingnr <- as.numeric(pred_one$metingnr) + 1
pred_one$age <- as.numeric(pred_one$age) + 1
pred_one$patnr <- as.numeric(pred_one$patnr)


newdata <- data.frame(patnr = pred_one$patnr,age = pred_one$age,  sexe = pred_one$sexe, diameter = pred_one$diameter)
X<-split(newdata, newdata$sexe)

male <- data.frame(X["1"])
female <- data.frame(X["0"])

names(male)[names(male) == 'X1.age'] <- "age"
names(male)[names(male) == 'X1.patnr'] <- "patnr"
names(male)[names(male) == 'X1.sexe'] <- "sexe"
names(male)[names(male) == 'X1.diameter'] <- "diameter"

names(female)[names(female) == 'X0.age'] <- "age"
names(female)[names(female) == 'X0.patnr'] <- "patnr"
names(female)[names(female) == 'X0.sexe'] <- "sexe"
names(female)[names(female) == 'X0.diameter'] <- "diameter"


#####new data hoort je voorspellingen te zijn
#####Model 3 lmer
#####ci_prediction je interval en predictie
newdata <- data.frame(patnr = pred_one$patnr,age = pred_one$age,  sexe = pred_one$sexe, diameter = pred_one$diameter)
model3 <- lmer(diameter~age +(1+age|patnr), data=marfan, REML = FALSE)
ci_prediction <- predictInterval(model3, newdata=newdata)

ci_prediction_male <- predictInterval(model3, newdata=male)
ci_prediction_female <- predictInterval(model3, newdata=female)





par(mfrow=c(1,1))
for(i in unique(marfan$patnr)){
  subject_id = i ### Ik heb hier alle patnr, maar je kan in de for loop dus 4 zetten voor range
  print(subject_id)
  plot(marfan$age[i], marfan$diameter[i], 
       xlab="Age", ylab="Diameter", type="n", 
       col="red",
       ylim=c(18, 65), xlim=c(0, 70),
       main = paste("aaa: ", subject_id))
  lines(marfan$age[i], 
        marfan$diameter[i],
        col=alpha(i, 0.5))
  lines(1, ci_prediction[newdata$patnr == subject_id,1], lty = 1, col = i) #####Deze 3 zijn van ci_predictie, dus je upper- lower- en voorspelling
  lines(1, ci_prediction[newdata$patnr == subject_id,2], lty = 2, col= i)
  lines(1, ci_prediction[newdata$patnr == subject_id,3], lty = 2, col= i)
}










#par(mfrow=c(2,2))
for(i in female$patnr){
  subject_id = i
  print(subject_id)
  plot(female$age[i], female$diameter[i], 
       xlab="Age", ylab="Diameter", type="n", 
       col="red",
       ylim=c(18, 65), xlim=c(0, 140),
       main = paste("aaa: ", subject_id))
  lines(female$age[i], 
        female$diameter[i],
        col=alpha(i, 0.5))
  lines(female$age, ci_prediction_female[female$patnr == subject_id,1], lty = 1, col = i)
  lines(female$age, ci_prediction_female[female$patnr == subject_id,2], lty = 2, col= i)
  lines(female$age, ci_prediction_female[female$patnr == subject_id,3], lty = 2, col= i)
}







