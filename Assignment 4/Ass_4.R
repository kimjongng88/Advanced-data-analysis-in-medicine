library("mice")
imputed=mice(dlong)
data <- complete(imputed)

imputedd=mice(d)
dd<- complete(imputedd)


#install.packages("table1")
library("table1")

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

###### Multi variable
table1(~years + gfr + creat + map + sex_pat + age_at_tx + bmi + type_dia + duur_dia + sexdon + agedon +retrans + stat_gra|stat_pat,data = data, overall = F, extra.col = list('P-value' = pvalue))

#####Verdeling
library('dplyr')
library('funModeling')
data %>%
  dplyr::summarise(meanLE=mean(time_to_graft_failure,na.rm=TRUE),
                   medLE=median(time_to_graft_failure,na.rm=TRUE),
                   sd=sd(time_to_graft_failure,na.rm=TRUE),
                   iqr=IQR(time_to_graft_failure,na.rm=TRUE),
                   Q1=quantile(time_to_graft_failure,probs=0.25,na.rm=TRUE),
                   Q3=quantile(time_to_graft_failure,probs=0.75),
                   n=n())

funModeling::profiling_num(data)

funModeling::plot_num(data)

library("skimr")
varlist <- c("n_missing","complete_rate")
data %>% 
  select(time_to_graft_failure) %>% 
  skimr::skim_without_charts() %>%
  skimr::yank("numeric") %>%
  dplyr::select(-one_of(varlist))



########NEW########
#install.packages("cmprsk")
#install.packages("JMbayes2")
library("JMbayes2")
library("cmprsk")
library("survival")
library("survminer")
library("nlme")
library("survival")
library("JM")   
library("JMbayes")

#############################################
#############################################
####################COX######################
#############################################
#############################################

# COX model for death
m1=coxph(Surv(time_to_death,as.numeric(stat_pat)-1)~age_at_tx+agedon+type_dia+duur_dia,data=d,x=TRUE)
summary(m1)  
qqnorm(resid(m1))

ggsurvplot(survfit(m1), data = data, color = "#2E9FDF",
           ggtheme = theme_minimal())


# COX model for graft_status
m2=coxph(Surv(time_to_graft_failure,as.numeric(stat_gra)-1)~age_at_tx+agedon+type_dia+duur_dia,data=d,x=TRUE)
summary(m2)  
qqnorm(resid(m2))


ggsurvplot(survfit(m2), data = data, color = "#2E9FDF",
           ggtheme = theme_minimal())



#Combine
d$dfs = d$time_to_graft_failure
d$dfs[d$stat_gra==0] = d$ time_to_death[d$stat_gra==0]
d$dfs_stat = d$stat_gra
d$dfs_stat[d$stat_pat==1] = 1
m_com = coxph(Surv(dfs, as.numeric(dfs_stat))~age_at_tx+agedon+type_dia+duur_dia,data=d,x=TRUE)

ggsurvplot(survfit(m_com), data = data, color = "#2E9FDF",
           ggtheme = theme_minimal())



cc1=d[,which(names(d) %in% c("sex_pat","agedon","screat1"))]
cc1$sex_pat = as.numeric(cc1$sex_pat)
d$status=0
d$status[as.integer(as.numeric(d$stat_gra))==1]=1
d$status[as.integer(as.numeric(d$stat_pat))==1]=2

plot(cuminc(d$dfs, d$status, cencode=0))
crr(d$dfs, d$status, failcode=1, cencode=0, cov1=as.matrix(cc1))



#############################################
#############################################
####################LME######################
#############################################
#############################################
# LME for gfr
dlong1=subset(dlong,is.na(gfr)==FALSE)      # verwijder records met missing gfr-metingen
dlong1$gfr[dlong1$gfr > 200]=NA             # verwijder 2 records met onwaarschijnlijke (onmogelijk) hoge gfr metingen
dlong1=subset(dlong1,is.na(gfr)==FALSE)
p3a=lme(gfr~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong1,method="REML")
summary(p3a)   
pred_p3a = predict(p3a, level = 1)

plot(dlong1$years[dlong1$ID==2],dlong1$gfr[dlong1$ID==2],
     type="s",lty=2,ylim=c(0,100),
     xlab="years",ylab="GFR",main="patient 2")
lines(dlong1$years[dlong1$ID==2], pred_p3a[dlong1$ID==2],
      lty=1,col=1,lwd=2)

plot(p3a)
qqnorm(p3a)
qqnorm(resid(p3a))



# LME for creat
dlong2=subset(dlong,is.na(creat)==FALSE)      # verwijder records met missing creat-metingen
dlong2$creat[dlong1$creat > 1000]=NA             # Outliers van Creat
dlong2=subset(dlong2,is.na(creat)==FALSE)
p3a_cr=lme(creat~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong2,method="REML")
summary(p3a_cr)   


pred_p3a_cr = predict(p3a_cr, level = 1)

plot(dlong1$years[dlong2$ID==2],dlong2$creat[dlong2$ID==2],
     type="s",lty=2,ylim=c(0,100),
     xlab="years",ylab="creat",main="patient 2")
lines(dlong1$years[dlong2$ID==2], pred_p3a_cr[dlong2$ID==2],
      lty=1,col=1,lwd=2)

plot(p3a_cr)
qqnorm(p3a_cr)
qqnorm(resid(p3a_cr))




# LME for map
dlong3=subset(dlong,is.na(map)==FALSE)      # verwijder records met missing gfr-metingen, geen outliers
dlong3=subset(dlong3,is.na(map)==FALSE)
p3a_map=lme(map~ns(years,df=2),random=~1|ID,data=dlong3,method="REML")
summary(p3a_map)   



pred_p3a_map = predict(p3a_map, level = 1)

plot(dlong3$years[dlong3$ID==2],dlong3$creat[dlong3$ID==2],
     type="s",lty=2,ylim=c(0,130),
     xlab="years",ylab="MAP",main="patient 2")
lines(dlong3$years[dlong3$ID==2], pred_p3a_map[dlong3$ID==2],
      lty=1,col=1,lwd=2)


plot(p3a_map)
qqnorm(p3a_map)
qqnorm(resid(p3a_map))

# you can see 838 patients with both analyses: that is an absolute must
# now the JOINT model
j1=jointModel(p3a,m1,timeVar="years")
summary(j1)
j2 = jointModel(p3a_cr,m1, timeVar="years")
summary(j2)
j3 = jointModel(p3a_map,m1, timeVar="years")
summary(j3)



######### Graft
j1g=jointModel(p3a,m2,timeVar="years")
summary(j1g)
j2g = jointModel(p3a_cr,m2, timeVar="years")
summary(j2g)
j3g = jointModel(p3a_map,m2, timeVar="years")
summary(j3g)


dlong1=subset(dlong,is.na(gfr)==FALSE & is.na(creat)==FALSE & is.na(map)==FALSE)  # verwijder alle records met missing
dlong1=subset(dlong1,(dlong1$years>=dlong1$time_to_death)==FALSE)
dlong1$gfr[dlong1$gfr > 200]=NA    # verwijder die 2 rare gfr metingen
dlong1=subset(dlong1,is.na(gfr)==FALSE)
dim(dlong1)   # check hoeveel records er in dlong1 zitten

lmemodel_gfr = lme(gfr ~ ns(years,2), random = ~1+ns(years,2)|ID, data=dlong1)  
summary(lmemodel_gfr)    # check hoeveel pats en hoeveel records in de analyse zitten (dit moet hetzelfde aantal als bij de 2 andere lme's zijn)
lmemodel_creat = lme(creat ~ years, random = ~1+years|ID, data=dlong1)
lmemodel_map = lme(map ~ years, random = ~1+years|ID, data=dlong1)
jmresult = jm(m1, list(lmemodel_gfr, lmemodel_creat, lmemodel_map), time_var="years")
summary(jmresult)

jmresult2 = jm(m2, list(lmemodel_gfr, lmemodel_creat, lmemodel_map), time_var="years")
summary(jmresult2)


#### https://www.researchgate.net/publication/47457473_JM_An_R_Package_for_the_Joint_Modelling_of_Longitudinal_and_Time-to-Event_Data




# We are interested in producing predictions of survival probabilities for Patient 155
dataP1 <- dlong1[dlong1$ID == 1, ]
len_id <- nrow(dataP1)


sfit3 <- survfitJM(j1, newdata = dataP1,idVar = "ID") 
sfit4 <- survfitJM(j2, newdata = dataP1,idVar = "ID") 
sfit5 <- survfitJM(j3, newdata = dataP1,idVar = "ID") 


par(mfrow=c(1,2))
plotfit3 <- plot(sfit3, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 1")
plotfit4 <- plot(sfit4, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 1")
plotfit5 <- plot(sfit5, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 1")



dataP1 <- dlong1[dlong1$ID == 1, ]
len_id <- nrow(dataP1)
sfit3g <- survfitJM(j1g, newdata = dataP1,idVar = "ID") 
sfit4g <- survfitJM(j2g, newdata = dataP1,idVar = "ID") 
sfit5g <- survfitJM(j3g, newdata = dataP1,idVar = "ID") 


par(mfrow=c(1,2))
plotfit3 <- plot(sfit3g, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 1")
plotfit4 <- plot(sfit4g, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 1")
plotfit5 <- plot(sfit5g, estimator="mean", include.y = TRUE, conf.int=0.95, fill.area=TRUE, col.area="lightblue", main="Patient 1")

#####test
patients_to_predict <- seq(1,100)
newdata <- dlong[dlong$ID %in% patients_to_predict, ]
newdata <- merge(newdata, d[d$ID%in% patients_to_predict,], by=c("ID", "years","sex_pat", "age_at_tx", "bmi", "type_dia", "duur_dia", "sexdon", "agedon", "retrans", "stat_gra", "stat_pat", "time_to_graft_failure", "time_to_death"))


predLong1 <- predict(jmresult, newdata = newdata, return_newdata = TRUE)
predSurv <- predict(jmresult, newdata = newdata, process = "event",times = seq(5, 12, length.out = 51),return_newdata = TRUE)
plot(predlong1)

newdata$stat_pat <- as.numeric(newdata$stat_pat) - 1
newdata$time_to_death <- as.numeric(newdata$time_to_death)


newdata$event <- as.numeric(newdata$stat_pat != "levend")
roc <- tvROC(jmresult, newdata = newdata, Tstart = 1, Dt = 3, idVar = "ID")



dlong2 <- dlong1
dlong2$stat_pat <- as.factor(dlong2$stat_pat)
dlong2$event <- as.numeric(dlong2$stat_pat != "levend")
roc <- tvROC(jmresult, newdata = dlong2, Tstart = 8, Dt = 3)



auc1 <- aucJM(j1, dlong1, Tstart = 5, Dt = 3, idVar = "ID")
auc2 <- aucJM(j2, dlong2, Tstart = 5, Dt = 3, idVar = "ID")
auc3 <- aucJM(j3, dlong3, Tstart = 5, Dt = 3, idVar = "ID")

auc1g <- aucJM(j1g, dlong1, Tstart = 5, Dt = 3, idVar = "ID")
auc2g <- aucJM(j2g, dlong2, Tstart = 5, Dt = 3, idVar = "ID")
auc3g <- aucJM(j3g, dlong3, Tstart = 5, Dt = 3, idVar = "ID")

roc <- rocJM(j1, dt = c(2, 4, 8), dlong1, idVar = "ID")
roc2 <- rocJM(j2, dt = c(2, 4, 8), dlong2, idVar = "ID")
roc3 <- rocJM(j3, dt = c(2, 4, 8), dlong3, idVar = "ID")


