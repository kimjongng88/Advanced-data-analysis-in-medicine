library(haven)
data <- read_sav("renaltx.sav")

table(data$dgf)
data$dgf = as.factor(as.character(data$dgf))
data$uprotein = as.factor(as.character(data$uprotein))


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

table1(~dgf + acclft + aantalre + creat + predias + prac + uprotein + cregsh |gstatus,data = data, overall = F, extra.col = list('P-value' = pvalue))

########  UNIVARIABLE REGRESSION ########
acclft <- summary(glm(gstatus ~ acclft, family = "binomial", data = data, na.action = na.omit))
dgf <- summary(glm(gstatus ~ dgf, family = "binomial", data = data, na.action = na.omit))
aantalre <- summary(glm(gstatus ~ aantalre, family = "binomial", data = data, na.action = na.omit))
creat <- summary(glm(gstatus ~ creat, family = "binomial", data = data, na.action = na.omit))
predias <- summary(glm(gstatus ~ predias, family = "binomial", data = data, na.action = na.omit))
prac <- summary(glm(gstatus ~ prac, family = "binomial", data = data, na.action = na.omit))
uprotein <- summary(glm(gstatus ~ uprotein, family = "binomial", data = data, na.action = na.omit))
cregsh <- summary(glm(gstatus ~ cregsh, family = "binomial", data = data, na.action = na.omit))
gsurv <- summary(glm(gstatus ~ gsurv, family = "binomial", data = data, na.action = na.omit))

table1(~dgf + acclft + aantalre + creat + predias + prac + uprotein + cregsh|gsurv,data = data, overall = F, extra.col = list('P-value' = pvalue))

########  UNIVARIABLE REGRESSION ########
acclft_surv <- summary(lm(gsurv ~ acclft, data = data, na.action = na.omit))
dgf_surv <- summary(lm(gsurv ~ dgf, data = data, na.action = na.omit))
aantalre_surv <- summary(lm(gsurv ~ aantalre, data = data, na.action = na.omit))
creat_surv <- summary(lm(gsurv ~ creat,data = data, na.action = na.omit))
predias_surv <- summary(lm(gsurv ~ predias, data = data, na.action = na.omit))
prac_surv <- summary(lm(gsurv ~ prac,  data = data, na.action = na.omit))
uprotein_surv <- summary(lm(gsurv ~ uprotein, data = data, na.action = na.omit))
cregsh_surv <- summary(lm(gsurv ~ cregsh, data = data, na.action = na.omit))
gstatus_surv <- summary(lm(gsurv ~ gstatus, data = data, na.action = na.omit))

library(MASS)
library(rms)
model.correct <- lrm(gstatus~dgf + acclft + aantalre + creat + predias + prac + uprotein + cregsh,data = data, x=TRUE, y=TRUE)
summary(model.correct)
fastbw(model.correct)
plot(calibrate(model.correct,B=1000))
validate(model.correct, method = "boot", B=1000)


multiStatus <- glm(gstatus ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=data)
summary(multiStatus) # show results

multiSurv <- glm(gsurv ~ acclft + dgf + aantalre + creat + predias + prac + uprotein + cregsh, data=data)
summary(multiSurv) # show results


fastbw(multiStatus )
plot(calibrate(model.correct,B=1000))
validate(model.correct, method = "boot", B=1000)
