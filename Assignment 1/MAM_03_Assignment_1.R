data <- read_csv("costefficacydata.csv")

#####SURVIVAL RATES#####
table(data$trt)
surv_rate_1 <- length(which(data$trt == 1 & data$event == 1)) / (length(which(data$trt == 1))) # Survival rate group 1
surv_rate_2 <- length(which(data$trt == 2 & data$event == 1)) / (length(which(data$trt == 2))) # Survival rate group 2

p1 <- mean(data$event[data$trt==1])
varp1 <- p1*(1-p1)/n_1

p2 <- mean(data$event[data$trt==2])
varp2 <- p2*(1-p2)/n_2

e <- p1-p2
vare <- varp1 + varp2


n_1 <- length(which(data$trt == 1)) 
n_2 <- length(which(data$trt == 2)) 

diff_surv_rate <- surv_rate_1 - surv_rate_2 # Difference survival rate

log_RR <- log(surv_rate_1)-log(surv_rate_2) # LogRR

var_logRR <- (1 - surv_rate_1) / (n_1 * surv_rate_1) + (1 - surv_rate_2) / (n_2 * surv_rate_2)
root <- sqrt(var_logRR) #SE

log_RR_up <- log_RR + 1.96 * root #95% CI logRR
log_RR_low <- log_RR - 1.96 * root

log_CI_up <- exp(log_RR_up) #95 CI RR
log_CI_low <- exp(log_RR_low)

var_surv1 <- 1 / (206 * surv_rate_1 * (1 - surv_rate_1))

var_surv2 <- 1 / (203 * surv_rate_2 * (1 - surv_rate_2))

var_diff <- 1/52 + 1/154 + 1/74 + 1/129


#####MEAN COSTS######
aggregate(data$costs, list(data$trt), FUN=mean) # Calculating the mean costs of both treatment groups

X<-split(data, data$trt)
tr_1 <- do.call(rbind.data.frame, X[1])
tr_2 <- do.call(rbind.data.frame, X[2])

t.test(tr_1$costs)$"conf.int"
t.test(tr_2$costs)$"conf.int"

t.test(data$costs)$"conf.int"

mean_1 <- mean(tr_1$costs)
mean_2 <- mean(tr_2$costs)

stde_costs_1 <- std.error(tr_1$costs)
stde_costs_2 <- std.error(tr_2$costs)

costs_1 <- c(tr_1$costs)
costs_2 <- c(tr_2$costs)


s_var_1 <- var(data$costs[data$trt==1])/n_1
s_var_2 <- var(data$costs[data$trt==2])/n_2

mean_diff <- mean_1 - mean_2

pooled_variance <- (((n_1-1)*s_var_1 + (n_2-1)*s_var_2)) / (n_1+n_2-2)

CI_costs_up = mean_diff + 2.06 * sqrt((pooled_variance/n_1)+(pooled_variance/n_2))
CI_costs_low = mean_diff - 2.06 * sqrt((pooled_variance/n_1)+(pooled_variance/n_2))


var_costs <- s_var_1 - s_var_2



#####ICER
ICER <- (mean_diff) / (diff_surv_rate)



varlogCE <- var_costs/(mean_diff*mean_diff) + vare/(e*e) - 2*covCE/(mean_diff*e)
covCE <- cov(data$event[data$trt==1], data$costs[data$trt==1])/n_1 + cov(data$event[data$trt==2], data$costs[data$trt==2])/n_2



log(mean_diff/e)
log(mean_diff/e) - 1.96*sqrt(varlogCE)
log(mean_diff/e) + 1.96*sqrt(varlogCE)



(diff_surv_rate * mean_diff) - ((diff_surv_rate * mean_diff - covCE^2)^2 - (diff_surv_rate^2 - var_diff)*(mean_diff - var_costs))^0.5/(diff_surv_rate^2 - var_diff)

(diff_surv_rate * mean_diff) + ((diff_surv_rate * mean_diff - covCE^2)^2 - (diff_surv_rate^2 - var_diff)*(mean_diff - var_costs))^0.5/(diff_surv_rate^2 - var_diff)



library("ICEinfer")
library("dampack")

hund_icers <- calculate_icers(data$costs, data$event, data$trt)
