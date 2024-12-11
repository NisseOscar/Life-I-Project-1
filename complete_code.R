
############# MT7056 Assignment 1 ###############################

########## TASK 1 ###############################

# Import base packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(survival)
library(mvna)
library(survminer)
library(mice)
library(EnvStats)

# Parameters
a = 5.5
b = 22.5
t <- 1:350/10
plot_width = 8
plot_height = 10
y_lim = c(-0.5, 0.5)
# options(repr.plot.width = plot_width, repr.plot.height = plot_height)

set.seed(123)

# Set potential states for model (our case dead or alive)
trans_matrix <- matrix(nrow = 2, ncol = 2, FALSE)
trans_matrix[1, 2] <- TRUE
trans_matrix


# plot theoretical survival function
p <- ggplot(df, aes(x = t, y = exp(-(t / b)^a)), size = 1) +
    geom_line(show.legend = FALSE) +
    xlab("t") +
    ylab("S(t)") +
    theme_bw()
ggsave(paste("plots-relative/task1/theoretical_survival.png"), plot = p, width = plot_width, height = plot_height, units = "in", dpi = 600, device = "png")

# Plot theoretical Hazard function
p <- ggplot(df, aes(x = t, y = (t / b)^a), size = 1) +
    geom_line(show.legend = FALSE) +
    xlab("t") +
    ylab("S(t)") +
    theme_bw()
ggsave(paste("plots-relative/task1/theoretical_hazard.png"), plot = p, width = plot_width, height = plot_height, units = "in", dpi = 600, device = "png")


################### Wihtout Cencoring  ############################"
# Create side by side ggplt plot for n=10,100,500,1000 and save in plots folder
n_values = c(10,100,500,1000)
for (n in n_values){
    T = rweibull(n, a, b)
    dN = rep(1, n)
    t_tmp = t[t <= max(T)]
    theoretical <- exp(-(t_tmp / b)^a)

    df_tmp <- data.frame(id=1:n,from=0,to=1,time=T)
    na_fit = mvna(df_tmp, state.names = c(0, 1), tra = trans_matrix, cens.name = "NULL")
    na_pred = predict(na_fit, times = t_tmp, var.type = "aalen", level = 0.95)$"0 1"
    S_na_values <- exp(-na_pred$na) - theoretical
    conf_na_low = exp(-na_pred$lower) -theoretical
    conf_na_high = exp(-na_pred$upper) - theoretical

    # Calculate Kepler-Meier estimators
    kp_fit = survfit(Surv(T, dN) ~ 1,type="kaplan-meier")
    kp_fit = summary(kp_fit, times = t_tmp)

    # Create dataset
    df <- rbind(
        data.frame(t=t_tmp, S=S_na_values, Method="Nelson-Aalen", c_low = conf_na_low,c_high=conf_na_high, Confidence_Interval = "95% Nelson-Aalen"), 
        data.frame(t = t_tmp, S = kp_fit$surv-theoretical, Method = "Kepler-Meier", c_low = kp_fit$lower-theoretical,c_high= kp_fit$upper-theoretical, Confidence_Interval = "95% Kepler-Meier"),
        data.frame(t = t_tmp, S = 0, Method = "Theoretical", c_low = 0, c_high = 0, Confidence_Interval = "")
    )

    # Plot and save results
    p <- ggplot(df,aes(x=t, y=S,color=Method),size=1) + 
        geom_line(show.legend = FALSE) +
        geom_ribbon(aes(ymin=c_low, ymax=c_high, fill=Confidence_Interval), show.legend = FALSE, alpha = 0.2,linetype=2) +
        xlab("t") + ylab("S(t) minus theoretical true value") + theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5))
    ggsave(paste("plots-relative/task1/survival_simulation_n", n, ".png"), plot = p, width = plot_width, height = plot_height, units = "cm", device = "png")
}

# Create intervals for cumulative hazard function
n_values <- c(10, 100, 500, 1000)
for (n in n_values) {
    T <- rweibull(n, a, b)
    dN <- rep(1, n)
    t_tmp = t[t <= max(T)]
    theoretical <- (t_tmp / b)^a


    # Calculate Nelson-Aalen  estimators
    df_tmp <- data.frame(id = 1:n, from = 0, to = 1, time = T)
    na_fit = mvna(df_tmp, state.names = c(0, 1), tra = trans_matrix, cens.name = "NULL")
    na_pred <- predict(na_fit, times = t_tmp, var.type = "aalen", level = 0.95)$"0 1"

    # Calculate Kepler-Meier estimators
    kp_fit <- survfit(Surv(T, dN) ~ 1, type = "kaplan-meier")
    kp_fit = summary(kp_fit, times = t_tmp)
    A_kp_values <- -log(kp_fit$surv) - theoretical
    conf_kp_low <- -log(kp_fit$lower) - theoretical
    conf_kp_high <- -log(kp_fit$upper) - theoretical
    
    # Create dataset
    df <- rbind(
        data.frame(t = t_tmp, A = na_pred$na-theoretical, Method = "Nelson-Aalen", c_low = na_pred$lower-theoretical, c_high = na_pred$upper-theoretical, Confidence_Interval = "95% Nelson-Aalen"),
        data.frame(t = t_tmp, A = A_kp_values, Method = "Kepler-Meier", c_low = conf_kp_low, c_high = conf_kp_high, Confidence_Interval = "95% Kepler-Meier"),
        data.frame(t = t_tmp, A = 0, Method = "Theoretical", c_low = 0, c_high = 0, Confidence_Interval = "")
    )

    # Plot and save results
    p <- ggplot(df, aes(x = t, y = A, color = Method), size = 1) +
        geom_line(show.legend = FALSE) +
        geom_ribbon(aes(ymin = c_low, ymax = c_high, fill = Confidence_Interval), alpha = 0.2,show.legend = FALSE,,linetype=2) +
        xlab("t") +
        ylab("A(t) minus theoretical true value") +
        theme_bw() + coord_cartesian(ylim = c(-1, 3))
    ggsave(paste("plots-relative/task1/hazard_simulation_n", n, ".png"), plot = p, width = plot_width, height = plot_height, units = "cm", device = "png")
}

################### With Cencoring  ############################"
# Create side by side ggplt plot for n=10,100,500,1000 and save in plots folder
n_values <- c(10, 100, 500, 1000)
for (n in n_values) {
    T <- rweibull(n, a, b)
    # Censor observations uniformly on the interval 20-60
    C <- runif(n, 20, 60)
    dN <- sapply(1:n, function(i) ifelse(T[i] <= C[i], 1, 0))
    t_tmp

    # Theoretical true
    theoretical <- exp(-(t_tmp / b)^a)

    # Calculate the weibull mle estimates ignoring all censored data
    subset = which(dN == 1)
    mle_fit <- eweibull(T[subset])
    mle_fit$parameters
    mle_est_theoretical <- exp(-(t_tmp / mle_fit$parameters[2])^mle_fit$parameters[1]) - theoretical

    # Calculate Nelson-Aalen  estimators
    df_tmp <- data.frame(id = 1:n, from = 0, to = 1,time=T, status = dN)
    na_fit <- mvna(df_tmp, state.names = c(0, 1), tra = trans_matrix, cens.name = "status")
    na_pred <- predict(na_fit, times = t_tmp, var.type = "aalen", level = 0.95)$"0 1"
    S_na_values <- exp(-na_pred$na) - theoretical
    conf_na_low <- exp(-na_pred$lower) - theoretical
    conf_na_high <- exp(-na_pred$upper) - theoretical

    # Calculate Kepler-Meier estimators
    kp_fit <- survfit(Surv(T, dN) ~ 1, type = "kaplan-meier")
    kp_fit <- summary(kp_fit, times = t_tmp)

    # Create dataset
    df <- rbind(
        data.frame(t=t_tmp , S=S_na_values, Method="Nelson-Aalen", c_low = conf_na_low,c_high=conf_na_high, Confidence_Interval = "95% Nelson-Aalen"), 
        data.frame(t = t_tmp, S = kp_fit$surv-theoretical, Method = "Kepler-Meier", c_low = kp_fit$lower-theoretical,c_high= kp_fit$upper-theoretical, Confidence_Interval = "95% Kepler-Meier"),
        data.frame(t = t_tmp, S = 0, Method = "Weibull", c_low = 0, c_high = 0, Confidence_Interval = ""),
        data.frame(t = t_tmp, S = mle_est_theoretical, Method = "MLE Weibull", c_low = mle_est_theoretical, c_high = mle_est_theoretical, Confidence_Interval = "")
    )

    # Plot and save results
    p <- ggplot(df, aes(x = t, y = S, color = Method), size = 1) +
        geom_line(show.legend = FALSE) +
        geom_ribbon(aes(ymin = c_low, ymax = c_high, fill = Confidence_Interval), alpha = 0.2,show.legend = FALSE,linetype=2) +
        xlab("t") +
        ylab("S(t) minus theoretical true value") +
        theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5))
    ggsave(paste("plots-relative/task1/censored_survival_simulation_n", n, ".png"), plot = p, width = plot_width, height = plot_height, units = "cm", device = "png")
}


# Create intervals for cumulative hazard function
n_values <- c(10, 100, 500, 1000)
for (n in n_values) {
    T <- rweibull(n, a, b)
    # Censor observations
    C <- runif(n, 20, 60)
    T <- pmin(T, C)
    dN <- rep(1, n)
    t_tmp = t[t <= max(T)]
    theoretical <- (t_tmp / b)^a

    # MLE estimation
    subset = which(dN == 1)
    mle_fit <- eweibull(T[subset])
    mle_est_theoretical <- (t_tmp / mle_fit$parameters[2])^mle_fit$parameters[1] - theoretical


    # Calculate Nelson-Aalen  estimators
    df_tmp <- data.frame(id = 1:n, from = 0, to = 1, time = T, status = dN)
    na_fit <- mvna(df_tmp, state.names = c(0, 1), tra = trans_matrix, cens.name = "status")
    na_pred <- predict(na_fit, times = t_tmp, var.type = "aalen", level = 0.95)$"0 1"

    # Calculate Kepler-Meier estimators
    kp_fit <- survfit(Surv(T, dN) ~ 1, type = "kaplan-meier")
    kp_fit <- summary(kp_fit, times = t_tmp)
    A_kp_values <- -log(kp_fit$surv) - theoretical
    conf_kp_low <- -log(kp_fit$lower) - theoretical
    conf_kp_high <- -log(kp_fit$upper) - theoretical

    # Create dataset
    df <- rbind(
        data.frame(t = t_tmp, A = na_pred$na-theoretical, Method = "Nelson-Aalen", c_low = na_pred$lower-theoretical, c_high = na_pred$upper-theoretical, Confidence_Interval = "95% Nelson-Aalen"),
        data.frame(t = t_tmp, A = A_kp_values, Method = "Kepler-Meier", c_low = conf_kp_low, c_high = conf_kp_high, Confidence_Interval = "95% Kepler-Meier"),
        data.frame(t = t_tmp, A = 0, Method = "Theoretical", c_low = 0, c_high = 0, Confidence_Interval = ""),
        data.frame(t = t_tmp, A = mle_est_theoretical, Method = "MLE Weibull", c_low = mle_est_theoretical, c_high = mle_est_theoretical, Confidence_Interval = "")
    )

    # Plot and save results
    p <- ggplot(df, aes(x = t, y = A, color = Method), size = 1) +
        geom_line(show.legend = FALSE) +
        geom_ribbon(aes(ymin = c_low, ymax = c_high, fill = Confidence_Interval), alpha = 0.2, ,linetype=2,show.legend = FALSE) +
        xlab("t") +
        ylab("A(t) minus theoretical true value") +
        theme_bw() + coord_cartesian(ylim = c(-1, 3))
    ggsave(paste("plots-relative/task1/censored_simulation_n", n, ".png"), plot = p, width = plot_width, height = plot_height, units = "cm", device = "png")
}


# Plot function to show the legend underneath graph
# order according to df set MLE weibull last
order <- c('Nelson-Aalen','Kepler-Meier', 'Theoretical')
p <- ggplot(df, aes(x = t, y = S, color = Method), size = 1) +
    geom_line(show.legend = TRUE) +
    geom_ribbon(aes(ymin = c_low, ymax = c_high, fill = Method,color=Method), alpha = 0.2, ,linetype=2,show.legend = TRUE) +
    xlab("t") +
    ylab("A(t) minus theoretical true value") +
    theme_bw() + coord_cartesian(ylim = c(-1, 3))+
    theme(legend.position = "bottom")
    # scale_fill_discrete(breaks = order, labels = order, name = "Method")
ggsave(paste("plots-relative/task1/legend_plot.png"), plot = p, width = plot_width, height = plot_height, units = "cm", device = "png")
p

##### Estimate Means ###################
set.seed(123)

# Generate random Weibull distributed numbers
n10 <- rweibull(10, shape = 5.5, scale = 22.5)
n100 <- rweibull(100, shape = 5.5, scale = 22.5)
n500 <- rweibull(500, shape = 5.5, scale = 22.5)
n1000 <- rweibull(1000, shape = 5.5, scale = 22.5) 

# Kaplan-Meier estimate
KM10 <- survfit(Surv(n10, rep(c(1), each = 10))~1, stype = 1)
KM100 <- survfit(Surv(n100, rep(c(1), each = 100))~1, stype = 1)
KM500 <- survfit(Surv(n500, rep(c(1), each = 500))~1, stype = 1)
KM1000 <- survfit(Surv(n1000, rep(c(1), each = 1000))~1, stype = 1)

medianKM = matrix(c(as.numeric(quantile(KM10, probs = c(0.5))), as.numeric(quantile(KM100, probs = c(0.5))), as.numeric(quantile(KM500, probs = c(0.5))), as.numeric(quantile(KM1000, probs = c(0.5)))), ncol=3, byrow=TRUE)
rownames(medianKM) <- c('n=10','n=100','n=500','n=1000')
colnames(medianKM) = c('Median','Lower CI','Upper CI')
medianKM_table <- as.table(medianKM)
medianKM_table

# Nelson-Aalen estimate
NA10 <- survfit(Surv(n10, rep(c(1), each = 10))~1, stype = 2, ctype = 1)
NA100 <- survfit(Surv(n100, rep(c(1), each = 100))~1, stype = 2, ctype = 1)
NA500 <- survfit(Surv(n500, rep(c(1), each = 500))~1, stype = 2, ctype = 1)
NA1000 <- survfit(Surv(n1000, rep(c(1), each = 1000))~1, stype = 2, ctype = 1)

medianNA = matrix(c(as.numeric(quantile(NA10, probs = c(0.5))), as.numeric(quantile(NA100, probs = c(0.5))), as.numeric(quantile(NA500, probs = c(0.5))), as.numeric(quantile(NA1000, probs = c(0.5)))), ncol=3, byrow=TRUE)
rownames(medianNA) <- c('n=10','n=100','n=500','n=1000')
colnames(medianNA) = c('Median','Lower CI','Upper CI')
medianNA_table <- as.table(medianNA)
medianNA_table


# PART 1, task 2
# Censoring times that are uniform on (20,60)
rand_cens10 <- runif(10, min=20, max=60)
rand_cens100 <- runif(100, min=20, max=60)
rand_cens500 <- runif(500, min=20, max=60)
rand_cens1000 <- runif(1000, min=20, max=60)

# Event
event10 <- ifelse(n10 <= rand_cens10, 1, 0)
event100 <- ifelse(n100 <= rand_cens100, 1, 0)
event500 <- ifelse(n500 <= rand_cens500, 1, 0)
event1000 <- ifelse(n1000 <= rand_cens1000, 1, 0)

# Observed values
n10_cens <- pmin(n10, rand_cens10)
n100_cens <- pmin(n100, rand_cens100)
n500_cens <- pmin(n500, rand_cens500)
n1000_cens <- pmin(n1000, rand_cens1000)

# Kaplan-Meier estimate
KM10_cens <- survfit(Surv(n10_cens, event10)~1, stype = 1)
KM100_cens <- survfit(Surv(n100_cens, event100)~1, stype = 1)
KM500_cens <- survfit(Surv(n500_cens, event500)~1, stype = 1)
KM1000_cens <- survfit(Surv(n1000_cens, event1000)~1, stype = 1)

medianKM_cens = matrix(c(as.numeric(quantile(KM10_cens, probs = c(0.5))), as.numeric(quantile(KM100_cens, probs = c(0.5))), as.numeric(quantile(KM500_cens, probs = c(0.5))), as.numeric(quantile(KM1000_cens, probs = c(0.5)))), ncol=3, byrow=TRUE)
rownames(medianKM_cens) <- c('n=10','n=100','n=500','n=1000')
colnames(medianKM_cens) = c('Median','Lower CI','Upper CI')
medianKM_table_cens <- as.table(medianKM_cens)
medianKM_table_cens

# Nelson-Aalen estimate
NA10_cens <- survfit(Surv(n10_cens, event10)~1, stype = 2, ctype = 1)
NA100_cens <- survfit(Surv(n100_cens, event100)~1, stype = 2, ctype = 1)
NA500_cens <- survfit(Surv(n500_cens, event500)~1, stype = 2, ctype = 1)
NA1000_cens <- survfit(Surv(n1000_cens, event1000)~1, stype = 2, ctype = 1)

medianNA_cens = matrix(c(as.numeric(quantile(NA10_cens, probs = c(0.5))), as.numeric(quantile(NA100_cens, probs = c(0.5))), as.numeric(quantile(NA500_cens, probs = c(0.5))), as.numeric(quantile(NA1000_cens, probs = c(0.5)))), ncol=3, byrow=TRUE)
rownames(medianNA_cens) <- c('n=10','n=100','n=500','n=1000')
colnames(medianNA_cens) = c('Median','Lower CI','Upper CI')
medianNA_table_cens <- as.table(medianNA_cens)
medianNA_table_cens

################# TASK 2 ###############################
T1 <- n500_cens
eventT1 <- event500 # Using censored data from part 1
T2_non_cens <- rweibull(100, shape = 4.5, scale = 28)
eventT2 <- ifelse(T2_non_cens <= rand_cens100, 1, 0)
T2 <- pmin(T2_non_cens, rand_cens100)

# Logrank test
combined_T <- c(T1, T2)
combined_eventT <- c(eventT1, eventT2)
group <- c(rep(c(0), each = 500), rep(c(1), each = 100))
logrank_res <- survdiff(Surv(combined_T, combined_eventT) ~ group, rho = 0) # rho=0 indicate that it is a Logrank test
print(logrank_res)

# Fit a Cox regression model
cox_model <- coxph(Surv(combined_T, combined_eventT) ~ group, method = "breslow")
summary(cox_model)

# Kaplan-Meier
KMpart2 <- survfit(Surv(combined_T, combined_eventT) ~ group, stype = 1)
# KM_plot <- plot(KMpart2)

# Plot difference between Kaplan-Meier curves and Cox regression curves
plot(KMpart2, col = c("blue", "red"), lty = 1:1, xlab = "Time", ylab = "Survival Probability")
cox_surv <- survfit(cox_model, newdata = data.frame(group = c(0, 1)))
lines(cox_surv, col = c("blue", "red"), lty = 2:2)
legend("topright",
    legend = c("Group 1 (KM)", "Group 2 (KM)", "Group 1 (Cox)", "Group 2 (Cox)"),
    col = c("blue", "red", "blue", "red"), lty = c(1, 1, 2, 2)
)
