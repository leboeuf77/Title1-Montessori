# LeBoeuf, Goldstein-Greenwood, & Lillard
# Racial discipline disproportionality in Montessori and non-Montessori schools
# This script include only the analysis using multilevel negative binomial models
# It uses data created/modified in the script "1_RD-RRR-logit analysis.R". 

# This code includes pseudo-random processes (bootstrapping). Values reported in the paper
# are based on running the code from top to bottom with the seed set below (091395).

# Load packages
invisible(lapply(list('dplyr', 'tidyr', 'stringi', 'ggeffects', 
                      'ggplot2', 'glmmTMB', 'DHARMa'),
                 function(pkg) library(pkg, character.only = TRUE)))

# Initialize function for calculating bias-corrected CIs
# See Bootstrap Methods (Chernick, 2008) Ch. 3 for info on bias-corrected confidence intervals
# Hat tip to the infer R package (part of tidymodels; made available under an MIT License) 
# for exemplifying a function for calculating bias-corrected confidence intervals
bc_ci <- function(x, level, point_estimate) {
  p <- mean(x <= point_estimate)
  z0 <- stats::qnorm(p)
  # z_alpha_2 is z_(alpha/2)
  z_alpha_2 <- stats::qnorm((1 + c(-level, level)) / 2)
  new_probs <- stats::pnorm(2 * z0 + z_alpha_2)
  # CIs
  ci_vec <- stats::quantile(x, probs = new_probs)
  ci_vec
}

# Read data and set seed for random-process reproducibility 
dat <- read.csv("DiscAnalysisDatNBmlm.csv")
set.seed(091395)

################################################################################################################################################
# Black/Hispanic/White_num = the number of students in each racial group in a school
# Need to add to get the total suspension counts for each racial group 

# Total Black OSS & ISS
dat$Black_count_oss <- (dat$only.one.out.of.school_Black + dat$more.than.one.out.of.school_Black)
dat$Black_count_iss <- dat$in.school.suspensions_Black

# Total White OSS & ISS 
dat$White_count_oss <- (dat$only.one.out.of.school_White + dat$more.than.one.out.of.school_White)
dat$White_count_iss <- dat$in.school.suspensions_White

# Total Hispanic OSS & ISS 
dat$Hispanic_count_oss <- (dat$only.one.out.of.school_Hispanic + dat$more.than.one.out.of.school_Hispanic)
dat$Hispanic_count_iss <- dat$in.school.suspensions_Hispanic

# Overall OSS & ISS count for all races 
dat$tot_OSSCount <- (dat$Black_count_oss + dat$White_count_oss + dat$Hispanic_count_oss)
dat$tot_ISSCount <- (dat$Black_count_iss + dat$White_count_iss + dat$Hispanic_count_iss)

########## Overall ISS and OSS counts Bootstrapping ##########################################################################################################

# OSS ##########
observedmodossAll <- glmmTMB(tot_OSSCount ~ Montessori + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                               offset(log(No.students)), family=nbinom2, data = dat)
summary(observedmodossAll)

fixef(observedmodossAll)$cond
exp(fixef(observedmodossAll)$cond)

b <- 5000
coefs1 <- matrix(ncol = 9, nrow = b)
# No observations from dat are missing values on the following, so bootstrapping can occur on data as-is: 'tot_OSSRate', 'Montessori', 'charter', 'magnet', 'percBlack', 'percWhite', 'percHispanic', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(dat), replace = T)
  nb_dat <- dat[indexes, ]
  mod <- glmmTMB(tot_OSSCount ~ Montessori + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                   offset(log(No.students)), family=nbinom2, data = nb_dat)
  coefs1[i, ] <- fixef(mod)$cond # modify to be coefficients of interest
}
colnames(coefs1) <- names(fixef(mod)$cond)
# 95% bias-corrected CIs for each coefficient
osscis <- paste0('bc_ci(coefs1[, ', 1:9, '], .95, fixef(observedmodossAll)$cond[', 1:9, '])')
lapply(osscis, function(x) eval(parse(text = x)))

# ISS ##########
observedmodissAll <-glmmTMB(tot_ISSCount ~ Montessori + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                              offset(log(No.students)), family=nbinom2, data = dat)
summary(observedmodissAll)

fixef(observedmodissAll)$cond
exp(fixef(observedmodissAll)$cond)

coefs2 <- matrix(ncol = 9, nrow = b)
# No observations from dat are missing values on the following, so bootstrapping can occur on data as-is: 'tot_ISSRate', 'Montessori', 'charter', 'magnet', 'percBlack', 'percWhite', 'percHispanic', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(dat), replace = T)
  nb_dat <- dat[indexes, ]
  mod <- glmmTMB(tot_ISSCount ~ Montessori + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                   offset(log(No.students)), family=nbinom2, data = nb_dat)
  coefs2[i, ] <- fixef(mod)$cond # modify to be coefficients of interest
}
colnames(coefs2) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Blackperc', 'Whiteperc', 'Hispanicperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
isscis <- paste0('bc_ci(coefs2[, ', 1:9, '], .95, fixef(observedmodissAll)$cond[', 1:9, '])')
lapply(isscis, function(x) eval(parse(text = x)))

########## Racial disparities in OSS MLM  ##########################################################################################################

long_count_oss <- pivot_longer(dat, cols = c('Black_num', 'Hispanic_num', 'White_num',
                                             'Black_count_oss', 'White_count_oss', 'Hispanic_count_oss'),
                               names_to = c('group','.value'),
                               names_pattern = '(Black|White|Hispanic)_(.+)')

# Set White as reference level
long_count_oss$group <- factor(long_count_oss$group)
long_count_oss$group <- relevel(long_count_oss$group, ref = 'White')

# Setup
nboot <- 5000
# Matrix to stash bootstrapped coefficients in
boot_coefs_oss <- matrix(nrow = nboot, ncol = 6)
# Data frame for bootstrap predictions
pred_df_oss <- data.frame(white_OSS_rate0 = NA, white_OSS_rate1 = NA, black_OSS_rate0 = NA,
                          black_OSS_rate1 = NA, hispanic_OSS_rate0 = NA, hispanic_OSS_rate1 = NA)

# Only bootstrap and run model on complete (usable) data
# Have to remove rows where there are zero students in a given racial group -- 
# can't use these in the offset because you can't take the log of zero 
complete_long_dat_oss <- long_count_oss[which(long_count_oss$num > 0 & is.na(long_count_oss$Montessori) == F & is.na(long_count_oss$count_oss) == F), ]

# Identify unique schools in order to bootstrap clusters
school_names_oss <- unique(complete_long_dat_oss$schoolname)
# Vector in which to stash Ns of bootstrap replicates
boot_ns_oss <- vector(length = nboot)

# Negative binomial model
fit_mod_oss <- glmmTMB(count_oss ~ Montessori*group + offset(log(num)) + (1|schoolname),
                       data=complete_long_dat_oss,
                       family=nbinom2)
summary(fit_mod_oss)

fixef(fit_mod_oss)$cond
exp(fixef(fit_mod_oss)$cond)

# Doubly robust model 
mlm_mod_dr_oss <- glmmTMB(count_oss ~ Montessori*group + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                            offset(log(num)) + (1|schoolname),
                          data=complete_long_dat_oss,
                          family=nbinom2)
summary(mlm_mod_dr_oss)

# Maximum absolute coefficient change between core model and doubly robust model:
abs(fixef(fit_mod_oss)$cond[c('(Intercept)', 'groupBlack', 'groupHispanic', 'Montessori',
                              'Montessori:groupBlack', 'Montessori:groupHispanic')] - 
      fixef(mlm_mod_dr_oss)$cond[c('(Intercept)', 'groupBlack', 'groupHispanic', 'Montessori',
                                   'Montessori:groupBlack', 'Montessori:groupHispanic')])

## Checking residuals and dispersion 
simulationOutputOSS <- simulateResiduals(fittedModel = fit_mod_oss)
plot(simulationOutputOSS)

testDispersion(fit_mod_oss)

# Bootstrap
for (i in 1:nboot) {
  cat('\rOSS bootstrap:', i, '...')
  
  boot_schools <- sample(school_names_oss, replace = T)
  boot_dat <- data.frame()
  for (j in 1:length(boot_schools)) {
    to_add <- complete_long_dat_oss[complete_long_dat_oss$schoolname == boot_schools[j], ]
    to_add$schoolname <- paste0(to_add$schoolname, '_', j)
    boot_dat <- rbind(boot_dat, to_add)
  }
  colnames(boot_dat) <- colnames(complete_long_dat_oss)
  boot_ns_oss[i] <- nrow(boot_dat)
  
  boot_mod <- glmmTMB(count_oss ~ Montessori*group + offset(log(num)) + (1|schoolname), data=boot_dat, family=nbinom2)
  boot_coefs_oss[i, ] <- fixef(boot_mod)$cond
  
  temp_preds <- ggpredict(boot_mod, terms = c('group', 'Montessori'), condition = c('num' = 100))
  temp_preds$x <- paste0(temp_preds$x, temp_preds$group)
  temp_preds <- t(temp_preds[, c('x', 'predicted')])
  pred_df_oss[i, ] <- temp_preds[2, ]
}
colnames(boot_coefs_oss) <- names(fixef(boot_mod)$cond)

# 95% bias-corrected CIs for each coefficient
cis_oss <- paste0('bc_ci(boot_coefs_oss[, ', 1:6, '], .95, fixef(fit_mod_oss)$cond[', 1:6, '])')
lapply(cis_oss, function(x) eval(parse(text = x)))
# Observed coefficients
fixef(fit_mod_oss)

# Bootstrap N info
range(boot_ns_oss)
median(boot_ns_oss)
IQR(boot_ns_oss)
quantile(boot_ns_oss, 3/4)
quantile(boot_ns_oss, 1/4)

########## Racial disparities in ISS MLM ##########################################################################################################

long_count_iss <- pivot_longer(dat, cols = c('Black_num', 'Hispanic_num', 'White_num',
                                             'Black_count_iss', 'White_count_iss', 'Hispanic_count_iss'),
                               names_to = c('group','.value'),
                               names_pattern = '(Black|White|Hispanic)_(.+)')

# Set White as reference level
long_count_iss$group <- factor(long_count_iss$group)
long_count_iss$group <- relevel(long_count_iss$group, ref = 'White')

# Setup
nboot <- 5000
# Matrix to stash bootstrapped coefficients in
boot_coefs_iss <- matrix(nrow = nboot, ncol = 6)
# Data frame for bootstrap predictions
pred_df_iss <- data.frame(white_iss_rate0 = NA, white_iss_rate1 = NA, black_iss_rate0 = NA,
                          black_iss_rate1 = NA, hispanic_iss_rate0 = NA, hispanic_iss_rate1 = NA)

# Only bootstrap and run model on complete (usable) data
# Have to remove rows where there are zero students in a given racial group -- 
# can't use these in the offset because you can't take the log of zero 
complete_long_dat_iss <- long_count_iss[which(long_count_iss$num > 0 & is.na(long_count_iss$Montessori) == F & is.na(long_count_iss$count_iss) == F), ]

# Identify unique schools in order to bootstrap clusters
school_names_iss <- unique(complete_long_dat_iss$schoolname)
# Vector in which to stash Ns of bootstrap replicates
boot_ns_iss <- vector(length = nboot)

# Negative binomial model
fit_mod_iss <- glmmTMB(count_iss ~ Montessori*group + offset(log(num)) + (1|schoolname),
                       data=complete_long_dat_iss,
                       family=nbinom2)
summary(fit_mod_iss)

fixef(fit_mod_iss)$cond
exp(fixef(fit_mod_iss)$cond)

# Doubly robust model 
mlm_mod_dr_iss <- glmmTMB(count_iss ~ Montessori*group + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                            offset(log(num)) + (1|schoolname),
                          data=complete_long_dat_iss,
                          family=nbinom2)
summary(mlm_mod_dr_iss)

# Maximum absolute coefficient change between core model and doubly robust model:
abs(fixef(fit_mod_iss)$cond[c('(Intercept)', 'groupBlack', 'groupHispanic', 'Montessori',
                              'Montessori:groupBlack', 'Montessori:groupHispanic')] - 
      fixef(mlm_mod_dr_iss)$cond[c('(Intercept)', 'groupBlack', 'groupHispanic', 'Montessori',
                                   'Montessori:groupBlack', 'Montessori:groupHispanic')])

## Checking residuals and dispersion 
simulationOutputISS <- simulateResiduals(fittedModel = fit_mod_iss)
plot(simulationOutputISS)

testDispersion(fit_mod_iss)

# Bootstrap
for (i in 1:nboot) {
  cat('\riss bootstrap:', i, '...')
  
  boot_schools <- sample(school_names_iss, replace = T)
  boot_dat <- data.frame()
  for (j in 1:length(boot_schools)) {
    to_add <- complete_long_dat_iss[complete_long_dat_iss$schoolname == boot_schools[j], ]
    to_add$schoolname <- paste0(to_add$schoolname, '_', j)
    boot_dat <- rbind(boot_dat, to_add)
  }
  colnames(boot_dat) <- colnames(complete_long_dat_iss)
  boot_ns_iss[i] <- nrow(boot_dat)
  
  boot_mod <- glmmTMB(count_iss ~ Montessori*group + offset(log(num)) + (1|schoolname), data=boot_dat, family=nbinom2)
  boot_coefs_iss[i, ] <- fixef(boot_mod)$cond
  
  temp_preds <- ggpredict(boot_mod, terms = c('group', 'Montessori'), condition = c('num' = 100))
  temp_preds$x <- paste0(temp_preds$x, temp_preds$group)
  temp_preds <- t(temp_preds[, c('x', 'predicted')])
  pred_df_iss[i, ] <- temp_preds[2, ]
}
colnames(boot_coefs_iss) <- names(fixef(fit_mod_iss)$cond)

# 95% bias-corrected CIs for each coefficient
cis_iss <- paste0('bc_ci(boot_coefs_iss[, ', 1:6, '], .95, fixef(fit_mod_iss)$cond[', 1:6, '])')
lapply(cis_iss, function(x) eval(parse(text = x)))
# Observed coefficients
fixef(fit_mod_iss)

# Bootstrap N info
range(boot_ns_iss)
median(boot_ns_iss)
IQR(boot_ns_iss)
quantile(boot_ns_iss, 3/4)
quantile(boot_ns_iss, 1/4)

########## Plots  ##########################################################################################################

## OSS ---------- ---------- 
# 95% bias-corrected CIs for predictions
pred_df_oss <- apply(pred_df_oss, 2, as.numeric)
cis_preds_oss <- paste0("bc_ci(pred_df_oss[, ", 1:6, "], .95, ggpredict(fit_mod_oss, terms = c('group', 'Montessori'), condition = c('num' = 100))$predicted[", 1:6, "])")
lapply(cis_preds_oss, function(x) eval(parse(text = x))) # NOTE: Output here is White, White, Black, Black, Hisp., Hisp.
pred_cis_oss <- sapply(cis_preds_oss, function(x) eval(parse(text = x))) # NOTE: Output here is White, White, Black, Black, Hisp., Hisp.
colnames(pred_cis_oss) <- c('WhiteNonMont', 'WhiteMont', 'BlackNonMont', 'BlackMont', 'HispanicNonMont', 'HispanicMont')
pred_cis_oss <- pred_cis_oss[, c(1, 3, 5, 2, 4, 6)] # Reorder columns to match bar order in plot

group_pred_nonmont_oss <- ggpredict(fit_mod_oss, terms = 'group', condition = c(Montessori = 0, 'num' = 100))
group_pred_nonmont_oss$mont <- 'Non-Montessori'
group_pred_mont_oss <- ggpredict(fit_mod_oss, terms = 'group', condition = c(Montessori = 1, 'num' = 100))
group_pred_mont_oss$mont <- 'Montessori'
mlm_pred_oss <- rbind(group_pred_nonmont_oss, group_pred_mont_oss)
mlm_pred_oss$group <- mlm_pred_oss$x

mlm_pred_oss <- cbind(mlm_pred_oss, t(pred_cis_oss))
colnames(mlm_pred_oss)[(ncol(mlm_pred_oss)-1):ncol(mlm_pred_oss)] <- c('lower', 'upper')

ossplot <- ggplot(mlm_pred_oss, aes(x = mont, y = predicted, color = group)) +
  geom_point(stat = 'identity', position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = .5), width = .25) +
  labs(x = "School type", y = "Model-estimated OSS counts per 100 students") +
  scale_color_manual('Racial group', labels = c('White', 'Black', 'Hispanic'),
                     values = c('#3d405b', '#e07a5f', '#81b29a')) + ylim(0, 5.5) + theme_bw() + 
  theme(text = element_text(size = 15))

## ISS ---------- ---------- 
# 95% bias-corrected CIs for predictions
pred_df_iss <- apply(pred_df_iss, 2, as.numeric)
cis_preds_iss <- paste0("bc_ci(pred_df_iss[, ", 1:6, "], .95, ggpredict(fit_mod_iss, terms = c('group', 'Montessori'), condition = c('num' = 100))$predicted[", 1:6, "])")
lapply(cis_preds_iss, function(x) eval(parse(text = x))) # NOTE: Output here is White, White, Black, Black, Hisp., Hisp.
pred_cis_iss <- sapply(cis_preds_iss, function(x) eval(parse(text = x))) # NOTE: Output here is White, White, Black, Black, Hisp., Hisp.
colnames(pred_cis_iss) <- c('WhiteNonMont', 'WhiteMont', 'BlackNonMont', 'BlackMont', 'HispanicNonMont', 'HispanicMont')
pred_cis_iss <- pred_cis_iss[, c(1, 3, 5, 2, 4, 6)] # Reorder columns to match bar order in plot

group_pred_nonmont_iss <- ggpredict(fit_mod_iss, terms = 'group', condition = c(Montessori = 0, 'num' = 100))
group_pred_nonmont_iss$mont <- 'Non-Montessori'
group_pred_mont_iss <- ggpredict(fit_mod_iss, terms = 'group', condition = c(Montessori = 1, 'num' = 100))
group_pred_mont_iss$mont <- 'Montessori'
mlm_pred_iss <- rbind(group_pred_nonmont_iss, group_pred_mont_iss)
mlm_pred_iss$group <- mlm_pred_iss$x

mlm_pred_iss <- cbind(mlm_pred_iss, t(pred_cis_iss))
colnames(mlm_pred_iss)[(ncol(mlm_pred_iss)-1):ncol(mlm_pred_iss)] <- c('lower', 'upper')

issplot <- ggplot(mlm_pred_iss, aes(x = mont, y = predicted, color = group)) +
  geom_point(stat = 'identity', position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = .5), width = .25) +
  labs(x = "School type", y = "Model-estimated ISS counts per 100 students") +
  scale_color_manual('Racial group', labels = c('White', 'Black', 'Hispanic'),
                     values = c('#3d405b', '#e07a5f', '#81b29a')) + ylim(0, 5.5) + theme_bw() + 
  theme(text = element_text(size = 15)) 


# ggsave(filename = 'oss-NBmlm.png', ossplot)
# ggsave(filename = 'iss-NBmlm.png', issplot)
