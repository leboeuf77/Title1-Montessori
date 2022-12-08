# LeBoeuf, Goldstein-Greenwood, & Lillard
# Racial discipline disproportionality in Montessori and non-Montessori schools
# This script include only the analysis using a logit model to predict 
# whether a school gave zero suspensions or not as well as regressions comparing risk differences
# and relative rate ratios in suspensions. 
# It also writes out a data frame that is used in subsequent analyses with multilevel negative binomial models

# This code includes pseudo-random processes (bootstrapping). Values reported in the paper
# are based on running the code from top to bottom with the seed set below (091395).

# Load packages
invisible(lapply(list('dplyr', 'tidyr','stringi', 'stargazer', 'ggplot2'),
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
dat <- read.csv("DiscAnalysisDat.csv")
set.seed(091395)

# Calculating the number of Black, White, and Hispanic students from the percentages
num_per_group_cols <- c('Black_num', 'Hispanic_num', 'White_num')
dat[num_per_group_cols] <- NA
dat[, num_per_group_cols] <- lapply(dat[, colnames(dat)[8:10]], function(x) round((x / 100) * dat$No.students))
# Spot check: Are the numbers correct?
dat[17, c(4, 8:10, 72:74)]

################################################################################################################################################
# The CRDC 2017-2018 reports ISS counts for each racial group for "students receiving one or more in-school suspensions", 
# and these counts are disaggregated by race and gender.
# (see documentation for more details: https://ocrdata.ed.gov/assets/downloads/2017-18%20CRDC%20Overview%20Changes%20Data%20Elements.pdf) 
# We aggregated all students who received an ISS within each racial group. 
# That means, e.g., that the number of Black students who received an ISS includes all Black students (male and female, with and without a disability). 

# One school, Charlevoix Montessori Academy for the Arts, reported having only one Black student (based on the percentage of students who
# were Black and the total number of students) in the school but reported two Black students who received an ISS, for a Black ISS rate of 2.0.
# For Charlevoix Montessori Academy for the Arts, we can reasonably infer, then, that this value outside the unit interval may be in error.
# This value does not pose issues for any analysis involving RRRs or RDs because, for reasons explained in the paper, schools that do
# not have at least 5% Black or 5% Hispanic students (depending on the groups being compared) are removed from those analyses
# (and Charlevoix does not satisfy that criterion). Nor does the value pose issues for any analysis involving OSS. The one analysis
# where the value is an issue is the MLM predicting ISS rates. There, the inclusion of the value leads to estimation based on a rate
# outside the unit interval. Because the rates we calculate are based on the number of students from each racial group that receive one or
# more ISS, we truncated the number of Black students who received one or more ISS at Charlevoix to 1, because we assume that the 
# one black student in that school received at least one ISS, though the reporting seems to be in error.
dat$in.school.suspensions_Black[dat$schoolname == "charlevoixmontessoriacademyforthearts"] <- 1
# The Black ISS rate for Charlevoix, then, shifts from 2.0 (impossible) to 1.0.

# OSS counts in the CRDC data are disaggregated by race, gender, and disability status, and they report separate counts for students
# who (1) received only 1 OSS or (2) received more than one OSS. Again, we aggregated the OSS counts such that our total represents
# the number of students from each racial group who received an OSS (male and female, with or without a disability, and whether
# they received 1 or more OSS). 
################################################################################################################################################

# Calculating ISS and OSS Rates
# All OSS
dat$tot_OSSsus <- (dat$only.one.out.of.school_Black + dat$only.one.out.of.school_White + dat$only.one.out.of.school_Hisp
                   + dat$more.than.one.out.of.school_Black + dat$more.than.one.out.of.school_White + dat$more.than.one.out.of.school_Hispanic)
dat$tot_OSSRate <- (dat$tot_OSSsus / (dat$Black_num + dat$Hispanic_num + dat$White_num))
aggregate(tot_OSSRate ~ Montessori, data = dat, FUN = mean)
aggregate(tot_OSSRate ~ Montessori, data = dat, FUN = sd)
# ISS
dat$tot_ISSsus <- (dat$in.school.suspensions_Black + dat$in.school.suspensions_Hispanic + dat$in.school.suspensions_White)
dat$tot_ISSRate <- (dat$tot_ISSsus / (dat$Black_num + dat$Hispanic_num + dat$White_num))
aggregate(tot_ISSRate ~ Montessori, data = dat, FUN = mean)
aggregate(tot_ISSRate ~ Montessori, data = dat, FUN = sd)
# OSS Rates for Black, White, and Hispanic students
dat$white_OSS_rate <- (dat$only.one.out.of.school_White + dat$more.than.one.out.of.school_White) / dat$White_num
dat$black_OSS_rate <- (dat$only.one.out.of.school_Black + dat$more.than.one.out.of.school_Black) / dat$Black_num
dat$hispanic_OSS_rate <- (dat$only.one.out.of.school_Hispanic + dat$more.than.one.out.of.school_Hispanic) / dat$Hispanic_num
# ISS Rates for Black, White, and Hispanic students
dat$white_ISS_rate <- dat$in.school.suspensions_White / dat$White_num
dat$black_ISS_rate <- dat$in.school.suspensions_Black / dat$Black_num
dat$hispanic_ISS_rate <- dat$in.school.suspensions_Hispanic / dat$Hispanic_num

# Check rate ranges
range(dat$white_OSS_rate, na.rm = T)
range(dat$black_OSS_rate, na.rm = T)
range(dat$hispanic_OSS_rate, na.rm = T)
range(dat$white_ISS_rate, na.rm = T)
range(dat$black_ISS_rate, na.rm = T)
range(dat$hispanic_ISS_rate, na.rm = T)

########## Logit for zero ISS and OSS ##########################################################################################################

# ISS
dat$zeroISS <- ifelse(dat$tot_ISSsus == 0, "zero", "nonZero")
table(dat$zeroISS)
aggregate(zeroISS ~ Montessori, data = dat, FUN = table)

dat$zeroISS <- factor(dat$zeroISS)
isslogit <- glm(zeroISS ~ Montessori + No.students + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc, data = dat, family = "binomial")
summary(isslogit)
exp(isslogit$coefficients)

# OSS
dat$zeroOSS <- ifelse(dat$tot_OSSsus == 0, "zero", "nonZero")
table(dat$zeroOSS)
aggregate(zeroOSS ~ Montessori, data = dat, FUN = table)

dat$zeroOSS <- factor(dat$zeroOSS)
osslogit <- glm(zeroOSS ~ Montessori + No.students + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc, data = dat, family = "binomial")
summary(osslogit)
exp(osslogit$coefficients)

# Creating table
# stargazer(isslogit, osslogit, type = "text",
#           dep.var.labels = c("Zero ISS", "Zero OSS"),
#           out = "logitISSOSSmodels.htm", star.char = c("*", "**", "***"), star.cutoffs = c(0.05, 0.01, 0.001))

########## Black and White Student Disparities ##########################################################################################################
b <- 5000
##### Black and White Risk Differences ------------------------------------------------------------
### Black and White Risk Difference subsets
blWhOSSrdsub <- dat[which(dat$percBlack >= 5),]
table(blWhOSSrdsub$Montessori)
blWhISSrdsub <- dat[which(dat$percBlack >= 5),]
table(blWhISSrdsub$Montessori)

blWhOSSrdsub$blWhOSSrdcol <- blWhOSSrdsub$black_OSS_rate - blWhOSSrdsub$white_OSS_rate
aggregate(blWhOSSrdcol ~ Montessori, data = blWhOSSrdsub, FUN = mean)
aggregate(blWhOSSrdcol ~ Montessori, data = blWhOSSrdsub, FUN = sd)
blWhISSrdsub$blWhISSrdcol <- blWhISSrdsub$black_ISS_rate - blWhISSrdsub$white_ISS_rate
aggregate(blWhISSrdcol ~ Montessori, data = blWhISSrdsub, FUN = mean)
aggregate(blWhISSrdcol ~ Montessori, data = blWhISSrdsub, FUN = sd)

### Black and White Risk Differences Bootstrapping
# OSS RD
observedmodssRD <- lm(blWhOSSrdcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = blWhOSSrdsub)
summary(observedmodssRD)
shapiro.test(resid(observedmodssRD))
plot(observedmodssRD, which = 1)
plot(observedmodssRD, which = 2)
coefs3 <- matrix(ncol = 8, nrow = b)
# No observations from blWhOSSrdsub are missing values on the following, so bootstrapping can occur on data as-is: 'blWhOSSrdcol', 'Montessori', 'charter', 'magnet', 'percBlack', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(blWhOSSrdsub), replace = T)
  lm_dat <- blWhOSSrdsub[indexes, ]
  mod <- lm(blWhOSSrdcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs3[i, ] <- mod$coefficients[1:8] # modify to be coefficients of interest
}
colnames(coefs3) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Blackperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
ossRDcis <- paste0('bc_ci(coefs3[, ', 1:8, '], .95, observedmodssRD$coefficients[', 1:8, '])')
lapply(ossRDcis, function(x) eval(parse(text = x)))

# ISS RD
observedmodissRD <- lm(blWhISSrdcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = blWhISSrdsub)
summary(observedmodissRD)
shapiro.test(resid(observedmodissRD))
plot(observedmodissRD, which = 1)
plot(observedmodissRD, which = 2)
coefs4 <- matrix(ncol = 8, nrow = b)
# No observations from blWhISSrdsub are missing values on the following, so bootstrapping can occur on data as-is: 'blWhISSrdcol', 'Montessori', 'charter', 'magnet', 'percBlack', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(blWhISSrdsub), replace = T)
  lm_dat <- blWhISSrdsub[indexes, ]
  mod <- lm(blWhISSrdcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs4[i, ] <- mod$coefficients[1:8] # modify to be coefficient of interest
}
colnames(coefs4) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Blackperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
issRDcis <- paste0('bc_ci(coefs4[, ', 1:8, '], .95, observedmodissRD$coefficients[', 1:8, '])')
lapply(issRDcis, function(x) eval(parse(text = x)))

##### Black and White Relative Rate Ratios ------------------------------------------------------------
### Black and White RRR subsets
blWhossRRRsub <- dat[which(dat$percBlack >= 5 & dat$white_OSS_rate > 0 & dat$black_OSS_rate > 0),]
table(blWhossRRRsub$Montessori)
blWhissRRRsub <- dat[which(dat$percBlack >= 5 & dat$white_ISS_rate > 0 & dat$black_ISS_rate > 0),]
table(blWhissRRRsub$Montessori)

blWhossRRRsub$blWhossRRRcol <- blWhossRRRsub$black_OSS_rate / blWhossRRRsub$white_OSS_rate
aggregate(blWhossRRRcol ~ Montessori, data = blWhossRRRsub, FUN = mean)
aggregate(blWhossRRRcol ~ Montessori, data = blWhossRRRsub, FUN = sd)
blWhissRRRsub$blWhissRRRcol <- blWhissRRRsub$black_ISS_rate / blWhissRRRsub$white_ISS_rate
aggregate(blWhissRRRcol ~ Montessori, data = blWhissRRRsub, FUN = mean)
aggregate(blWhissRRRcol ~ Montessori, data = blWhissRRRsub, FUN = sd)

### Black and White Relative Rate Ratios Bootstrapping
# OSS RRR
observedmodossRRR <- lm(blWhossRRRcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = blWhossRRRsub)
summary(observedmodossRRR)
shapiro.test(resid(observedmodossRRR))
plot(observedmodossRRR, which = 1)
plot(observedmodossRRR, which = 2)
coefs5 <- matrix(ncol = 8, nrow = b)
# No observations from blWhossRRRsub are missing values on the following, so bootstrapping can occur on data as-is: 'blWhossRRRcol', 'Montessori', 'charter', 'magnet', 'percBlack', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(blWhossRRRsub), replace = T)
  lm_dat <- blWhossRRRsub[indexes, ]
  mod <- lm(blWhossRRRcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs5[i, ] <- mod$coefficients[1:8] # modify to be coefficients of interest
}
colnames(coefs5) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Blackperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
ossRRRcis <- paste0('bc_ci(coefs5[, ', 1:8, '], .95, observedmodossRRR$coefficients[', 1:8, '])')
lapply(ossRRRcis, function(x) eval(parse(text = x)))

# ISS RRR
observedmodissRRR <- lm(blWhissRRRcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = blWhissRRRsub)
summary(observedmodissRRR)
shapiro.test(resid(observedmodissRRR))
plot(observedmodissRRR, which = 1)
plot(observedmodissRRR, which = 2)
coefs6 <- matrix(ncol = 8, nrow = b)
# No observations from blWhissRRRsub are missing values on the following, so bootstrapping can occur on data as-is: 'blWhissRRRcol', 'Montessori', 'charter', 'magnet', 'percBlack', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(blWhissRRRsub), replace = T)
  lm_dat <- blWhissRRRsub[indexes, ]
  mod <- lm(blWhissRRRcol ~ Montessori + charter + magnet + percBlack + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs6[i, ] <- mod$coefficients[1:8] # modify to be coefficients of interest
}
colnames(coefs6) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Blackperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
issRRRcis <- paste0('bc_ci(coefs6[, ', 1:8, '], .95, observedmodissRRR$coefficients[', 1:8, '])')
lapply(issRRRcis, function(x) eval(parse(text = x)))

########## Hispanic and White Student Disparities ##########################################################################################################

##### Hispanic and White Relative Risk Differences ------------------------------------------------------------
### Hispanic and White Risk Difference subsets
hisWhOSSrdsub <- dat[which(dat$percHispanic >= 5),]
table(hisWhOSSrdsub$Montessori)
hisWhISSrdsub <- dat[which(dat$percHispanic >= 5),]
table(hisWhISSrdsub$Montessori)

hisWhOSSrdsub$hisWhOSSrdcol <- hisWhOSSrdsub$hispanic_OSS_rate - hisWhOSSrdsub$white_OSS_rate
aggregate(hisWhOSSrdcol ~ Montessori, data = hisWhOSSrdsub, FUN = mean)
aggregate(hisWhOSSrdcol ~ Montessori, data = hisWhOSSrdsub, FUN = sd)
hisWhISSrdsub$hisWhISSrdcol <- hisWhISSrdsub$hispanic_ISS_rate - hisWhISSrdsub$white_ISS_rate
aggregate(hisWhISSrdcol ~ Montessori, data = hisWhISSrdsub, FUN = mean)
aggregate(hisWhISSrdcol ~ Montessori, data = hisWhISSrdsub, FUN = sd)

### Hispanic and White Risk Differences Bootstrapping
# OSS RD
observedmodossHISrd <- lm(hisWhOSSrdcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = hisWhOSSrdsub)
summary(observedmodossHISrd)
shapiro.test(resid(observedmodossHISrd))
plot(observedmodossHISrd, which = 1)
plot(observedmodossHISrd, which = 2)
coefs7 <- matrix(ncol = 8, nrow = b)
# No observations from hisWhOSSrdsub are missing values on the following, so bootstrapping can occur on data as-is: 'hisWhOSSrdcol', 'Montessori', 'charter', 'magnet', 'percHispanic', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(hisWhOSSrdsub), replace = T)
  lm_dat <- hisWhOSSrdsub[indexes, ]
  mod <- lm(hisWhOSSrdcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs7[i, ] <- mod$coefficients[1:8] # modify to be coefficients of interest
}
colnames(coefs7) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Hispanicperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
ossHISrdcis <- paste0('bc_ci(coefs7[, ', 1:8, '], .95, observedmodossHISrd$coefficients[', 1:8, '])')
lapply(ossHISrdcis, function(x) eval(parse(text = x)))

# ISS RD
observedmodissHISrd <- lm(hisWhISSrdcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = hisWhISSrdsub)
summary(observedmodissHISrd)
shapiro.test(resid(observedmodissHISrd))
plot(observedmodissHISrd, which = 1)
plot(observedmodissHISrd, which = 2)
coefs8 <- matrix(ncol = 8, nrow = b)
# No observations from hisWhISSrdsub are missing values on the following, so bootstrapping can occur on data as-is: 'hisWhISSrdcol', 'Montessori', 'charter', 'magnet', 'percHispanic', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(hisWhISSrdsub), replace = T)
  lm_dat <- hisWhISSrdsub[indexes, ]
  mod <- lm(hisWhISSrdcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs8[i, ] <- mod$coefficients[1:8] # modify to be coefficient of interest
}
colnames(coefs8) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Hispanicperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
issHISrdcis <- paste0('bc_ci(coefs8[, ', 1:8, '], .95, observedmodissHISrd$coefficients[', 1:8, '])')
lapply(issHISrdcis, function(x) eval(parse(text = x)))

##### Hispanic and White Relative Rate Ratios ------------------------------------------------------------
### Hispanic and White RRR subsets
hisWhossRRRsub <- dat[which(dat$percHispanic >= 5 & dat$white_OSS_rate > 0 & dat$hispanic_OSS_rate > 0),]
table(hisWhossRRRsub$Montessori)
hisWhissRRRsub <- dat[which(dat$percHispanic >= 5 & dat$white_ISS_rate > 0 & dat$hispanic_ISS_rate > 0),]
table(hisWhissRRRsub$Montessori)

hisWhossRRRsub$hisWhossRRRcol <- hisWhossRRRsub$hispanic_OSS_rate / hisWhossRRRsub$white_OSS_rate
aggregate(hisWhossRRRcol ~ Montessori, data = hisWhossRRRsub, FUN = mean)
aggregate(hisWhossRRRcol ~ Montessori, data = hisWhossRRRsub, FUN = sd)
hisWhissRRRsub$hisWhissRRRcol <- hisWhissRRRsub$hispanic_ISS_rate / hisWhissRRRsub$white_ISS_rate
aggregate(hisWhissRRRcol ~ Montessori, data = hisWhissRRRsub, FUN = mean)
aggregate(hisWhissRRRcol ~ Montessori, data = hisWhissRRRsub, FUN = sd)

### Hispanic and White Relative Rate Ratios Bootstrapping
# OSS RRR
observedmodossHISrrr <- lm(hisWhossRRRcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = hisWhossRRRsub)
summary(observedmodossHISrrr)
shapiro.test(resid(observedmodossHISrrr))
plot(observedmodossHISrrr, which = 1)
plot(observedmodossHISrrr, which = 2)
coefs9 <- matrix(ncol = 8, nrow = b)
# No observations from hisWhossRRRsub are missing values on the following, so bootstrapping can occur on data as-is: 'hisWhossRRRcol', 'Montessori', 'charter', 'magnet', 'percHispanic', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(hisWhossRRRsub), replace = T)
  lm_dat <- hisWhossRRRsub[indexes, ]
  mod <- lm(hisWhossRRRcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs9[i, ] <- mod$coefficients[1:8] # modify to be coefficient of interest
}
colnames(coefs9) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'hispanicperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
ossHISrrrcis <- paste0('bc_ci(coefs9[, ', 1:8, '], .95, observedmodossHISrrr$coefficients[', 1:8, '])')
lapply(ossHISrrrcis, function(x) eval(parse(text = x)))

# ISS RRR
observedmodissHISrrr <- lm(hisWhissRRRcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = hisWhissRRRsub)
summary(observedmodissHISrrr)
shapiro.test(resid(observedmodissHISrrr))
plot(observedmodissHISrrr, which = 1)
plot(observedmodissHISrrr, which = 2)
coefs10 <- matrix(ncol = 8, nrow = b)
# No observations from hisWhissRRRsub are missing values on the following, so bootstrapping can occur on data as-is: 'hisWhissRRRcol', 'Montessori', 'charter', 'magnet', 'percHispanic', 'percWhite', 'percIdea', 'FRLperc'
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(hisWhissRRRsub), replace = T)
  lm_dat <- hisWhissRRRsub[indexes, ]
  mod <- lm(hisWhissRRRcol ~ Montessori + charter + magnet + percHispanic + percWhite + percIdea + FRLperc, data = lm_dat)
  coefs10[i, ] <- mod$coefficients[1:8] # modify to be coefficient of interest
}
colnames(coefs10) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Hispanicperc', 'Whiteperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
issHISrrrcis <- paste0('bc_ci(coefs10[, ', 1:8, '], .95, observedmodissHISrrr$coefficients[', 1:8, '])')
lapply(issHISrrrcis, function(x) eval(parse(text = x)))


##### Writing out data frame to be used in subsequent analyses 
# write.csv(dat, "DiscAnalysisDatNBmlm.csv")




