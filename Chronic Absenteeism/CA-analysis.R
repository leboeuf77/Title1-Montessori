# Rates of Chronic Absenteeism in Title 1 Montessori and Non-Montessori schools
# Authors: Jacob Goldstein-Greenwood & Lee LeBoeuf
# Last revised: 3/03/2023

# Values from pseudo-random processes (bootstrapping) reported in the paper
# are based on running the code from top to bottom with a seed of 12

###################################### Load and read ######################################
# Libraries
invisible(lapply(list('dplyr', 'tidyr', 'stringi', 'rlist', 
                      'ggeffects', 'ggplot2', 'glmmTMB', 'DHARMa'),
                 function(pkg) library(pkg, character.only = TRUE)))

# Bias-corrected CI function
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

# Read
dat_raw <- read.csv("ChronAbsDat.csv")
# Seed
set.seed(12)

###################################### Cleaning and preparation ######################################
##### Identify and clean missing data and 0 values
# Identify columns that use '-' to represent 0/0 chronically absent students of a given group
hyphen <- data.frame(variable = apply(dat_raw, 2, function(x) any(stri_detect(x, fixed = '-'))))
# Select columns where '-' needs to be converted to 0 so that vars can be treated as numeric; then convert and make numeric
cols_numerify <- rownames(hyphen)[45:55]
dat_raw[, cols_numerify] <- apply(dat_raw[, cols_numerify], 2, function(x) ifelse(x == '-', 0, x))
dat_raw[, cols_numerify] <- apply(dat_raw[, cols_numerify], 2, function(x) as.numeric(as.character(x)))

##### Sum counts of chronically absent students across the male/female row splits for each school
dat <- dat_raw %>%
  group_by(schoolname) %>%
  mutate(AIAN_chron_ab_total = no.AmericanIndianAbs + lag(no.AmericanIndianAbs),
         asian_chron_ab_total = no.AsianAbs + lag(no.AsianAbs),
         HPI_chron_ab_total = no.HawaiianPIAbs + lag(no.HawaiianPIAbs),
         hispanic_chron_ab_total = no.HispanicAbs + lag(no.HispanicAbs),
         black_chron_ab_total = no.BlackAbs + lag(no.BlackAbs),
         white_chron_ab_total = no.WhiteAbs + lag(no.WhiteAbs),
         twoplus_chron_ab_total = no.Two.or.more.racesAbs + lag(no.Two.or.more.racesAbs),
         total_chron_ab = TotalAbs + lag(TotalAbs)) %>%
  ungroup()
# Keep only rows that that contain combined M+F numbers; keep only columns of analytic interest
dat <- dat[-which(is.na(dat$total_chron_ab)), -which(colnames(dat) == 'Sex')]
table(dat$Montessori)

# Check: Are the total chronic absenteeism counts for each student subgroup the same in the M+F summed data set and in the original M+F split data set?
# all(t(unlist(lapply(dat_raw[, 45:52], function(x) sum(x, na.rm = T)))) == t(lapply(dat[, 52:ncol(dat)], function(x) sum(x, na.rm = T))))

##### 
# Calculate number of students in each racial group 
num_per_group_cols <- c('AIAN_num', 'asian_num', 'HPI_num', 'black_num', 'hispanic_num', 'white_num', 'twoplus_num')
dat[num_per_group_cols] <- NA
dat[, num_per_group_cols] <- lapply(dat[, colnames(dat)[5:11]], function(x) round((x / 100) * dat$No.students))

# On average, how close are the sums of the by-group student counts to the reported number total students?
dat$total_est_students <- rowSums(dat[, 63:69])
dat$total_reported <- dat$No.students
mean(abs(dat$total_est_students - dat$No.students))
# Average difference of <1 student

#####
# If a school had *no* students for a subgroup, indicate NA for that school's chronic absenteeism
# (I.e., 0/0 should be NA; 0/0 is not a meaningful zero that should be included in models)
dat$hispanic_chron_ab_total <- ifelse(dat$percHispanic == 0, NA, dat$hispanic_chron_ab_total)
dat$black_chron_ab_total <- ifelse(dat$percBlack == 0, NA, dat$black_chron_ab_total)
dat$white_chron_ab_total <- ifelse(dat$percWhite == 0, NA, dat$white_chron_ab_total)

# Below, we check that the number of chronically absent students in each racial group
# does not exceed the number of students in that racial group. 
# We assume these cases are due to data entry errors with the CRDC, so in those cases,
# we truncate the number of chronically absent students to equal the number of 
# students in the racial group.

dat[which(dat$black_chron_ab_total > dat$black_num), ]
# Three schools are likely in error:
# charlevoixmontessoriacademyforthearts -- 2.2% Black; 1 Black student; 4 Black abs
# mariamontessoriacademy -- 0.5% Black; 3 Black students; 4 Black abs
# prescottvalleyschool -- 1.2% Black; 3 Black students; 5 Black abs
dat$black_chron_ab_total[dat$schoolname == "charlevoixmontessoriacademyforthearts"] <- dat$black_num[dat$schoolname == "charlevoixmontessoriacademyforthearts"]
dat$black_chron_ab_total[dat$schoolname == "mariamontessoriacademy"] <- dat$black_num[dat$schoolname == "mariamontessoriacademy"]
dat$black_chron_ab_total[dat$schoolname == "prescottvalleyschool"] <- dat$black_num[dat$schoolname == "prescottvalleyschool"]
# Post-correction check:
dat[which(dat$black_chron_ab_total > dat$black_num), ]

dat[which(dat$hispanic_chron_ab_total > dat$hispanic_num), ]
# Two schools are likely in error:
# acornmontessoricharterschoolinc.-west -- 7.7% Hispanic; 9 Hispanic students; 15 Hispanic abs
# rocklinindependentcharteracademy -- 14.6% Hispanic; 22 Hispanic students; 23 Hispanic abs 
dat$hispanic_chron_ab_total[dat$schoolname == "acornmontessoricharterschoolinc.-west"] <- dat$hispanic_num[dat$schoolname == "acornmontessoricharterschoolinc.-west"]
dat$hispanic_chron_ab_total[dat$schoolname == "rocklinindependentcharteracademy"] <- dat$hispanic_num[dat$schoolname == "rocklinindependentcharteracademy"]
# Post-correction check:
dat[which(dat$hispanic_chron_ab_total > dat$hispanic_num), ]

dat[which(dat$white_chron_ab_total > dat$white_num), ]
# None
#####

# Calculate total number of chronically absent students among the focal racial groups
dat$total_chron_ab <- rowSums(dat[, c('black_chron_ab_total', 'white_chron_ab_total', 'hispanic_chron_ab_total')], na.rm = T)
dat$total_students <- rowSums(dat[, c('black_num', 'white_num', 'hispanic_num')], na.rm = T)

# Calculating average percents of chronically absent students
dat$percChronicAbs <- dat$total_chron_ab / dat$total_students * 100
aggregate(percChronicAbs ~ Montessori, data = dat, FUN = mean)

# Average percents for each racial group
percBlackca <- dat$black_chron_ab_total / dat$black_num * 100
aggregate(percBlackca ~ Montessori, data = dat, FUN = mean)
percHispca <- dat$hispanic_chron_ab_total / dat$hispanic_num * 100
aggregate(percHispca ~ Montessori, data = dat, FUN = mean)
percWhiteca <- dat$white_chron_ab_total / dat$white_num * 100
aggregate(percWhiteca ~ Montessori, data = dat, FUN = mean)


###################################### Main analyses ######################################
##### Overall rates  ------------------------------------------------------------
##### Overall difference between school types in the number of students who were chronically absent from the combined groups of interest in the matched sample
total_diff <- glmmTMB(total_chron_ab ~ Montessori + offset(log(total_students)), family = nbinom2, data = dat)
summary(total_diff)

fixef(total_diff)$cond
exp(fixef(total_diff)$cond)

# Doubly robust model accounting for some residual omitted variable bias by including school-level controls
total_diff_dr <- glmmTMB(total_chron_ab ~ Montessori + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                           offset(log(total_students)), family = nbinom2, data = dat)
summary(total_diff_dr)

fixef(total_diff_dr)$cond
exp(fixef(total_diff_dr)$cond)

# Checking residuals and dispersion 
simOutput_overallCA <- simulateResiduals(fittedModel = total_diff_dr)
plot(simOutput_overallCA) # "No significant problems detected"
testDispersion(total_diff_dr)

# Maximum absolute coefficient change between core model and doubly robust model:
abs(fixef(total_diff)$cond[c('(Intercept)', 'Montessori')] - 
      fixef(total_diff_dr)$cond[c('(Intercept)', 'Montessori')])

# Bootstrapping
b <- 5000
coefs1 <- matrix(ncol = 9, nrow = b)
for (i in 1:b) {
  cat(i, '... ')
  indexes <- sample(nrow(dat), replace = T)
  nb_dat <- dat[indexes, ]
  mod <- glmmTMB(total_chron_ab ~ Montessori + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                   offset(log(total_students)), family = nbinom2, data = nb_dat)
  coefs1[i, ] <- fixef(mod)$cond # modify to be coefficients of interest
}
colnames(coefs1) <- c('Intercept', 'Montessori', 'Charter', 'Magnet', 'Blackperc', 'Whiteperc', 'Hispanicperc', 'Idea', 'FRLperc')
# 95% bias-corrected CIs for each coefficient
overallcis <- paste0('bc_ci(coefs1[, ', 1:9, '], .95, fixef(total_diff_dr)$cond[', 1:9, '])')
lapply(overallcis, function(x) eval(parse(text = x)))

##### MLM  ------------------------------------------------------------
##### Multilevel model: How do disparities in chronic absenteeism rates between subgroups change at Montessori vs. non-Montessori schools?

# Convert data from wide to long
long_dat <- pivot_longer(dat,
                         cols = c('black_chron_ab_total', 'hispanic_chron_ab_total', 'white_chron_ab_total',
                                  'black_num', 'hispanic_num', 'white_num'),
                         names_to = c('group','.value'),
                         names_pattern = '(black|white|hispanic)_(.+)')

# Set White as reference level
long_dat$group <- factor(long_dat$group)
long_dat$group <- relevel(long_dat$group, ref = 'white')

# Only bootstrap and run model on complete (usable) data
# Have to remove rows where there are zero students in a given racial group -- 
# can't use these in the offset because you can't take the log of zero 
complete_long_dat <- long_dat[which(long_dat$num > 0 & is.na(long_dat$Montessori) == F & is.na(long_dat$chron_ab_total) == F), ]

# Check: None of the variables used in the negative binomial model (doubly robust or not)
#        now have NAs in the long data set; all rows can be bootstrapped
apply(complete_long_dat[, c('chron_ab_total', 'Montessori', 'group', 'charter', 'magnet',
                            'percBlack', 'percWhite', 'percHispanic', 'percIdea', 'FRLperc', 'num')], 2, function(x) any(is.na(x)))

# Negative binomial model
mlmmodca <- glmmTMB(chron_ab_total ~ Montessori*group + offset(log(num)) + (1|schoolname),
                    family = nbinom2,
                    data=complete_long_dat)
summary(mlmmodca)

fixef(mlmmodca)$cond
exp(fixef(mlmmodca)$cond)

# Checking residuals and dispersion 
simulationOutput1 <- simulateResiduals(fittedModel = mlmmodca)
plot(simulationOutput1)
testDispersion(mlmmodca)
testZeroInflation(mlmmodca)

# Doubly robust model 
mlmmodca_dr <- glmmTMB(chron_ab_total ~ Montessori*group + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                            offset(log(num)) + (1|schoolname),
                       family = nbinom2,
                       data = complete_long_dat)
summary(mlmmodca_dr)

fixef(mlmmodca_dr)$cond
exp(fixef(mlmmodca_dr)$cond)

# Checking residuals and dispersion 
simulationOutput2 <- simulateResiduals(fittedModel = mlmmodca_dr)
plot(simulationOutput2)
testDispersion(mlmmodca_dr)
testZeroInflation(mlmmodca_dr)

# Maximum absolute coefficient change between core model and doubly robust model:
abs(fixef(mlmmodca)$cond[c('(Intercept)', 'Montessori', 'groupblack', 'grouphispanic', 'Montessori:groupblack', 'Montessori:grouphispanic')] - 
      fixef(mlmmodca_dr)$cond[c('(Intercept)', 'Montessori', 'groupblack', 'grouphispanic', 'Montessori:groupblack', 'Montessori:grouphispanic')])

# Getting the median values for each covariate in the doubly robust model to generate predicted values during bootstrapping and creating the plot below
apply(dat[, c('percBlack', 'percWhite', 'percHispanic', 'percIdea', 'FRLperc')], 2, function(x) mean(x))

# Seeing whether most schools in the dataset are either magnet or charter schools to generate predicted values during bootstrapping and creating the plot below
apply(dat[, c('charter', 'magnet')], 2, function(x) table(x))


##### MLM bootstrapping ----
# Setup
nboot <- 5000
# Matrix to stash bootstrapped coefficients in
boot_coefs <- matrix(nrow = nboot, ncol = 13)
# Data frame for bootstrap predictions
pred_df <- data.frame(white_abs_rate0 = NA, white_abs_rate1 = NA, black_abs_rate0 = NA,
                      black_abs_rate1 = NA, hispanic_abs_rate0 = NA, hispanic_abs_rate1 = NA)
# Identify unique schools in order to bootstrap clusters
school_names_ca <- unique(complete_long_dat$schoolname)
# Vector in which to stash Ns of bootstrap replicates
boot_ns <- vector(length = nboot)

# Bootstrap
for (i in 1:nboot) {
  cat('\rchronicabs bootstrap:', i, '...')
  
  boot_schools <- sample(school_names_ca, replace = T)
  boot_dat <- data.frame()
  for (j in 1:length(boot_schools)) {
    to_add <- complete_long_dat[complete_long_dat$schoolname == boot_schools[j], ]
    to_add$schoolname <- paste0(to_add$schoolname, '_', j)
    boot_dat <- rbind(boot_dat, to_add)
  }
  colnames(boot_dat) <- colnames(complete_long_dat)
  boot_ns[i] <- nrow(boot_dat)
  
  boot_mod <- glmmTMB(chron_ab_total ~ Montessori*group + charter + magnet + percBlack + percWhite + percHispanic + percIdea + FRLperc + 
                        offset(log(num)) + (1|schoolname),
                      family = nbinom2,
                      data = boot_dat)
  boot_coefs[i, ] <- fixef(boot_mod)$cond
  
  temp_preds <- ggpredict(boot_mod, terms = c('group', 'Montessori'), condition = c('num' = 100, 'charter' = 1, 'magnet' = 0, 'percBlack' = 26.44, 
                                                                                    'percWhite' = 40.22, 'percHispanic' = 22.51, 'percIdea' = 12.04, 'FRLperc' = 51.84))
  temp_preds$x <- paste0(temp_preds$x, temp_preds$group)
  temp_preds <- t(temp_preds[, c('x', 'predicted')])
  pred_df[i, ] <- temp_preds[2, ]
}
colnames(boot_coefs) <- names(fixef(boot_mod)$cond)

# 95% bias-corrected CIs for each coefficient
cis_ca <- paste0('bc_ci(boot_coefs[, ', 1:13, '], .95, fixef(mlmmodca_dr)$cond[', 1:13, '])')
lapply(cis_ca, function(x) eval(parse(text = x)))
# Observed coefficients
fixef(mlmmodca_dr)$cond
exp(fixef(mlmmodca_dr)$cond)

# Bootstrap N info
range(boot_ns)
median(boot_ns)
IQR(boot_ns)
quantile(boot_ns, 3/4)
quantile(boot_ns, 1/4)

##### MLM Plot ----
# 95% bias-corrected CIs for predictions
pred_df_ca <- apply(pred_df, 2, as.numeric)
cis_preds_ca <- paste0("bc_ci(pred_df_ca[, ", 1:6, "], .95, ggpredict(mlmmodca_dr, terms = c('group', 'Montessori'), condition = c('num' = 100, 'charter' = 1, 'magnet' = 0, 'percBlack' = 26.44, 
                                                                               'percWhite' = 40.22, 'percHispanic' = 22.51, 'percIdea' = 12.04, 'FRLperc' = 51.84))$predicted[", 1:6, "])")
lapply(cis_preds_ca, function(x) eval(parse(text = x))) # NOTE: Output here is White, White, Black, Black, Hisp., Hisp.
pred_cis_ca <- sapply(cis_preds_ca, function(x) eval(parse(text = x))) # NOTE: Output here is White, White, Black, Black, Hisp., Hisp.
colnames(pred_cis_ca) <- c('WhiteNonMont', 'WhiteMont', 'BlackNonMont', 'BlackMont', 'HispanicNonMont', 'HispanicMont')
pred_cis_ca <- pred_cis_ca[, c(1, 3, 5, 2, 4, 6)] # Reorder columns to match bar order in plot

group_pred_nonmont_ca <- ggpredict(mlmmodca_dr, terms = 'group', condition = c(Montessori = 0, 'num' = 100, 'charter' = 1, 'magnet' = 0, 'percBlack' = 26.44, 
                                                                               'percWhite' = 40.22, 'percHispanic' = 22.51, 'percIdea' = 12.04, 'FRLperc' = 51.84))
group_pred_nonmont_ca$mont <- 'Non-Montessori'
group_pred_mont_ca <- ggpredict(mlmmodca_dr, terms = 'group', condition = c(Montessori = 1, 'num' = 100, 'charter' = 1, 'magnet' = 0, 'percBlack' = 26.44, 
                                                                            'percWhite' = 40.22, 'percHispanic' = 22.51, 'percIdea' = 12.04, 'FRLperc' = 51.84))
group_pred_mont_ca$mont <- 'Montessori'
mlm_pred_ca <- rbind(group_pred_nonmont_ca, group_pred_mont_ca)
mlm_pred_ca$group <- mlm_pred_ca$x

mlm_pred_ca <- cbind(mlm_pred_ca, t(pred_cis_ca))
colnames(mlm_pred_ca)[(ncol(mlm_pred_ca)-1):ncol(mlm_pred_ca)] <- c('lower', 'upper')

caplot <- ggplot(mlm_pred_ca, aes(x = mont, y = predicted, fill = group)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 1), width = .25) +
  labs(x = "School type", y = "Model-estimated Chronic Absenteeism counts per 100 students") +
  scale_fill_manual('Racial group', labels = c('White', 'Black', 'Hispanic'),
                     values = c('navajowhite', 'thistle', 'salmon')) + ylim(0, 25) + theme_bw() + 
  theme(text = element_text(size = 15))

ggsave(filename = 'ca-NBmlm.png', caplot)



