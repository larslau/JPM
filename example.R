library(TMB)
library(ggplot2)
library(optimx)

#################
# Simulate data #
#################

# Number of different tracers
n_tracer <- 4

# Number of subjects per anchor point group (0 and 100 CenTauR)
n_anchor <- 20

# Head-to-head combinations
h2h_combi <- rbind(c(1, 2),
                   c(2, 3),
                   c(3, 4))

# Number of subjects per head-to-head combination
n_h2h <- c(25, 20, 18)

# Simulate slopes and intercepts mapping 
# true CenTauR values to SUVR for the three tracers
slopes <- runif(n_tracer, min = 0.02, max = 0.04)
intercepts <- runif(n_tracer, min = 0.5, max = 1.5)

# Standard deviation on CenTauR scale of 0 anchor point group
sigma_0 <- 3
# Standard deviation on CenTauR scale of 100 anchor point group
sigma_100 <- 35 
# Standard deviation on CenTauR scale of residual variation for all tracers
sigma_r <- 5

#
# Create simulated dataset
#

dat <- data.frame(id = NA, tracer = NA, group = NA, SUVR = NA, CTR_true = NA)

# Add anchor point data
for (tracer in 1:n_tracer) {
  # 0 anchor
  CTRs <- rnorm(n_anchor, mean = 0, sd = sigma_0)
  dat <- rbind(dat,
               data.frame(id = paste0('CTR0', 'tracer', tracer, '_', 1:n_anchor), 
                          tracer = paste('Tracer', tracer), 
                          group = 'CTR0',
                          SUVR = slopes[tracer] * 
                            (CTRs + rnorm(n_anchor, mean = 0, sd = sigma_r)) + 
                            intercepts[tracer],
                          CTR_true = CTRs))
  
  # 100 anchor
  CTRs <- rnorm(n_anchor, mean = 100, sd = sigma_100)
  dat <- rbind(dat,
               data.frame(id = paste0('CTR100', 'tracer', tracer, '_', 1:n_anchor), 
                          tracer = paste('Tracer', tracer), 
                          group = 'CTR100', 
                          SUVR = slopes[tracer] * 
                            (CTRs + rnorm(n_anchor, mean = 0, sd = sigma_r)) + 
                            intercepts[tracer],
                          CTR_true = CTRs))
}

# Add head-to-head data
for (i in 1:nrow(h2h_combi)) {
  # Simulate true CenTauRs
  CTRs <- runif(n_h2h[i], min = -10, max = 250)
  
  dat <- rbind(dat,
               data.frame(id = paste0('H2H', 'tracers', h2h_combi[i, 1], h2h_combi[i, 2], '_', 1:n_h2h[i]), 
                          tracer = paste('Tracer', h2h_combi[i, 1]), 
                          group = 'H2H',
                          SUVR = slopes[h2h_combi[i, 1]] * 
                            (CTRs + rnorm(n_h2h[i], mean = 0, sd = sigma_r)) + 
                            intercepts[h2h_combi[i, 1]],
                          CTR_true = CTRs),
               data.frame(id = paste0('H2H', 'tracers', h2h_combi[i, 1], h2h_combi[i, 2], '_', 1:n_h2h[i]), 
                          tracer = paste('Tracer', h2h_combi[i, 2]), 
                          group = 'H2H', 
                          SUVR = slopes[h2h_combi[i, 2]] * 
                            (CTRs + rnorm(n_h2h[i], mean = 0, sd = sigma_r)) + 
                            intercepts[h2h_combi[i, 2]],
                          CTR_true = CTRs))
}

# Remove temporary variables
rm(CTRs, i, tracer)

# Remove missing observation
dat <- na.omit(dat)

# Plot underlying CenTauR data
ggplot(dat, aes(x = interaction(tracer, group), y = CTR_true, color = tracer)) +
  geom_point() +
  geom_line(aes(group = id)) +
  scale_x_discrete('', labels = rep(unique(dat$group), each = n_tracer)) +
  scale_color_viridis_d('')

# Plot simulated SUVR data
ggplot(dat, aes(x = interaction(tracer, group), y = SUVR, color = tracer)) +
  geom_point() +
  geom_line(aes(group = id)) +
  scale_x_discrete('', labels = rep(unique(dat$group), each = n_tracer)) +
  scale_color_viridis_d('')


################
# JPM analysis #
################

compile('JPM.cpp')
dyn.load(dynlib('JPM'))

# Make dummy variables for analysis
dat$h2h <- as.numeric(dat$group == 'H2H')
dat$c0 <- as.numeric(dat$group == 'CTR0')
dat$c100 <- as.numeric(dat$group == 'CTR100')

# Fixed mean centaur for different groups
dat$centaur_numeric <- 0 + 100 * dat$c100

# Make group factor
dat$group <- factor(dat$group)

# Construct required design matrices

# h2h subjects where a fixed effect will be used to estimate stage
X_h2h_subj <- model.matrix(~ h2h : id + 0, data = dat)
# Remove zero rows
X_h2h_subj <- X_h2h_subj[, colMeans(X_h2h_subj != 0) != 0]

# Anchor subjects where a random effect will be used to estimate stage
Z_subj <- model.matrix(~ (c0 + c100) : id + 0, data = dat) 
# Remove zero rows
Z_subj <- Z_subj[, colMeans(Z_subj != 0) != 0]

# Reorder columns to match data order
ids <- gsub('.*:id', '', colnames(Z_subj))
Z_subj <- Z_subj[, match(dat$id[dat$id %in% ids], ids)]

fit_dat <- list(SUVR = dat$SUVR, 
                X_tracer = model.matrix(~ tracer + 0, data = dat),
                x_centaur = model.matrix(~ centaur_numeric + 0, data = dat),
                X_h2h_subj = X_h2h_subj,
                Z_subj = Z_subj,
                V = model.matrix(~ 1, data = dat),
                V_subj = model.matrix(~ droplevels(group) + 0, 
                                      data = subset(dat, h2h == 0 & !duplicated(id)))) 

rm(X_h2h_subj, Z_subj)

parameters <- list(slope_tracer = rep(0.01, ncol(fit_dat$X_tracer)),
                   intercept_tracer = rep(1, ncol(fit_dat$X_tracer)),
                   centaur_h2h_subj = rep(0, ncol(fit_dat$X_h2h_subj)),
                   centaur_subj = rep(0, ncol(fit_dat$Z_subj)),
                   log_sigma_tracer = rep(1, ncol(fit_dat$V)),
                   log_sigma_subj_group = c(log(5), log(30))) 

model_CTR <- MakeADFun(data = fit_dat, 
                       parameters = parameters,
                       DLL = 'JPM',
                       random = c('centaur_subj'))

# Initial likelihood optimization
fit_CTR0 <- optimx(par = model_CTR$par,
                   fn = function(par) as.numeric(model_CTR$fn(par)),
                   gr = model_CTR$gr,
                   method = c('BFGS', 'Nelder-Mead', 'nlminb'),
                   itnmax = 10000)

# Restart optimization from best parameter set
fit_CTR <- nlminb(start = fit_CTR0[which.min(fit_CTR0$value), 1:length(fit_CTR$par)], 
                  objective = model_CTR$fn, 
                  gradient = model_CTR$gr,
                  control = list(eval.max = 1e5,
                                 iter.max = 1e5)) 

#
# Inspect results
#

# Estimated variance parameters
exp(fit_CTR$par[grepl('sigma', names(fit_CTR$par))])
# vs. true
c(sigma_r, sigma_0, sigma_100)

# Estimated slopes
fit_CTR$par[grepl('slope', names(fit_CTR$par))]
# vs. true
slopes

# Estimated intercepts
fit_CTR$par[grepl('intercept', names(fit_CTR$par))]
# vs. true
intercepts

#
# Compute CenTauRs
#

# CenTauR estimation based on fixed effects
dat$CTR_estimated <- (dat$SUVR - 
                        fit_dat$X_tracer %*% fit_CTR$par[grepl('intercept', names(fit_CTR$par))]) /
  fit_dat$X_tracer %*% fit_CTR$par[grepl('slope', names(fit_CTR$par))]

# Plot
ggplot(dat, aes(x = CTR_true, y = CTR_estimated, color = tracer)) +
  geom_point() +
  geom_abline(slope = 1) +
  coord_fixed() +
  scale_color_viridis_d('') +
  xlab('True CenTauRs') +
  ylab('Estimated CenTauRs')
