# Functions for simulating unbalanced data with group membership
# 
# Unlike average_treatment_effect, uses only the AIPW methodology. Adds in one-way standard error clustering.
#
# Parameters:
#   num.ids   The number of groups
#   n.obs     Total number of simulated observations across all groups
#   sim.seed  A random seed for replicability
#   te.fun    Function that computes treatment effect for each observation
#               One of: fixed.treatment.effect, discrete.treatment.effect, continuous.treatment.effect
#   y.fun     Function that computes the outcome for each observation
#               One of: linear.outcome or complex.outcome
#
#   The functions used by te.fun and y.fun are also included in this file

simulate.data <- function(num.ids, n.obs, sim.seed,
                          te.fun, y.fun) {
  set.seed(sim.seed)
  
  # set number per group
  mean.num.periods = n.obs / num.ids
  n.per.id = rnbinom(num.ids, (mean.num.periods - 3), 0.5) + 3
  n = sum(n.per.id)
  if (n < n.obs) {
    n = n.obs - n
    if (n < num.ids) {
      n.per.id[1:n] = n.per.id[1:n] + 1
    } else {
      n1 = floor(n / num.ids) # to add to everyone
      n2 = n %% num.ids # first n2 ids get one more
      n <- n - n1 * num.ids
      n.per.id = n.per.id + n1
      n.per.id[1:n] = n.per.id[1:n] + 1
    }
  } else if (n > n.obs) {
    n = n - n.obs
    n.ids.3 = sum(n.per.id > 3) # how many ids have more than 3 observations (so we can subtract)
    if (n < n.ids.3) {
      n.per.id[n.per.id > 3][1:n] = n.per.id[n.per.id > 3][1:n] - 1
    } else {
      n1 = floor(n / n.ids.3) # to subtract from as many as we cann
      n2 = n %% n.ids.3 # a few more will get one less
      n.per.id[n.per.id > 3] = n.per.id[n.per.id > 3] - n1
      n.per.id[n.per.id > 3][1:n2] = n.per.id[n.per.id > 3][1:n2] - 1
    }
  }
  n = sum(n.per.id)
  
  # id-level characteristics
  fixed.effect.y = rnorm(num.ids) # y fixed effect
  fixed.effect.x1 = rnorm(num.ids) # x1 fixed effect
  fixed.effect.x2 = rnorm(num.ids) # x2 fixed effect
  fixed.effect.x3 = rnorm(num.ids) # x3 fixed effect
  fixed.effect.x4 = rnorm(num.ids) # x4 fixed effect
  id.y.sig = runif(num.ids, 0.5, 1.5) # heteroscedastic residual variance
  id.x1.sig = runif(num.ids, 0.5, 1.5) # heteroscedastic residual variance: x1
  id.x2.sig = runif(num.ids, 0.5, 1.5) # heteroscedastic residual variance: x2
  id.x3.sig = runif(num.ids, 0.5, 1.5) # heteroscedastic residual variance: x3
  id.x4.sig = runif(num.ids, 0.5, 1.5) # heteroscedastic residual variance: x4
  fixed.treated.prob = runif(num.ids, 0, 0.5) # fixed treated probability
  
  # empty dataframe
  sim.data = data.frame(id = rep(0, n),
                        t = rep(0, n),
                        y = rep(0, n),
                        x1 = rep(0, n),
                        x2 = rep(0, n),
                        x3 = rep(0, n),
                        x4 = rep(0, n),
                        treated = rep(0,n),
                        true_effect = rep(0,n))
  
  # starting values for xs and other randoms
  x1.start = rnorm(num.ids)
  x2.start = runif(num.ids, -2, 2)
  x3.start = rnorm(num.ids)
  x4.start = runif(num.ids, -2, 2)
  
  # residual draws for the entire sample
  y.eps = rnorm(n)
  x1.eps = rnorm(n)
  x2.eps = rnorm(n)
  x3.eps = rnorm(n)
  x4.eps = rnorm(n)
  
  # used to determine treatment
  treated.eps = runif(n, 0, 1)
  
  # simulate data
  counter = 1
  for (i in 1:num.ids) {
    sim.data$id[counter] = i
    sim.data$t[counter] = 1
    sim.data$x1[counter] = x1.start[i]
    sim.data$x2[counter] = x2.start[i]
    sim.data$x3[counter] = x3.start[i]
    sim.data$x4[counter] = x4.start[i]
    sim.data$treated[counter] = 1 * (treated.eps[counter] < fixed.treated.prob[i])
    sim.data$true_effect[counter] = te.fun(fixed.effect = fixed.effect.y[i], X = sim.data[counter,])
    sim.data$y[counter] = y.fun(X = sim.data[counter,], fixed.effect = fixed.effect.y[i], eps = y.eps[counter], eps.sig = id.y.sig[i])
    counter = counter + 1
    for (t in 2:n.per.id[i]) {
      sim.data$id[counter] = i
      sim.data$t[counter] = t
      sim.data$x1[counter] = fixed.effect.x1[i] + 0.2 * sim.data$x2[counter-1] - 0.1 * sim.data$x3[counter-1] + x1.eps[counter] * id.x1.sig[i]
      sim.data$x2[counter] = fixed.effect.x2[i] + 0.2 * sim.data$x3[counter-1] - 0.1 * sim.data$x4[counter-1] + x2.eps[counter] * id.x2.sig[i]
      sim.data$x3[counter] = fixed.effect.x3[i] + 0.2 * sim.data$x4[counter-1] - 0.1 * sim.data$x1[counter-1] + x3.eps[counter] * id.x3.sig[i]
      sim.data$x4[counter] = fixed.effect.x4[i] + 0.2 * sim.data$x1[counter-1] - 0.1 * sim.data$x2[counter-1] + x4.eps[counter] * id.x4.sig[i]
      if (sim.data$treated[counter-1] == 0) {
        sim.data$treated[counter] = 1 * (treated.eps[counter] < (fixed.treated.prob[i] + 0.01 * fixed.effect.y[i] - 0.02 * sim.data$x1[counter] - 0.15))
      } else {
        sim.data$treated[counter] = 1 * (treated.eps[counter] < (fixed.treated.prob[i] + 0.01 * fixed.effect.y[i] - 0.02 * sim.data$x1[counter]))
      }
      sim.data$true_effect[counter] = te.fun(fixed.effect = fixed.effect.y[i], X = sim.data[counter,])
      sim.data$y[counter] = y.fun(X = sim.data[counter,], fixed.effect = fixed.effect.y[i], eps = y.eps[counter], eps.sig = id.y.sig[i])
      counter = counter + 1
    }
  }
  return(sim.data)
}

### functions that determine treatment effect and outcome
fixed.treatment.effect <- function(fixed.effect, X) {
  return(0.50)
}
discrete.treatment.effect <- function(fixed.effect, X) {
  return(0.25 + (fixed.effect < 0) * (0.25 + 0.25 * (X$x2 > 0)))
}
continuous.treatment.effect <- function(fixed.effect, X) {
  return(0.50 - 0.10 * X$x1 + 0.25 * (X$x2 > 0) - 0.05 * X$x3 * (X$x4 < 0) + 0.05 * (fixed.effect < 0))
}
linear.outcome <- function(X, fixed.effect, eps, eps.sig) {
  return(fixed.effect + X$treated * X$true_effect + 0.025 * X$t + X$x1 - X$x2 + X$x3 - X$x4 + eps * eps.sig)
}
complex.outcome <- function(X, fixed.effect, eps, eps.sig) {
  return(fixed.effect + X$treated * X$true_effect + 0.05 * X$t - 0.001 * X$t^2 +
           0.50 * (X$x1 < 0) * (0.25 + X$x2 + X$x3 * X$x4 - X$x4^2) +
           0.50 * (1 * (X$x2 > 0) - X$x3 * (X$x4 < 0) + 0.50 * (X$x1 > 1)) +
           eps * eps.sig)
}
