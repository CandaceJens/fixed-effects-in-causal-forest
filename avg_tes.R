# Estimate doubly-robust average treatment effects as in the average_treament_effect function from the grf package
# 
# Unlike average_treatment_effect, uses only the AIPW methodology. Adds in one-way standard error clustering.
#
# Parameters:
#   forest    The trained forest
#   cluster   Vector containing identifiers for the groups for clustering standard errors. Optional.
#   subset    Specifies the subset of the training example over which to compute average effects. Optional.

avg.tes <- function(forest, cluster = NULL, subset = NULL) {
  require(grf)
  require(sandwich)
  require(Matrix)

  tau.hat.pointwise <- forest$predictions
  n = length(tau.hat.pointwise)

    if (is.null(cluster)) { # if errors are not clustered, use the base function from grf
      if (is.null(subset)) { # estimate effects for the whole sample
        temp.ate <- average_treatment_effect(forest, target.sample = "all")
        ate <- c(estimate = temp.ate[1], 
                 se = temp.ate[2], 
                 t = temp.ate[1] / temp.ate[2], p = 2 * (1 - pt(abs(temp.ate[1]) / temp.ate[2], df = n - 1))
                )
        temp.atc <- average_treatment_effect(forest, target.sample = "control")
        atc <- c(estimate = temp.atc[1], 
                 se = temp.atc[2], 
                 t = temp.atc[1] / temp.atc[2], p = 2 * (1 - pt(abs(temp.atc[1]) / temp.atc[2], df = n - 1))
        )
        temp.att <- average_treatment_effect(forest, target.sample = "treated")
        att <- c(estimate = temp.att[1], 
                 se = temp.att[2], 
                 t = temp.att[1] / temp.att[2], p = 2 * (1 - pt(abs(temp.att[1]) / temp.att[2], df = n - 1))
        )
        temp.ato <- average_treatment_effect(forest, target.sample = "overlap")
        ato <- c(estimate = temp.ato[1], 
                 se = temp.ato[2], 
                 t = temp.ato[1] / temp.ato[2], p = 2 * (1 - pt(abs(temp.ato[1]) / temp.ato[2], df = n - 1))
        )
      } else { # estimate effects for a specified subset of the data
        n = sum(subset)
        temp.ate <- average_treatment_effect(forest, target.sample = "all", subset = subset)
        ate <- c(estimate = temp.ate[1], 
                 se = temp.ate[2], 
                 t = temp.ate[1] / temp.ate[2], p = 2 * (1 - pt(abs(temp.ate[1]) / temp.ate[2], df = n - 1))
        )
        temp.atc <- average_treatment_effect(forest, target.sample = "control", subset = subset)
        atc <- c(estimate = temp.atc[1], 
                 se = temp.atc[2], 
                 t = temp.atc[1] / temp.atc[2], p = 2 * (1 - pt(abs(temp.atc[1]) / temp.atc[2], df = n - 1))
        )
        temp.att <- average_treatment_effect(forest, target.sample = "treated", subset = subset)
        att <- c(estimate = temp.att[1], 
                 se = temp.att[2], 
                 t = temp.att[1] / temp.att[2], p = 2 * (1 - pt(abs(temp.att[1]) / temp.att[2], df = n - 1))
        )
        temp.ato <- average_treatment_effect(forest, target.sample = "overlap", subset = subset)
        ato <- c(estimate = temp.ato[1], 
                 se = temp.ato[2], 
                 t = temp.ato[1] / temp.ato[2], p = 2 * (1 - pt(abs(temp.ato[1]) / temp.ato[2], df = n - 1))
        )
      }
      out.df <- data.frame(estimate = c(ate[1],atc[1],att[1],ato[1]),
                           std.err =  c(ate[2],atc[2],att[2],ato[2]),
                           t.stat =   c(ate[3],atc[3],att[3],ato[3]),
                           p.value =  c(ate[4],atc[4],att[4],ato[4]))
      
    } else { # clustered standard errors
      n.clust = length(unique(cluster))
      
      y.hat = forest$Y.hat
      y.orig = forest$Y.orig
      w.hat = forest$W.hat
      w.orig = forest$W.orig
      
      if (is.null(subset)) {
        subset.idx = 1:n
      } else {
        subset.idx = which(subset)
        n = sum(subset)
        y.hat = forest$Y.hat[subset.idx]
        y.orig = forest$Y.orig[subset.idx]
        w.hat = forest$W.hat[subset.idx]
        w.orig = forest$W.orig[subset.idx]
        tau.hat.pointwise = tau.hat.pointwise[subset.idx]
      }
      
      control.idx = which(w.orig == 0)
      treated.idx = which(w.orig == 1)
      Y.hat.0 = y.hat - w.hat * tau.hat.pointwise
      Y.hat.1 = y.hat + (1 - w.hat) * tau.hat.pointwise
      
      # ATE
      tau.raw = mean(tau.hat.pointwise)
      gamma.control.raw = 1 / (1 - w.hat[control.idx])
      gamma.treated.raw = 1 / w.hat[treated.idx]
      gamma = rep(0, n)
      gamma[control.idx] = gamma.control.raw / sum(gamma.control.raw) * n
      gamma[treated.idx] = gamma.treated.raw / sum(gamma.treated.raw) * n
      dr.correction.all = gamma * (w.orig * (y.orig - Y.hat.1) - (1 - w.orig) * (y.orig - Y.hat.0))
      dr.correction = mean(dr.correction.all)
      sigma2 = mean(dr.correction.all^2) / (n - 1)
      correction.clust = sparse.model.matrix(~factor(cluster[subset.idx]) + 0, transpose = TRUE) %*% dr.correction.all
      sigma2.clust = sum(correction.clust^2) / n^2 * n.clust / (n.clust - 1)
      ate = c(estimate = tau.raw + dr.correction, se = sqrt(sigma2), 
              t = (tau.raw + dr.correction) / sqrt(sigma2), p = 2 * (1 - pt(abs(tau.raw + dr.correction) / sqrt(sigma2), df = n - 1)),
              se.clust = sqrt(sigma2.clust), 
              t.clust = (tau.raw + dr.correction) / sqrt(sigma2.clust), p.clust = 2 * (1 - pt(abs(tau.raw + dr.correction) / sqrt(sigma2.clust), df = n.clust - 1)))
      
      # ATC
      tau.raw = mean(tau.hat.pointwise[control.idx])
      gamma.control.raw = rep(1, length(control.idx))
      gamma.treated.raw = (1 - w.hat[treated.idx]) / w.hat[treated.idx]
      gamma = rep(0, n)
      gamma[control.idx] = gamma.control.raw / sum(gamma.control.raw) * n
      gamma[treated.idx] = gamma.treated.raw / sum(gamma.treated.raw) * n
      dr.correction.all = gamma * (w.orig * (y.orig - Y.hat.1) - (1 - w.orig) * (y.orig - Y.hat.0))
      dr.correction = mean(dr.correction.all)
      sigma2 = mean(dr.correction.all^2) / (n - 1)
      correction.clust = sparse.model.matrix(~factor(cluster[subset.idx]) + 0, transpose = TRUE) %*% dr.correction.all
      sigma2.clust = sum(correction.clust^2) / n^2 * n.clust / (n.clust - 1)
      atc = c(estimate = tau.raw + dr.correction, se = sqrt(sigma2), 
              t = (tau.raw + dr.correction) / sqrt(sigma2), p = 2 * (1 - pt(abs(tau.raw + dr.correction) / sqrt(sigma2), df = n - 1)),
              se.clust = sqrt(sigma2.clust), 
              t.clust = (tau.raw + dr.correction) / sqrt(sigma2.clust), p.clust = 2 * (1 - pt(abs(tau.raw + dr.correction) / sqrt(sigma2.clust), df = n.clust - 1)))
      
      # ATT
      tau.raw = mean(tau.hat.pointwise[treated.idx])
      gamma.control.raw = w.hat[control.idx] / (1 - w.hat[control.idx])
      gamma.treated.raw = rep(1, length(treated.idx))
      gamma = rep(0, n)
      gamma[control.idx] = gamma.control.raw / sum(gamma.control.raw) * n
      gamma[treated.idx] = gamma.treated.raw / sum(gamma.treated.raw) * n
      dr.correction.all = gamma * (w.orig * (y.orig - Y.hat.1) - (1 - w.orig) * (y.orig - Y.hat.0))
      dr.correction = mean(dr.correction.all)
      sigma2 = mean(dr.correction.all^2) / (n - 1)
      correction.clust = sparse.model.matrix(~factor(cluster[subset.idx]) + 0, transpose = TRUE) %*% dr.correction.all
      sigma2.clust = sum(correction.clust^2) / n^2 * n.clust / (n.clust - 1)
      att = c(estimate = tau.raw + dr.correction, se = sqrt(sigma2), 
              t = (tau.raw + dr.correction) / sqrt(sigma2), p = 2 * (1 - pt(abs(tau.raw + dr.correction) / sqrt(sigma2), df = n - 1)),
              se.clust = sqrt(sigma2.clust), 
              t.clust = (tau.raw + dr.correction) / sqrt(sigma2.clust), p.clust = 2 * (1 - pt(abs(tau.raw + dr.correction) / sqrt(sigma2.clust), df = n.clust - 1)))
      
      # ATO
      Y.residual = y.orig - y.hat
      W.residual = w.orig - w.hat
      tau.ols = lm(Y.residual ~ W.residual)
      tau.est = coef(tau.ols)[2]
      names(tau.est) = NULL
      sigma2 = vcovHC(tau.ols)[2,2]
      sigma2.clust = vcovCL(tau.ols, cluster = cluster[subset.idx])[2,2]
      ato = c(estimate = tau.est, se = sqrt(sigma2), 
              t = tau.est / sqrt(sigma2), p = 2 * (1 - pt(abs(tau.est) / sqrt(sigma2), df = n - 1)),
              se.clust = sqrt(sigma2.clust), 
              t.clust = tau.est / sqrt(sigma2.clust), p.clust = 2 * (1 - pt(abs(tau.est) / sqrt(sigma2.clust), df = n.clust - 1)))
      
      out.df <- data.frame(estimate = c(ate[1],atc[1],att[1],ato[1]),
                           std.err =  c(ate[2],atc[2],att[2],ato[2]),
                           t.stat =   c(ate[3],atc[3],att[3],ato[3]),
                           p.value =  c(ate[4],atc[4],att[4],ato[4]),
                           clust.std.err =  c(ate[5],atc[5],att[5],ato[5]),
                           clust.t.stat =   c(ate[6],atc[6],att[6],ato[6]),
                           clust.p.value =  c(ate[7],atc[7],att[7],ato[7]))
    }
    
  rownames(out.df) <- c("ATE","ATC","ATT","ATO")
  return(out.df)
}
