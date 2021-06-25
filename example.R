### example use of fe_extract function

## Simulate test data
source("simulate_data.R")
data <- simulate.data(num.ids = 50, n.obs = 800, sim.seed = 549731, 
                      te.fun = discrete.treatment.effect, y.fun = linear.outcome)

# extract fixed effects for one identifier
source("fe_extract.R")
data <- cbind(data,
              fe.extract(X = data[c("x1","x2","x3","x4")], Y = data$y, W = data$treated, FE = data$id,
                         other.FE = data$t, degree = 3, orthogonalize = FALSE))
        

# estimate a causal forest with and without fixed effects (for comparison)
require(grf)
no.fe.cf <- causal_forest(X = data[c("x1","x2","x3","x4","t")], Y = data$y, W = data$treated)
yes.fe.cf <- causal_forest(X = data[c("x1","x2","x3","x4","t","id.est")], Y = data$y, W = data$treated)

# standard errors
source("avg_tes.R")
avg.tes(no.fe.cf)
avg.tes(yes.fe.cf, cluster = data$id)


# extract fixed effects for two identifiers
newdata <- cbind(data,
              fe.extract(X = data[c("x1","x2","x3","x4")], Y = data$y, W = data$treated, FE = data[c("id","t")],
                         degree = 3, orthogonalize = FALSE))


# estimate a causal forest with both fixed effects (for comparison)
yes.2fe.cf <- causal_forest(X = newdata[c("x1","x2","x3","x4","t.est","id.est")], Y = newdata$y, W = newdata$treated)

# standard errors
avg.tes(yes.2fe.cf, cluster = data$id)
avg.tes(yes.2fe.cf, cluster = data$t)
