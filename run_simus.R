# source helper functions
source("helpers.R")


# examples
# generate femdom matrix and check
target_val <- 0.7
m <- mat_fem_dom(n_ids = 15, n_fem = 10, fem_dom_val = target_val, tol = 0.05)
get_fem_dom(m = m, fem_index = which(substr(colnames(m), 1, 1) == "F"))

# thin matrix (increase prunks)
target_val <- 0.65
tm <- thin_matrix(m, total_prunk = target_val)
split_prunks(tm)

# analyse matrix
analyze_matrix(tm)$f_ranks



# do simus (takes less than a minute with the settings below)
n <- 200
simdata <- data.frame(n_ids = sample(10:40, n, TRUE), 
                      n_fem = NA, 
                      fem_dom = runif(n), 
                      prunk_tot = runif(n, 0.0, 0.8)
                      )
# make sure that there are always at least 3 males and 3 females
simdata$n_fem <- round(simdata$n_ids * runif(n, 0.1, 0.9))
simdata$n_fem[simdata$n_fem < 3] <- 3
simdata$n_fem[simdata$n_fem > simdata$n_ids - 3] <- simdata$n_ids[simdata$n_fem > simdata$n_ids - 3] - 3

# 3 measures of interest
#   correlation between true ranks and DS-based ranks (group level)
#   proportion of individuals that changed rank (group level)
#   estimated proportion of individuals that changed rank (individual level
#     but effectively the same as the proportion before in this simulation 
#     because only the intercept is estimated from a model)
simdata$full_cor <- NA
simdata$full_prop_changed <- NA
simdata$full_model_prob <- NA

simdata$fem_cor <- NA
simdata$fem_prop_changed <- NA
simdata$fem_model_prob <- NA

for (i in seq_len(n)) {
  m <- mat_fem_dom(n_ids = simdata$n_ids[i], 
                   n_fem = simdata$n_fem[i], 
                   fem_dom_val = simdata$fem_dom[i], 
                   tol = 0.05, 
                   max_tries = 100000)
  if (!is.null(m)) {
    tm <- thin_matrix(m, 
                      total_prunk = simdata$prunk_tot[i], 
                      tol = 0.05, 
                      max_tries = 100000)
    if (!is.null(tm)) {
      ares <- analyze_matrix(tm)
      
      simdata$full_cor[i] <- cor(ares$full$truth, ares$full$ds, method = "s")
      simdata$full_prop_changed[i] <- mean(ares$full$changed)
      # modelled via logistic regresssion
      # this will fail in cases with NO (or very few) changes
      r <- try(glm(changed ~ 1, data = ares$full), silent = TRUE)
      if (any(class(r) == "glm")) {
        if (any(r$residuals != 0)) {
          simdata$full_model_prob[i] <- coef(r)[1]
        }
      }
      
      simdata$fem_cor[i] <- cor(ares$f_ranks$truth, ares$f_ranks$ds, method = "s")
      simdata$fem_prop_changed[i] <- mean(ares$f_ranks$changed)
      # modelled via logistic regresssion
      # this will fail in cases with NO (or very few) changes
      r <- try(glm(changed ~ 1, data = ares$f_ranks), silent = TRUE)
      if (any(class(r) == "glm")) {
        if (any(r$residuals != 0)) {
          simdata$fem_model_prob[i] <- coef(r)[1]
        }
      }
    }
  }
  cat(i, "\r")
}



plot(simdata$prunk_tot, simdata$full_cor, xlim = c(0, 1))
plot(simdata$prunk_tot, simdata$full_prop_changed, xlim = c(0, 1), ylim = c(0, 1))
plot(simdata$prunk_tot, simdata$full_model_prob, xlim = c(0, 1), ylim = c(0, 1))

plot(simdata$prunk_tot, simdata$fem_cor, xlim = c(0, 1))
plot(simdata$prunk_tot, simdata$fem_prop_changed, xlim = c(0, 1), ylim = c(0, 1))
plot(simdata$prunk_tot, simdata$fem_model_prob, xlim = c(0, 1), ylim = c(0, 1))


