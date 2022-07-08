#' Bernoulli hierarchical Beyesian model
#'
#' Specify the parameter for hierarchical Beyesian model
#' @param dat input for the fucntion, should be two colums, first one is read
#'   counts of HapI, second is sum of HapI and HapII
#' @param betaA  Shape value A for prior beta distribution
#' @param betaB  Shape value B for prior beta distribution
#' @param numSavedSteps total simulation steps
#' @param burnInSteps number of burnin steps in MCMC
#' @param adaptSteps number of adapting steps
#' @param thinSteps thin used in estimating posterior distribution
#' @param parallel should MCMC run in parallel?
#' @param nChains number of chains used. Maximum is 4
#' @param saveName prefix for saved Rdata file
#' @return MCMC simulations results as a coda object
#' @note Scripts have been modified from "Doing Bayesian Data Analysis (2nd)" by
#'   John K. Kruschke.
#' @examples
#' byesRes.A <- BernHierModel(F1.TypeA[, -1], saveName = "Family1.TypeA")
#' @export
BernHierModel <- function(dat, betaA = 2, betaB = 2, numSavedSteps = 100000,
                           burnInSteps = 4000, adaptSteps = 1000,
                           thinSteps = 1, parallel = TRUE,
                           nChains = 4, saveName = "HierModel") {
  z <- as.numeric(dat[, 1])
  z.2 <- as.numeric(dat[, 2])
  N <- z + z.2
  Nsubj <- dim(dat)[1]
  dataList = list(z = z , N = N , Nsubj = Nsubj)

  modelString = "
  model {
  for (s in 1:Nsubj) {
  z[s] ~ dbin(theta[s], N[s])
  theta[s] ~ dbeta(omega*(kappa - 2) + 1, (1 - omega)*(kappa - 2) + 1)
  }
  omega ~ dbeta(%f, %f)
  kappaMinusTwo ~ dgamma(1.105125, 0.1051249)
  kappa <- kappaMinusTwo + 2
  }
  "

  # modelString <- sprintf(modelString, betaA, betaB, shape, rate)
  modelString <- sprintf(modelString, betaA, betaB)
  writeLines(modelString, con = "TEMPmodel.txt" )

  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  initsList <- function() {
    thetaInit <- rep(0, Nsubj)
    for (i in 1:Nsubj) {
      hap1 <- as.numeric(dat[i, 1])
      hap2 <- as.numeric(dat[i, 2])
      resampleD <- sample(rep(c(1, 0), times = c(hap1, hap2)), replace = TRUE)
      thetaInit[i] <- sum(resampleD)/length(resampleD)
    }
    thetaInit <- 0.001 + 0.998 * thetaInit # keep away from 0,1
    meanThetaInit <- mean(thetaInit)
    kappaInit <- 100 # lazy, start high and let burn-in find better value
    return(list(theta = thetaInit, omega = meanThetaInit ,
                kappaMinusTwo = kappaInit - 2))
  }

  nCores <- parallel::detectCores()
  if (nCores > nChains) {
    n.sims <- nChains
  } else {
    n.sims <- nChains - 1
  }
  if (n.sims < 1) {
    n.sims <- 1
  }

  parameters <- c("theta", "omega", "kappa") # The parameters to be monitored

  if (parallel & (n.sims > 1)) {
    runJagsOut <- runjags::run.jags(method = "parallel",
                           model = "TEMPmodel.txt",
                           monitor = parameters,
                           data = dataList,
                           inits = initsList,
                           n.chains = nChains,
                           adapt = adaptSteps,
                           burnin = burnInSteps,
                           sample = ceiling(numSavedSteps/nChains),
                           thin = thinSteps,
                           summarise = FALSE,
                           plots = FALSE,
                           n.sims = n.sims)
  } else {
    runJagsOut <- runjags::run.jags(method = "rjags",
                           model = "TEMPmodel.txt",
                           monitor = parameters,
                           data = dataList,
                           inits = initsList,
                           n.chains = nChains,
                           adapt = adaptSteps,
                           burnin = burnInSteps,
                           sample = ceiling(numSavedSteps/nChains),
                           thin = thinSteps,
                           summarise = FALSE,
                           plots = FALSE)
  }

  codaSamples <- coda::as.mcmc.list(runJagsOut)

  if ( !is.null(saveName) ) {
    file1 <- paste(saveName, ".Mcmc.Rdata", sep = "")
    save(codaSamples, file = file1)
    message("--- Hierarchical Bayesian model results saved in ", file1, "! ---")
  }
  return(codaSamples)
}


