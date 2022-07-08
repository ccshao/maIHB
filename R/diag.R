#' diagnose plots

#' Various plots to test the effiency and convergence of MCMC
#' @importFrom coda varnames
#' @param codaObject output from \code{\link{BernHierModel}}
#' @param parName parameter to be dignosised
#' @param outputF output folder
#' @return save diagnosis files in the output folder
#' @note Scripts have been modified from "Doing Bayesian Data Analysis (2nd)" by
#'   John K. Kruschke.
#' @examples
#' byesRes.A <- BernHierModel(F1.TypeA[, -1], saveName = "Family1.TypeA")
#' parameterNames.A <- coda::varnames(byesRes.A)[1]
#' diagMCMC(byesRes.A, parameterNames.A, outputF = "Family1.TypeA.parameters")
#' @export
diagMCMC <- function(codaObject, parName,
                     outputF = "diagPlot") {
  outputF <- file.path(".", outputF)
  if (!dir.exists(outputF)) {
    message("--- ", outputF, " created for storing dignosis plots. ---")
    dir.create(outputF)
  }

  DBDAplColors <- c("skyblue", "black", "royalblue", "steelblue")

  fileName <- file.path(outputF, paste0("MCMCdiag.", parName, ".pdf"))
  if (file.exists(fileName)) {
    message("--- ", fileName, " overwritten! ---")
  }

  pdf(fileName)
  par(mar = 0.5 + c(3, 4, 1, 0), oma = 0.1 + c(0, 0, 2, 0),
      mgp = c(2.25, 0.7, 0), cex.lab <- 1.5)
  layout(matrix(1:4, nrow = 2))
  # traceplot and gelman.plot are from CODA package:
  coda::traceplot(codaObject[, c(parName)], main = "", ylab = "Param. Value",
                   col = DBDAplColors)
  tryVal  <-  try(
    coda::gelman.plot(codaObject[, c(parName)], main = "", auto.layout = FALSE,
                       col = DBDAplColors)
 )
  # if it runs, gelman.plot returns a list with finite shrink values:
  if (class(tryVal)  == "try-error") {
    plot.new()
    print(paste0("Warning: coda::gelman.plot fails for ", parName))
  } else {
    if (class(tryVal) == "list" & !is.finite(tryVal$shrink[1])) {
      plot.new()
      print(paste0("Warning: coda::gelman.plot fails for ", parName))
    }
  }
  DbdaAcfPlot(codaObject, parName, plColors = DBDAplColors)
  DbdaDensPlot(codaObject, parName, plColors = DBDAplColors)
  mtext(text = parName, outer = TRUE, adj = c(0.5, 0.5), cex = 2.0)
  dev.off()
}

DbdaAcfPlot <- function(codaObject, parName,
                        plColors = NULL) {
  if (all(parName != varnames(codaObject))) {
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject)
  if (is.null(plColors)) {
    plColors = 1:nChain
  }
  xMat = NULL
  yMat = NULL
  for (cIdx in 1:nChain) {
    acfInfo = acf(codaObject[,c(parName)][[cIdx]], plot = FALSE)
    xMat = cbind(xMat,acfInfo$lag)
    yMat = cbind(yMat,acfInfo$acf)
  }
  matplot(xMat, yMat, type = "o", pch = 20, col = plColors, ylim = c(0,1),
           main = "", xlab = "Lag", ylab = "Autocorrelation")
  abline(h = 0,lty = "dashed")
  EffChnLngth = coda::effectiveSize(codaObject[, c(parName)])
  text(x = max(xMat), y = max(yMat), adj = c(1.0, 1.0), cex = 1.25,
        labels = paste("ESS  = ", round(EffChnLngth, 1)))
}

DbdaDensPlot  =  function(codaObject, parName,
                          plColors = NULL) {
  if (all(parName != varnames(codaObject))) {
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if (is.null(plColors)) {
    plColors = 1:nChain
  }
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for (cIdx in 1:nChain) {
    densInfo = density(codaObject[,c(parName)][[cIdx]])
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims, BEST::hdi(codaObject[, c(parName)][[cIdx]]))
  }
  matplot(xMat, yMat, type = "l", col = plColors,
           main = "", xlab = "Param. Value", ylab = "Density")
  abline(h = 0)
  points(hdiLims[1,], rep(0,nChain), col = plColors, pch = "|")
  points(hdiLims[2,], rep(0,nChain), col = plColors, pch = "|")
  text(mean(hdiLims), 0, "95% HDI", adj = c(0.5, -0.2))
}


