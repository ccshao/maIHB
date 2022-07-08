#' Summarize the posterior distribution
#'
#' Summarize the posterior distribution for HapBias, kappa and other parameters.
#'   The posterior distributs will be saved into a pdf file under current folder.
#' @param codaSamples output from \code{\link{BernHierModel}}
#' @param indiffZone limits of indiffZone
#' @param credMass  probility mass for highest density interval (HDI)
#' @param usexlim manually set the xlim of histogram?
#' @param xlimLow low end of xlim
#' @param xlimHigh high end of xlim
#' @param saveCSV save the results as a csv file?
#' @param saveName prefix for the saved pdf file and csv file
#' @return put the histograms from the posterior distributions into a pdf file
#'   and return a data.frame contain the summary info
#' @note Scripts have been modified from "Doing Bayesian Data Analysis (2nd)" by
#'   John K. Kruschke.
#' @examples
#' byesRes.A <- BernHierModel(F1.TypeA[, -1], saveName = "Family1.TypeA")
#' summaryRes.A <- smryMCMC(byesRes.A, saveName = "Family1.TypeA")
#' @export
smryMCMC <- function(codaSamples, indiffZone = c(0.49, 0.51),
                     credMass = 0.95, usexlim = FALSE,
                     xlimLow = 0.4, xlimHigh = 0.6, saveCSV = TRUE,
                     saveName = "HieraraBayesian") {
  mcmcMat <- data.frame(as.matrix(codaSamples, chains = TRUE), check.names = FALSE)
  Ntheta <- length(grep("theta", colnames(mcmcMat)))
  tem_1 <- gtools::mixedsort(colnames(mcmcMat))
  tem_1 <- c(tem_1[1], tem_1[3], tem_1[2], tem_1[-c(1:3)])
  mcmcMat.1 <- mcmcMat[, tem_1]
  summaryInfo <- data.frame(t(apply(mcmcMat.1[, -1], 2, summarizePost, ROPE = indiffZone)))

  postPdf <- paste0(saveName, ".Posterior.pdf")
  pdf(postPdf)
  tem_1 <- colnames(mcmcMat.1)[-1]
  if (usexlim == FALSE) {
    lapply(tem_1, function(X) {BEST::plotPost(mcmcMat.1[, X], ROPE = indiffZone,
                                              credMass = credMass,
                                              showMode = FALSE, xlab = X)})
  } else {
    lapply(tem_1, function(X) {BEST::plotPost(mcmcMat.1[, X], ROPE = indiffZone,
                                              credMass = credMass,
                                              showMode = FALSE, xlab = X,
                                              xlim = c(xlimLow, xlimHigh))})
  }
  dev.off()
  message("--- Posterior distribution saved in ", postPdf, "! ---")

  if(saveCSV) {
    postSum <- paste0(saveName, ".SummaryInfo.csv")
    write.table(summaryInfo, file = postSum, sep = "\t", quote = FALSE)
    message("--- Summary file saved in ", postSum, "! ---")
  }
  return(summaryInfo)
}


#' plot the posterior distribution
#'
#' plot the posterior distribution for a single parameter
#' @param codaSamples output from \code{\link{BernHierModel}}
#' @param indiffZone limits of indiffZone
#' @param credMass  probility mass for highest density interval (HDI)
#' @param usexlim manually set the xlim of histogram?
#' @param xlimLow low end of xlim
#' @param xlimHigh high end of xlim
#' @param pref prefix of saved pdf file
#' @note A wrapper of plotPost fucntion in BEST.
#' @examples
#' byesRes.A <- BernHierModel(F1.TypeA[, -1], saveName = "Family1.TypeA")
#' singlePoster(byesRes.A, "kappa")
#' @export
singlePoster <- function(codaSamples, para, indiffZone = c(0.49, 0.51),
                         credMass = 0.95, showMode = FALSE, usexlim = FALSE,
                         xlimLow = 0.4, xlimHigh = 0.6, pref = "single" ) {

  mcmcMat <- data.frame(as.matrix(codaSamples, chains = TRUE), check.names = FALSE)

  postPdf <- paste0(pref, ".", para, ".Posterior.pdf")
  pdf(postPdf)
  if (usexlim == TRUE) {
    BEST::plotPost(mcmcMat[, para], ROPE = indiffZone, credMass = credMass,
                   xlim = c(xlimLow, xlimHigh), showMode = FALSE, xlab = para)
  } else {
    BEST::plotPost(mcmcMat[, para], ROPE = indiffZone, credMass = credMass,
                   showMode = FALSE, xlab = para)
  }
  dev.off()
}


summarizePost <- function(paramSampleVec, ROPE = NULL, credMass = 0.95) {
  meanParam <- mean(paramSampleVec)
  medianParam <- median(paramSampleVec)
  dres <- density(paramSampleVec)
  modeParam <- dres$x[which.max(dres$y)]
  mcmcEffSz <- round(coda::effectiveSize(paramSampleVec))
  names(mcmcEffSz) <- NULL

  hdiLim <- BEST::hdi(paramSampleVec, credMass = credMass)
  if(!is.null(ROPE)) {
    pcltRope <- (100 * sum( paramSampleVec < ROPE[1])
                 / length(paramSampleVec))
    pcgtRope <- (100 * sum( paramSampleVec > ROPE[2])
                 / length(paramSampleVec))
    pcinRope <- 100 - (pcltRope + pcgtRope)
  } else {
    ROPE <- c(NA, NA)
    pcltRope <- NA
    pcgtRope <- NA
    pcinRope <- NA
  }
  return(c(Mean = meanParam, Median = medianParam, Mode = modeParam,
             ESS = mcmcEffSz, HDImass = credMass, HDIlow = as.numeric(hdiLim[1]),
             HDIhigh = as.numeric(hdiLim[2]), ROPElow = ROPE[1], ROPEhigh = ROPE[2],
             PcntLtROPE = pcltRope, PcntInROPE = pcinRope, PcntGtROPE = pcgtRope))
}









