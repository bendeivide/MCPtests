# References:
# Bhering, L. L. et. al. Alternative methodology for Scott-Knott test. Brazilian Society
#            of PLant Breeding, v.8, n.1, 2008, p. 9-16.

# # Temporario
# library(foreach)
# library(doParallel)
# library(parallel)
#
#
# # Response variable
# y <- rv <- c(100.08, 105.66, 97.64, 100.11, 102.60, 121.29, 100.80,
#              99.11, 104.43, 122.18, 119.49, 124.37, 123.19, 134.16,
#              125.67, 128.88, 148.07, 134.27, 151.53, 127.31)
#
# # Treatments
# trt <- treat <- factor(rep(LETTERS[1:5], each = 4))
# #trt <- factor(rep(c("boi", "vaca", "bode", "cabra", "pato"), each = 4))
#
# # ExpDes.pt::dic(trt, y, quali = TRUE, mcomp = "sk")
#
#
# # dados <- data.frame(trt, y)
# # write.table(dados, "dados.csv", sep = ";")
#
# # Anova
# res     <- anova(aov(rv~treat))
# dferror <- DFerror <- res$Df[2]
# mserror <- MSerror <- res$`Mean Sq`[2]
# replication <- 4
# alpha <- 0.05
#
# #ExpDes.pt::scottknott(y = rv, trt = treat, DFerror = dferror, SSerror = mserror, alpha = 0.05, group = TRUE, main = NULL)
#
# y <- rv
#
# trt <- treat
#
# alpha <- 0.05
#
# # Exemplo do livro de Daniel
# y <- Ybar <- c(21.16,22.72,23.58,28.14)
# trt <- c("A", "B", "c", "D")
# mserror <- 10.19
# dferror <- 20
# replication <- 6
#
# # Exemplo Bhering
# y <- c(255.5, 259.3, 271.6, 290.6, 298.8, 334.9, 341.0, 348.7, 384.3, 495.5)
# trt <- factor(1:10)
# mserror <- 1254.327
# dferror <- 18
# replication <- 3
# alpha <- 0.05


# Scott-Knott-Bhering.etal.'s test

sktest_bhering2008 <- function(y, trt, dferror, mserror, replication, alpha) {
  # Parallel activate
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)

  # Ordered means
  Ybar <- sort(tapply(y, trt, mean), decreasing = TRUE)
  # Length of means
  n <- length(Ybar)
  # Results of Scott-Knott's test
  groups <- rep(0, times = length(Ybar))

  # Temporary files
  pvalueext <- tempfile(pattern = "pvalues.", tmpdir = tempdir())
  breakgroupsext <- tempfile(pattern = "breakgroups.", tmpdir = tempdir())
  maxboext <- tempfile(pattern = "maxbo.", tmpdir = tempdir())
  s2ext <- tempfile(pattern = "s2.", tmpdir = tempdir())
  nutestext <- tempfile(pattern = "nutest.", tmpdir = tempdir())
  stattestext <- tempfile(pattern = "stattest.", tmpdir = tempdir())
  resultsext <- tempfile(pattern = "results.", tmpdir = tempdir())
  g2ext <- tempfile(pattern = "g2.", tmpdir = tempdir())


  # Function for the separation of the groups
  sg <- function(means, ...) {
    means <- sort(means, decreasing = TRUE)
    posmeans <- pmatch(means, Ybar)
    names(means) <- names(Ybar)[na.omit(posmeans)]


    tm <- length(means) - 1
    # Calculate b0
    b0 <- foreach(i = 1:tm, .combine = c) %dopar%  {
      g1 <- means[1:i]
      g2 <- means[(i + 1):length(means)]
      gall <- c(g1, g2)
      sum(g1)^2/length(g1) + sum(g2)^2/length(g2) - (sum(gall))^2/length(gall)
    }
    # break of groups
    corte <- which.max(b0)

    # two groups
    g1 <- means[1:corte]
    g2 <- means[(corte + 1):length(means)]
    tg <- c(g1,g2)

    # Calculate ML estimate of sigma²
    sig2 <- (1 / (length(tg) + dferror)) * (sum(tg^2) - (sum(tg))^2 / length(tg) + dferror * mserror / replication)

    # Test statistic
    ts <- pi / (2 * (pi - 2)) * max(b0) / sig2

    # Degrees of freedom
    nu <- length(tg) / (pi - 2)

    # P-value
    pvalue <- pchisq(ts, nu, lower.tail = FALSE)


    # Separation of the groups
    if (pvalue > alpha) {
      # Condition for pvalue > alpha
      # Classification of Scott-Knott's test
      cat(pvalue,"\n", file = pvalueext, append = TRUE)
      cat(substr(names(g1), 1, 3), "_vs_", substr(names(g2), 1, 3), ";", "\n", sep = " ", file = breakgroupsext, append = TRUE)
      cat(max(b0),"\n", file = maxboext, append = TRUE)
      cat(sig2,"\n", file = s2ext, append = TRUE)
      cat(nu,"\n", file = nutestext, append = TRUE)
      cat(ts,"\n", file = stattestext, append = TRUE)

      # for (i in 1:length(c(g1, g2))) {
      #   cat(names(c(g1, g2)[i]),"\n", file = "results",append = TRUE)
      # }
      # cat("*","\n", file = "results", append = TRUE)


      # Test into test
      if (file.exists("g2")) {
        g1 <- as.vector(read.table("g2")[[1]])
        file.remove("g2")
        if (length(g1) > 1) sg(g1)
        if (length(g1) == 1) {
          cat(names(g1),"\n", file = resultsext,append = TRUE)
        }
      }
    }

    if (pvalue <= alpha) {
      # Classification of Scott-Knott's test
      for (i in 1:length(g1)) {
        cat(names(g1[i]),"\n", file = resultsext,append = TRUE)
      }
      cat("*","\n", file = resultsext, append = TRUE)
      # Auxiliar results
      cat(pvalue,"\n", file = pvalueext, append = TRUE)
      cat(substr(names(g1), 1, 3), "_vs_", substr(names(g2), 1, 3), ";", "\n", sep = " ", file = breakgroupsext, append = TRUE)
      cat(max(b0),"\n", file = maxboext, append = TRUE)
      cat(sig2,"\n", file = s2ext, append = TRUE)
      cat(nu,"\n", file = nutestext, append = TRUE)
      cat(ts,"\n", file = stattestext, append = TRUE)
      if (length(g2) > 1) {
        for (i in 1:length(g2)) {
          cat(g2[i],"\n", file = g2ext, append = TRUE)
        }
      } else cat(g2,"\n", file = g2ext, append = TRUE)
      # Test into test
      if (length(g1) == 1) {
        g1 <- as.vector(read.table(g2ext)[[1]])
        file.remove(g2ext)
      }
      if (length(g1) > 1) Recall(g1)
    }
  }

  # Loading the separation of the groups and generating external files
  sg(Ybar)

  # Result of Scott-Knott's test
  if (file.exists(resultsext) == FALSE) {
    stop("Missing data entry!", call. = FALSE)
  }
  else{
    # Loading external file of results
    xx <- read.table(resultsext)
    # Remove external file
    file.remove(resultsext)
    x <- as.vector(xx[[1]])
    z <- 1

    # Results of Scott-Knott's test
    for (j in 1:length(x)) {
      if (x[j] == "*")	{z <- z + 1}
      for (i in 1:n) {
        if (names(Ybar)[i] == x[j]) {
          groups[i] <- z
        }
      }
    }
  }

  # Stop parallel
  stopCluster(cl)

  # Complete results
  breakgroups <- as.vector(read.table(breakgroupsext, header = FALSE, sep = ";")[,1])
  Bo <- round(as.vector(read.table(maxboext, header = FALSE)[,1]), 5)
  S2 <- round(as.vector(read.table(s2ext, header = FALSE)[,1]), 5)
  nutest <- round(as.vector(read.table(nutestext, header = FALSE)[,1]), 5)
  stattest <- round(as.vector(read.table(stattestext, header = FALSE)[,1]), 5)
  pvalues <- round(as.vector(read.table(pvalueext, header = FALSE)[,1]), 5)
  # Remove files
  file.remove(c(breakgroupsext, maxboext, s2ext, nutestext, stattestext, pvalueext))

  # Details of the results (invible)
  detres <- data.frame(Groups = breakgroups,
                       Bo = Bo,
                       Variance = S2,
                       DF = nutest,
                       Test = stattest,
                       "P-value" = pvalues)
  colnames(detres) <- c(gettext("Groups", domain = "R-MCP"),
                        "Bo",
                        gettext("Variance", domain = "R-MCP"),
                        gettext("DF", domain = "R-MCP"),
                        gettext("Test", domain = "R-MCP"),
                        gettext("P-value", domain = "R-MCP"))

  # Simple results
  result <- cbind(Ybar, groups)
  simple_results <- group.test2(result)

  # Output
  complete_results <- list("Details of results" = detres, "Simple results" = simple_results)

  return(complete_results)
}
#
# result <- sktest_bhering2008(y = y, trt = trt, dferror = dferror, mserror = mserror, replication = replication, alpha = alpha)
# result
