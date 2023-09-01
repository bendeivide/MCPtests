# References:
# Scott, A. J.; Knott, M. A cluster analysis method for grouping
#      means in the analyis of variance, v.30, n.3, 1974, p.507-512.

# Scott-Knott's test
sktest <- function(y, trt, dferror, mserror, replication, alpha,
                   parallel) {
  if (parallel) {
    # Parallel activate
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
  }

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

  # Function for the separation of the groups
  sg <- function(means, ...) {
    tm <- length(means) - 1
    # Calculate b0
    if (parallel) {
      b0 <- foreach(i = 1:tm, .combine = c) %dopar% {
        g1 <- means[1:i]
        g2 <- means[(i + 1):length(means)]
        gall <- c(g1, g2)
        sum(g1)^2/length(g1) + sum(g2)^2/length(g2) - (sum(gall))^2/length(gall)
      }
    } else {
      b0 <- foreach(i = 1:tm, .combine = c) %do% {
        g1 <- means[1:i]
        g2 <- means[(i + 1):length(means)]
        gall <- c(g1, g2)
        sum(g1)^2/length(g1) + sum(g2)^2/length(g2) - (sum(gall))^2/length(gall)
      }
    }

    # break of groups
    corte <- which.max(b0)

    # two groups
    g1 <- means[1:corte]
    g2 <- means[(corte + 1):length(means)]
    tg <- c(g1,g2)
    # Auxiliar name groups



    # Calculate ML estimate of sigma
    sig2 <- (1 / (length(tg) + dferror)) * (sum(tg^2) - (sum(tg))^2 / length(tg) + dferror * mserror / replication)

    # Test statistic
    ts <- pi / (2 * (pi - 2)) * max(b0) / sig2

    # Degrees of freedom
    nu <- length(tg) / (pi - 2)

    # P-value
    pvalue <- pchisq(ts, nu, lower.tail = FALSE)


    # Separation of the groups
    if (pvalue > alpha) {

      cat(pvalue,"\n", file = pvalueext, append = TRUE)
      cat(substr(names(g1), 1, 3), "_vs_", substr(names(g2), 1, 3), ";", "\n", sep = " ", file = breakgroupsext, append = TRUE)
      cat(max(b0),"\n", file = maxboext, append = TRUE)
      cat(sig2,"\n", file = s2ext, append = TRUE)
      cat(nu,"\n", file = nutestext, append = TRUE)
      cat(ts,"\n", file = stattestext, append = TRUE)
    }
    if (pvalue <= alpha) {
      # Classification of Scott-Knott's test
      for (i in 1:length(g1)) {
        cat(names(g1[i]),"\n", file = resultsext, append = TRUE)
      }
      cat("*","\n", file = resultsext, append = TRUE)
      # Auxiliar results
      cat(pvalue,"\n", file = pvalueext, append = TRUE)
      cat(substr(names(g1), 1, 3), "_vs_", substr(names(g2), 1, 3), ";", "\n", sep = " ", file = breakgroupsext, append = TRUE)
      cat(max(b0),"\n", file = maxboext, append = TRUE)
      cat(sig2,"\n", file = s2ext, append = TRUE)
      cat(nu,"\n", file = nutestext, append = TRUE)
      cat(ts,"\n", file = stattestext, append = TRUE)
    }
    if (length(g1) > 1) Recall(g1)
    if (length(g2) > 1) Recall(g2)
  }

  # Loading the separation of the groups and generating external files
  sg(Ybar)

  # Result of Scott-Knott's test
  if (file.exists(resultsext) == FALSE) {
    stop("Missing data entry!", call. = FALSE)
  } else{
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
  if (parallel) {
    # Stop parallel
    parallel::stopCluster(cl)
  }

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
  complete_results <- list("Details of results" = detres,
                           "Simple results" = simple_results)
  names(complete_results) <- c(gettext("Details of results", domain = "R-MCP"),
                               gettext("Simple results", domain = "R-MCP")
                               )

  return(complete_results)
}
