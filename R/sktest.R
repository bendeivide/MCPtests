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
      cat(pvalue,"\n", file = "pvalues", append = TRUE)
      cat(substr(names(g1), 1, 3), "_vs_", substr(names(g2), 1, 3), ";", "\n", sep = " ", file = "breakgroups", append = TRUE)
      cat(max(b0),"\n", file = "maxbo", append = TRUE)
      cat(sig2,"\n", file = "s2", append = TRUE)
      cat(nu,"\n", file = "nutest", append = TRUE)
      cat(ts,"\n", file = "stattest", append = TRUE)
    }
    if (pvalue <= alpha) {
      # Classification of Scott-Knott's test
      for (i in 1:length(g1)) {
        cat(names(g1[i]),"\n", file = "results",append = TRUE)
      }
      cat("*","\n", file = "results", append = TRUE)
      # Auxiliar results
      cat(pvalue,"\n", file = "pvalues", append = TRUE)
      cat(substr(names(g1), 1, 3), "_vs_", substr(names(g2), 1, 3), ";", "\n", sep = " ", file = "breakgroups", append = TRUE)
      cat(max(b0),"\n", file = "maxbo", append = TRUE)
      cat(sig2,"\n", file = "s2", append = TRUE)
      cat(nu,"\n", file = "nutest", append = TRUE)
      cat(ts,"\n", file = "stattest", append = TRUE)
    }
    if (length(g1) > 1) sg(g1)
    if (length(g2) > 1) sg(g2)
  }

  # Loading the separation of the groups and generating external files
  sg(Ybar)

  # Result of Scott-Knott's test
  if (file.exists("results") == FALSE) {
    stop("Missing data entry!", call. = FALSE)
  } else{
    # Loading external file of results
    xx <- read.table("results")
    # Remove external file
    file.remove("results")
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
  breakgroups <- as.vector(read.table("breakgroups", header = FALSE, sep = ";")[,1])
  Bo <- round(as.vector(read.table("maxbo", header = FALSE)), 5)
  S2 <- round(as.vector(read.table("s2", header = FALSE)), 5)
  nutest <- round(as.vector(read.table("nutest", header = FALSE)), 5)
  stattest <- round(as.vector(read.table("stattest", header = FALSE)), 5)
  pvalues <- round(as.vector(read.table("pvalues", header = FALSE)), 5)
  # Remove files
  file.remove(c("breakgroups", "maxbo", "s2", "nutest", "stattest", "pvalues"))

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
