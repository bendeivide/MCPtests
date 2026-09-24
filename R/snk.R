snk <- function(y, trt, n, dferror, mserror, alpha)
{
  # Number of treatments
  k <- length(unique(trt))

  # Mean of each treatment
  Ybar <- tapply(y, trt, mean, na.rm = TRUE)

  # Order the means
  Ybar <- sort(Ybar)

  # Number of observations per treatment
  rh <- table(trt)[names(Ybar)]


  # ==========================================================
  # SNK grouping algorithm
  # ==========================================================

  # Initialize groups
  groups <- rep(0, times = k)

  # First group
  groups[1] <- 1

  ng <- 1

  # Compare treatments
  for (i in 2:k) {

    # Number of means included in the comparison
    r <- i

    # Studentized range critical value
    Rq <- qtukey(
      1 - alpha,
      nmeans = r,
      df = dferror
    )

    # Tukey-Kramer / SNK critical difference
    dms <- Rq * sqrt(
      mserror / 2 *
        (1 / rh[1] + 1 / rh[i])
    )

    # Difference between means
    qobs <- Ybar[i] - Ybar[1]

    if (qobs <= dms) {

      groups[i] <- groups[1]

    } else {

      groups[i] <- max(groups) + 1
      ng <- ng + 1
    }
  }


  # ==========================================================
  # Check for overlapping groups
  # ==========================================================

  if (any(groups > 1)) {

    posI <- 1
    fim <- FALSE

    repeat {

      ini <- groups[posI]

      posF <- max(which(groups == ini))

      if ((posF - posI) > 0) {

        for (i in (posI + 1):posF) {

          # Number of means included in the comparison
          r <- i - posI + 1

          # Studentized range critical value
          Rq <- qtukey(
            1 - alpha,
            nmeans = r,
            df = dferror
          )

          # SNK critical difference
          dms <- Rq * sqrt(
            mserror / 2 *
              (1 / rh[posI] + 1 / rh[i])
          )

          # Difference between means
          qobs <- Ybar[i] - Ybar[posI]

          if (qobs > dms) {

            groups[i:posF] <- max(groups) + 1
            ng <- ng + 1

            break
          }
        }
      }

      posI <- posF + 1

      if (posI > k)
        fim <- TRUE

      if (fim)
        break
    }
  }


  # ==========================================================
  # Simple results
  # ==========================================================

  result <- cbind(Ybar, groups)

  simple_results <- group.test(result)


  # ==========================================================
  # Complete results
  # ==========================================================

  # Number of pairwise comparisons
  ncomp <- choose(k, 2)

  # Initialize vectors
  treatment1 <- character(ncomp)
  treatment2 <- character(ncomp)

  mean1 <- numeric(ncomp)
  mean2 <- numeric(ncomp)

  difference <- numeric(ncomp)
  r <- numeric(ncomp)

  se <- numeric(ncomp)
  qvalue <- numeric(ncomp)
  qcritical <- numeric(ncomp)

  pvalue <- numeric(ncomp)
  dms <- numeric(ncomp)

  lower <- numeric(ncomp)
  upper <- numeric(ncomp)

  significant <- logical(ncomp)

  # Counter
  cont <- 1


  # ==========================================================
  # Pairwise comparisons
  # ==========================================================

  for (i in 1:(k - 1)) {

    for (j in (i + 1):k) {

      # Treatment names
      treatment1[cont] <- names(Ybar)[i]
      treatment2[cont] <- names(Ybar)[j]

      # Means
      mean1[cont] <- Ybar[i]
      mean2[cont] <- Ybar[j]

      # Difference
      difference[cont] <- Ybar[j] - Ybar[i]

      # Number of means included in the SNK comparison
      r[cont] <- j - i + 1

      # Standard error
      se[cont] <- sqrt(
        mserror / 2 *
          (1 / rh[i] + 1 / rh[j])
      )

      # Observed studentized range statistic
      qvalue[cont] <- abs(difference[cont]) / se[cont]

      # Critical studentized range
      qcritical[cont] <- qtukey(
        1 - alpha,
        nmeans = r[cont],
        df = dferror
      )

      # SNK critical difference
      dms[cont] <- qcritical[cont] * se[cont]

      # Adjusted p-value
      pvalue[cont] <- ptukey(
        qvalue[cont],
        nmeans = r[cont],
        df = dferror,
        lower.tail = FALSE
      )

      # Confidence interval
      lower[cont] <- difference[cont] - dms[cont]
      upper[cont] <- difference[cont] + dms[cont]

      # Significance
      significant[cont] <- qvalue[cont] > qcritical[cont]

      cont <- cont + 1
    }
  }


  # ==========================================================
  # Complete results data frame
  # ==========================================================

  cresult <- data.frame(
    Treatment1 = treatment1,
    Treatment2 = treatment2,
    Mean1 = mean1,
    Mean2 = mean2,
    Difference = difference,
    Range = r,
    SE = se,
    q = qvalue,
    q_critical = qcritical,
    p_value = pvalue,
    DMS = dms,
    CI_lower = lower,
    CI_upper = upper,
    Significant = significant,
    stringsAsFactors = FALSE
  )


  # ==========================================================
  # Output
  # ==========================================================

  complete_results <- list(
    "Details of results" = cresult,
    "Simple results" = simple_results
  )

  return(complete_results)
}
