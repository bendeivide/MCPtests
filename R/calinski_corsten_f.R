# References
# [RAMOS,P.S.;VIEIRA,M.T.]. Bootstrap multiple comparison procedure based 
#                              on the F distribution.[ARTICLE].2014

# Exemplo do Caliski & Corsten
y <- Ybar <- c(97.7, 100.7, 111.3, 120.7, 124.3,
          128.7, 129.0, 131.0, 132, 141.7,
          150.7, 152.7, 176)
trt <- as.factor(1:13)
#length(trt)
mserror <- 124.29; dferror = 24; replication = 3


calinski_corsten_f <- function(y, trt, dferror, mserror, replication, alpha) {
  # Ordered means
  Ybar <- sort(tapply(y, trt, mean), decreasing = FALSE)
  
  # Length of means
  n <- length(Ybar)
  
  # Sum of squares between of treatments
  ssbt <- function(Ybar, replication, n) {
    r <- replication
    S <- matrix(0, n, n)
    for (i in 1:(n - 1)) {
      seq <- (i + 1):n
      for (j in seq) {
        g <- Ybar[i:j]
        S[i,j] <- r * (j - i) * var(g)
      }
    }
    return(S)  
  }
  # Results of the sum of squares between of treatments
  ssbt.result <- ssbt(Ybar, replication, n)
  
  # Position of the treatments
  pos <- 1:n
  
  # Result matrix
  complete.result <- foreach::foreach(p = 1:(n - 1), .combine = "rbind", .packages = "foreach") %do% {
    result <- stats::kmeans(Ybar, p, algorithm = "Hartigan-Wong", nstart = 10)$cluster
    Sp <- foreach::foreach(i = 1:p, .combine = '+') %do% {
      ssbt.result[min(pos[result == i]), max(pos[result == i])]
    }
    Sp <- round(Sp, 3)
    F_est <- round(Sp / (mserror * (n - 1)), 3)
    rbind(c(result, # Separation of groups
            Sp, # The smallest sum of square in p groups
            F_est, # F statistics
            round(pf(F_est, n-1, dferror, lower.tail = FALSE), 5) # P-value
    ))
  }

  # Counter
  i <- 1
  # Formation of the groups
  while(complete.result[i, n + 3] < alpha){
    # Results of Ramos-Vieira's test
    groups <- complete.result[i + 1, 1:n]
    i <- i + 1
  }
  
  # Simple result
  sresult <- data.frame(Ybar, groups)
  sresult <- sresult[order(Ybar, decreasing = c(TRUE, FALSE)),]
  simple_results <- group.test2(sresult)
  rownames(simple_results) <- row.names(sresult)
  
  # Complete result
  
  
}