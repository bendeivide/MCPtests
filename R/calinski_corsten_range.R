# References:
# [RAMOS, P.S.;FERREIRA,D.F.].Agrupamento de médias via bootstrap de populacoes 
# normais e nao-normais.[ARTICLE].2009

# [CALINSKI_CORSTEN].Clustering Means in ANOVA by Simultaneous Testing.[Article].1985

# Exemplo Steel & Torrie
# y <- c(19.4, 32.6, 27, 32.1, 33, 17.7,
#        24.8, 27.9, 25.2, 24.3, 17, 19.4, 9.1, 11.9, 15.8, 20.7, 21, 20.5, 18.8, 18.6, 14.3,
#        14.4, 11.8, 11.6, 14.2, 17.3, 19.4, 19.1, 16.9, 20.8)
# replication <- 5
# trt <- as.factor(rep(1:6, each = replication))
# dados <- data.frame(trt, y)
# anova(aov(y ~ trt))
# mserror <- 11.78867; dferror <- 28
# 
# 
# # Exemplo do Caliski & Corsten
# Ybar <- c(97.7, 100.7, 111.3, 120.7, 124.3, 
#           128.7, 129.0, 131.0, 132, 141.7, 
#           150.7, 152.7, 176)
# trt <- as.factor(letters[1:13])
# length(trt)
# mserror <- 124.29; dferror = 24; replication = 3
# 
# dados <- data.frame(y = Ybar, trt = 1:13)
# plot(cluster::agnes(dados))


calinski_corsten_range <- function(y, trt, dferror, mserror, replication, alpha) {
  # Parallel activate
  # cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # doParallel::registerDoParallel(cl)
  
  # Ordered means
  Ybar <- sort(tapply(y, trt, mean), decreasing = FALSE)
  
  # Length of means
  n <- length(Ybar)
  
  # Results of calinski-corsten's test
  groups <- rep(0, times = length(Ybar))
  
  # Matrix of distances (Ranges)
  dih.f <- function(Ybar) {
    k <- length(Ybar)
    D <- matrix(0, k, k)
    for (i in 1:(k - 1)) {
      seq = (i + 1):k
      for (j in seq)
      {
        D[i,j] = Ybar[j] - Ybar[i]
      }
    }
    return(D) 
  }
  
  # Grouping of means using the most distant neighbor
  # D <- Object of dih.f() function
  # n <- number of treatment means
  gdist = function(D, n){
    k <- n 
    v = matrix(1:k,k,1)
    res = as.data.frame(matrix(1:(k - 1), k - 1, 6))
    colnames(res) <- c("Step", "LL", "UL", "Range", "Stud. Range", "P-value")
    ct = 1
    
    # Search the minimum range
    fmin <- function(D, k) { 
      min <- D[1,2]; II = 1;JJ = 2
      for (ii in 1:(k - 1)) {
        for (jj in (ii + 1):k) {
          if (min > D[ii, jj]) {
            min = D[ii,jj]; II = ii; JJ = jj
          }
        }
      }  
      resmin <- list(min = min, II = II, JJ = JJ) 
      return(resmin)
    }
    
    # Auxiliar function
    arma.res <- function(ct, min, II, JJ) { 
      if ((v[II] > 0) & (v[JJ] > 0)) {
        res[ct,2] <- min(v[II],v[JJ])
        res[ct,3] <- max(v[II],v[JJ])
        res[ct,4] <- min }
      if ((v[II] < 0) & (v[JJ] > 0)) {
        L = res[abs(v[II]),2]
        U = res[abs(v[II]),3]
        min1 <- min(L,U); max1 <- max(L,U)
        res[ct,2] <- min(min1,v[JJ])
        res[ct,3] <- max(max1,v[JJ])
        res[ct,4] <- min }
      if ((v[II] > 0) & (v[JJ] < 0)) {
        L <- res[abs(v[JJ]),2]
        U <- res[abs(v[JJ]),3]
        min1 <- min(L,U); max1 <- max(L,U)
        res[ct,2] <- min(min1,v[II])
        res[ct,3] <- max(max1,v[II])
        res[ct,4] <- min }
      if ((v[II] < 0) & (v[JJ] < 0)) {
        L = res[abs(v[II]),2]
        U = res[abs(v[II]),3]
        min1 <- min(L,U); max1 <- max(L,U)
        L <- res[abs(v[JJ]),2]
        U <- res[abs(v[JJ]),3]
        min2 <- min(L,U); max2 <- max(L,U)
        res[ct,2] <- min(min1,min2)
        res[ct,3] <- max(max1,max2)
        res[ct,4] <- min }
      
      # Range Statistics
      R <- res[ct,4]
      
      # Studentized range
      Q <- R / ((mserror/replication)^0.5) 
      res[ct,5] <- round(Q, 5)
      
      # P-value of statistics
      res[ct,6] <- round(ptukey(Q, k, dferror, lower.tail = FALSE), 5)
      
      return(res) 
    }
    
    # Update of the matrix of distances
    monta.res <- function(D, k, II, JJ) {
      D1 <- matrix(0, k - 1, k - 1); iii = 2; jjj = iii + 1
      for (ii in 1:(k - 1)) {
        for (jj in (ii + 1):k) {
          if ((abs(ii - II) > 1e-6) & (abs(jj - JJ) > 1e-6) &
              (abs(ii - JJ) > 1e-6) & (abs(jj - II) > 1e-6)) {
            D1[iii,jjj] <- D[ii,jj]
            jjj <- jjj + 1
            if (jjj >= k) {
              iii = iii + 1; jjj = iii + 1 }
          }
        }
      }
      jj = 2;
      for (ii in 1:k) {
        if ((abs(ii - II) > 1e-6) & (abs(ii-JJ) > 1e-6)) {
          if ((ii > II) & (ii > JJ)) aux = max(D[II,ii], D[JJ,ii])
          if ((ii > II) & (ii < JJ)) aux = max(D[II,ii], D[jj,JJ])
          if ((ii < II) & (ii > JJ)) aux = max(D[ii,II], D[JJ,ii])
          if ((ii < II) & (ii < JJ)) aux = max(D[ii,II], D[ii,JJ])
          D1[1,jj] = aux
          jj = jj + 1 }
      } 
      return(D1) 
    } 
    
    # Create auxiliar vector
    monta.vet <- function(v, II, JJ) {
      v1 <- v[v != v[II]]
      v1 <- v1[v1 != v[JJ]]
      v1 <- c(-ct,v1)
      return(v1) 
    }
    
    # Grouping
    Dct <- D; kt <- k
    for (ct in 1:(k - 1)) {
      estmin <- fmin(Dct, kt)
      res <- arma.res(ct, estmin$min, estmin$II, estmin$JJ)
      Dct <- monta.res(Dct, kt, estmin$II, estmin$JJ)
      v <- monta.vet(v, estmin$II, estmin$JJ)
      kt <- kt - 1 
    }
    return(res) 
  }
  
  # Calculate
  D <- dih.f(Ybar)
  cresult <- gdist(D, n)
  
  # Assigning letters to groups
  z <- 1
  for (i in 1:dim(cresult)[1]) {
    if (cresult[i,6] > alpha) {
      groups[cresult[i,2]:cresult[i,3]] <- z
      z <- z + 1
    }
  }
  
  # Simple result
  sresult <- data.frame(Ybar, groups)
  sresult <- sresult[order(Ybar, decreasing = c(TRUE, FALSE)),]
  simple_results <- group.test2(sresult)
  rownames(simple_results) <- row.names(sresult)
  
  
  
  # Update cresult
  for (i in 1:dim(cresult)[1]) {
    cresult[i, 2] <- names(Ybar)[as.numeric(cresult[i, 2])]
    cresult[i, 3] <- names(Ybar)[as.numeric(cresult[i, 3])]
  }
  
  # # Stop parallel
  # stopCluster(cl)
  
  # Output
  complete_results <- list("Details of results" = cresult, "Simple results" = simple_results)
  return(complete_results)
}

#calinski_corsten_range(y, trt, dferror, mserror, replication, alpha = 0.05)
