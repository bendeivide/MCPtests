
# tukey <- function(y, trt, n, dferror, mserror, alpha)
# {
#   # Number of treatments
#   k <- length(unique(trt))
#
#   # Studentized range critical value
#   Rq <- qtukey(1 - alpha, k, dferror)
#
#   # Mean of each treatment
#   Ybar <- tapply(y, trt, mean, na.rm = TRUE)
#
#   # Order the means
#   Ybar <- sort(Ybar)
#
#   # Number of observations per treatment
#   rh <- table(trt)[names(Ybar)]
#
#   # Initialize groups
#   groups <- rep(0, times = k)
#
#   # First group
#   groups[1] <- 1
#   ng <- 1
#
#   # Compare each treatment with the first treatment
#   for (i in 2:k) {
#
#     # Tukey Kramer critical difference
#     dms <- Rq * sqrt(
#       mserror / 2 *
#         (1 / rh[1] + 1 / rh[i])
#     )
#
#     # Difference between means
#     qobs <- Ybar[i] - Ybar[1]
#
#     if (qobs <= dms) {
#       groups[i] <- groups[1]
#     } else {
#       groups[i] <- max(groups) + 1
#       ng <- ng + 1
#     }
#   }
#
#   # Check for overlapping groups
#   if (any(groups > 1)) {
#
#     posI <- 1
#     fim <- FALSE
#
#     repeat {
#
#       ini <- groups[posI]
#
#       posF <- max(which(groups == ini))
#
#       if ((posF - posI) > 0) {
#
#         for (i in (posI + 1):posF) {
#
#           # Tukey Kramer critical difference
#           dms <- Rq * sqrt(
#             mserror / 2 *
#               (1 / rh[posI] + 1 / rh[i])
#           )
#
#           # Difference between means
#           qobs <- Ybar[i] - Ybar[posI]
#
#           if (qobs > dms) {
#
#             groups[i:posF] <- max(groups) + 1
#             ng <- ng + 1
#
#             break
#           }
#         }
#       }
#
#       posI <- posF + 1
#
#       if (posI > k)
#         fim <- TRUE
#
#       if (fim)
#         break
#     }
#   }
#
#
#   # ==========================================================
#   # Simple results
#   # ==========================================================
#
#   result <- cbind(Ybar, groups)
#
#   simple_results <- group.test(result)
#
#
#   # ==========================================================
#   # Complete results
#   # ==========================================================
#
#   # Number of pairwise comparisons
#   ncomp <- choose(k, 2)
#
#   # Initialize result vectors
#   treatment1 <- character(ncomp)
#   treatment2 <- character(ncomp)
#
#   mean1 <- numeric(ncomp)
#   mean2 <- numeric(ncomp)
#   difference <- numeric(ncomp)
#   se <- numeric(ncomp)
#   qvalue <- numeric(ncomp)
#   pvalue <- numeric(ncomp)
#   dms <- numeric(ncomp)
#   lower <- numeric(ncomp)
#   upper <- numeric(ncomp)
#   significant <- logical(ncomp)
#
#   # Counter
#   cont <- 1
#
#   # Pairwise comparisons
#   for (i in 1:(k - 1)) {
#
#     for (j in (i + 1):k) {
#
#       # Treatment names
#       treatment1[cont] <- names(Ybar)[i]
#       treatment2[cont] <- names(Ybar)[j]
#
#       # Treatment means
#       mean1[cont] <- Ybar[i]
#       mean2[cont] <- Ybar[j]
#
#       # Difference between means
#       difference[cont] <- Ybar[j] - Ybar[i]
#
#       # Standard error
#       se[cont] <- sqrt(
#         mserror / 2 *
#           (1 / rh[i] + 1 / rh[j])
#       )
#
#       # Observed studentized range statistic
#       qvalue[cont] <- abs(difference[cont]) / se[cont]
#
#       # Tukey-Kramer adjusted p-value
#       pvalue[cont] <- ptukey(
#         qvalue[cont],
#         nmeans = k,
#         df = dferror,
#         lower.tail = FALSE
#       )
#
#       # Tukey-Kramer critical difference
#       dms[cont] <- Rq * se[cont]
#
#       # Confidence interval
#       lower[cont] <- difference[cont] - dms[cont]
#       upper[cont] <- difference[cont] + dms[cont]
#
#       # Significance
#       significant[cont] <- pvalue[cont] <= alpha
#
#       cont <- cont + 1
#     }
#   }
#
#
#   # Complete results data frame
#   cresult <- data.frame(
#     Treatment1 = treatment1,
#     Treatment2 = treatment2,
#     Mean1 = mean1,
#     Mean2 = mean2,
#     Difference = difference,
#     SE = se,
#     q = qvalue,
#     p_value = pvalue,
#     DMS = dms,
#     CI_lower = lower,
#     CI_upper = upper,
#     Significant = significant,
#     stringsAsFactors = FALSE
#   )
#
#
#   # ==========================================================
#   # Output
#   # ==========================================================
#
#   complete_results <- list(
#     "Details of results" = cresult,
#     "Simple results" = simple_results
#   )
#
#   return(complete_results)
# }
#
# tukey <- function(y, trt, dferror, mserror, alpha)
# {
#   # ==========================================================
#   # Number of treatments
#   # ==========================================================
#
#   k <- length(unique(trt))
#
#   if (k < 2)
#     stop("At least two treatments are required.")
#
#
#   # ==========================================================
#   # Treatment means
#   # ==========================================================
#
#   Ybar <- tapply(
#     y,
#     trt,
#     mean,
#     na.rm = TRUE
#   )
#
#   # Order means from largest to smallest
#   Ybar <- sort(
#     Ybar,
#     decreasing = TRUE
#   )
#
#   # Number of observations per treatment
#   rh <- table(trt)[names(Ybar)]
#
#
#   # ==========================================================
#   # Studentized range critical value
#   # ==========================================================
#
#   Rq <- qtukey(
#     1 - alpha,
#     k,
#     dferror
#   )
#
#
#   # ==========================================================
#   # Pairwise comparisons
#   # ==========================================================
#
#   ncomp <- choose(k, 2)
#
#   treatment1 <- character(ncomp)
#   treatment2 <- character(ncomp)
#
#   mean1 <- numeric(ncomp)
#   mean2 <- numeric(ncomp)
#
#   difference <- numeric(ncomp)
#   se <- numeric(ncomp)
#   qvalue <- numeric(ncomp)
#   pvalue <- numeric(ncomp)
#
#   dms <- numeric(ncomp)
#   lower <- numeric(ncomp)
#   upper <- numeric(ncomp)
#
#   significant <- logical(ncomp)
#
#   # Matrix of nonsignificant comparisons
#   #
#   # TRUE  = not significantly different
#   # FALSE = significantly different
#   #
#   nonsig <- matrix(
#     TRUE,
#     nrow = k,
#     ncol = k,
#     dimnames = list(
#       names(Ybar),
#       names(Ybar)
#     )
#   )
#
#   cont <- 1
#
#
#   for (i in 1:(k - 1)) {
#
#     for (j in (i + 1):k) {
#
#       # Treatment names
#       treatment1[cont] <- names(Ybar)[i]
#       treatment2[cont] <- names(Ybar)[j]
#
#       # Means
#       mean1[cont] <- Ybar[i]
#       mean2[cont] <- Ybar[j]
#
#       # Difference
#       difference[cont] <-
#         Ybar[i] - Ybar[j]
#
#       # Standard error
#       se[cont] <- sqrt(
#         mserror / 2 *
#           (1 / rh[i] + 1 / rh[j])
#       )
#
#       # Observed studentized range
#       qvalue[cont] <-
#         abs(difference[cont]) /
#         se[cont]
#
#       # Tukey adjusted p-value
#       pvalue[cont] <- ptukey(
#         qvalue[cont],
#         nmeans = k,
#         df = dferror,
#         lower.tail = FALSE
#       )
#
#       # Tukey-Kramer critical difference
#       dms[cont] <-
#         Rq * se[cont]
#
#       # Confidence interval
#       lower[cont] <-
#         difference[cont] - dms[cont]
#
#       upper[cont] <-
#         difference[cont] + dms[cont]
#
#       # Significance
#       significant[cont] <-
#         pvalue[cont] <= alpha
#
#       # Store comparison
#       if (significant[cont]) {
#
#         nonsig[i, j] <- FALSE
#         nonsig[j, i] <- FALSE
#       }
#
#       cont <- cont + 1
#     }
#   }
#
#
#   # ==========================================================
#   # Complete results
#   # ==========================================================
#
#   cresult <- data.frame(
#     Treatment1 = treatment1,
#     Treatment2 = treatment2,
#     Mean1 = mean1,
#     Mean2 = mean2,
#     Difference = difference,
#     SE = se,
#     q = qvalue,
#     p_value = pvalue,
#     DMS = dms,
#     CI_lower = lower,
#     CI_upper = upper,
#     Significant = significant,
#     stringsAsFactors = FALSE
#   )
#
#
#   # ==========================================================
#   # Build groups
#   # ==========================================================
#   #
#   # We construct groups using the complete matrix of
#   # nonsignificant comparisons.
#   #
#   # A group can contain treatments only when every pair
#   # within that group is nonsignificant.
#   #
#   # ==========================================================
#
#
#   # ----------------------------------------------------------
#   # Start with one group for each treatment
#   # ----------------------------------------------------------
#
#   group_list <- lapply(
#     seq_len(k),
#     function(i) i
#   )
#
#
#   # ----------------------------------------------------------
#   # Try to combine groups
#   # ----------------------------------------------------------
#
#   changed <- TRUE
#
#   while (changed) {
#
#     changed <- FALSE
#
#     if (length(group_list) > 1) {
#
#       for (a in seq_len(length(group_list) - 1)) {
#
#         for (b in (a + 1):length(group_list)) {
#
#           g1 <- as.integer(group_list[[a]])
#           g2 <- as.integer(group_list[[b]])
#
#           candidate <- sort(
#             unique(c(g1, g2))
#           )
#
#           # Check whether every pair in the candidate
#           # is nonsignificant
#           valid <- TRUE
#
#           if (length(candidate) > 1) {
#
#             for (u in 1:(length(candidate) - 1)) {
#
#               for (v in (u + 1):length(candidate)) {
#
#                 if (!nonsig[
#                   candidate[u],
#                   candidate[v]
#                 ]) {
#
#                   valid <- FALSE
#                   break
#                 }
#               }
#
#               if (!valid)
#                 break
#             }
#           }
#
#           if (valid) {
#
#             group_list[[a]] <- candidate
#
#             group_list <- group_list[-b]
#
#             changed <- TRUE
#
#             break
#           }
#         }
#
#         if (changed)
#           break
#       }
#     }
#   }
#
#
#   # ----------------------------------------------------------
#   # Add missing groups
#   # ----------------------------------------------------------
#
#   represented <- sort(
#     unique(
#       unlist(group_list)
#     )
#   )
#
#   missing <- setdiff(
#     seq_len(k),
#     represented
#   )
#
#   if (length(missing) > 0) {
#
#     for (i in missing) {
#
#       group_list[[length(group_list) + 1]] <-
#         i
#     }
#   }
#
#
#   # ----------------------------------------------------------
#   # Order groups according to the first treatment they
#   # contain
#   # ----------------------------------------------------------
#
#   first_position <- sapply(
#     group_list,
#     min
#   )
#
#   group_list <- group_list[
#     order(first_position)
#   ]
#
#
#   # ==========================================================
#   # Membership matrix
#   # ==========================================================
#
#   membership <- matrix(
#     0,
#     nrow = k,
#     ncol = length(group_list)
#   )
#
#   for (j in seq_along(group_list)) {
#
#     membership[
#       group_list[[j]],
#       j
#     ] <- 1
#   }
#
#
#   # ==========================================================
#   # Remove redundant groups
#   # ==========================================================
#
#   if (ncol(membership) > 1) {
#
#     keep <- rep(
#       TRUE,
#       ncol(membership)
#     )
#
#     for (j in seq_len(ncol(membership))) {
#
#       members_j <- which(
#         membership[, j] == 1
#       )
#
#       if (length(members_j) <= 1)
#         next
#
#       for (l in seq_len(ncol(membership))) {
#
#         if (l == j)
#           next
#
#         members_l <- which(
#           membership[, l] == 1
#         )
#
#         if (all(members_j %in% members_l)) {
#
#           keep[j] <- FALSE
#           break
#         }
#       }
#     }
#
#     membership <- membership[
#       ,
#       keep,
#       drop = FALSE
#     ]
#   }
#
#
#   # ==========================================================
#   # Make sure first treatment belongs to first group
#   # ==========================================================
#
#   first_group <- which(
#     membership[1, ] == 1
#   )
#
#   if (length(first_group) > 0 &&
#       first_group[1] != 1) {
#
#     membership <- membership[
#       ,
#       c(
#         first_group[1],
#         setdiff(
#           seq_len(ncol(membership)),
#           first_group[1]
#         )
#       ),
#       drop = FALSE
#     ]
#   }
#
#
#   # ==========================================================
#   # Prepare object for group.test()
#   # ==========================================================
#
#   result <- cbind(
#     Ybar,
#     membership
#   )
#
#   colnames(result)[-1] <- paste0(
#     "g",
#     seq_len(ncol(membership))
#   )
#
#
#   # ==========================================================
#   # Simple results
#   # ==========================================================
#
#   simple_results <- group.test(
#     result
#   )
#
#
#   # ==========================================================
#   # Output
#   # ==========================================================
#
#   complete_results <- list(
#     "Details of results" = cresult,
#     "Simple results" = simple_results
#   )
#
#   return(complete_results)
# }

tukey <- function(y, trt, dferror, mserror, alpha) {
  n <- length(unique(trt))
  data.na    <- subset(data.frame(y, trt), is.na(y) == FALSE)
  Ybar <- tapply(y, trt, mean)
  Ybar <- sort(Ybar)
  rn   <- tapply(data.na[, 1], data.na[, 2], "length")
  rh   <- 1/mean(1/rn)
  dms <- qtukey(1 - alpha, n, dferror) * sqrt(mserror / rh)
  std <- sqrt(mserror / rh)
  range <- n
  pos <- 1
  col <- 0
  qobs <- Ybar[n] - Ybar[1]
  aux <- rep(0, times = n)
  groups <- Ybar
  # Quando todas as medias sao estatisticamente iguais
  # comparando o maior gap entre as medias primeiro
  if ((qobs >= - dms) & (qobs <= dms)) aux[1:n] <- 1
  # Caso alguma media seja diferente, cria-se
  # o objeto groups para auxiliar nas letras
  if (!(any(aux == 0))) {
    groups <- cbind(Ybar, aux)
  }
  # Objeto logico continua1 para auxiliar
  # na composicao das letras
  if (any(aux == 0)) {
    continua1 = TRUE
  } else {
    continua1 = FALSE
  }
  # Em caso positivo de continua1
  if (continua1 == TRUE) {
    # Como as medias nao sao iguais, entre
    # Ybar[range] - Ybar[1]
    # faz-se a comparacao das range - 1
    # medias
    range <- range - 1
    # Inicialmente na posicao 0
    pos <- 0
    # Verificando o numero de comparacoes envolvidas
    # de medias no gap range - 1
    ncomp <- n - range + 1
    # Contador para controlar o num de comparacoes
    ct <- 1
    # Em caso de haver apenas 1 media envolvida
    # na comparacao
    if (range < 2) {
      # Cria-se um outro operador logico alternativo
      continua2 <- FALSE
    } else {
      continua2 <- TRUE
    }
    # Considerando que as comparacoes
    # apresentam mais de uma media envolvida
    if (continua2 == TRUE) {
      repeat {
        # Posicao de comparacao
        pos <- pos + 1
        # qobs baseada na pos+range-1 primeiras
        # medias
        qobs <- Ybar[pos+range-1] - Ybar[pos]
        # coluns aux a ser inserida em groups
        # identificando por 1 que as medias iguais
        aux[1:n] <- 0

        if (((qobs >= - dms) & (qobs <= dms))) {
          # Identifica por 1, as pos+range-1 medias
          # envolvidas sao iguais
          aux[pos:(pos + range - 1)] <- 1
        }

        dentro <- FALSE
        # No primeiro passo do repeat,
        # este if eh despresado. A partir,
        # do passo 2, se existi,
        # Verifica-se entre as medias
        # Ybar[pos] a Ybar[pos + range - 1]
        # se aux[pos:(pos + range - 1)] inserido
        # no passo anterior, apresenta valor 0.
        # Isso significa que existem medias diferentes
        if ((col > 0) & any(aux == 1)) {
          for (i in 1:col)
          {
            if (any(groups[aux == 1, i + 1] == 0)) {
              dentro <- FALSE
            } else {
              dentro <- TRUE
            }
            if (dentro == TRUE) break
          }
        }
        # Atualizar o contador col e o objeto groups
        if ((dentro == FALSE) & any(aux == 1)) {
          groups <- cbind(groups, aux)
          col <- col + 1
        }
        # Atualizacao do contador
        ct <- ct + 1
        # Arualizacao dos contadores
        # para auxiliar nas comparacoes
        # das medias envolvidas
        if (ct > ncomp)
        {
          range <- range - 1
          pos <- 0
          ncomp <- n - range + 1
          ct <- 1
        }
        # Se range < 2 eh pq esta
        # comparacao contem apenas uma
        # media envolvida, logo
        # se finaliza o loop do repeat
        if (range < 2) {
          continua2 <- FALSE
        } else {
          continua2 <- TRUE
        }
        if (continua2 == FALSE) break
      }
    }
  }

  # Details of the results (invible)
  stdgeral <- sqrt(mserror / rh)
  dms <- qtukey(1 - alpha, n, dferror) * stdgeral

  # Tabela de ICs

  trat <- levels(trt)
  ntrat <- length(trat)
  ncomp <- choose(ntrat, 2)
  Std.Error <- rep(stdgeral, ncomp)



  # Funcao para captar o nome das
  # combinacoes dos tratamentos
  nomes <- unlist(
    lapply(2:ntrat, function(i) {
      paste0(trat[i], " - ", trat[1:(i - 1)], " == 0")
    })
  )

  # Funcao para captar a diferenca
  # das medias envolvidas
  diferenca <- unlist(
    lapply(2:ntrat, function(i) {
      Ybar[trat[i]] - Ybar[trat[1:(i - 1)]]
    })
  )

  # P-values
  pvalue <- ptukey(
    q = abs(diferenca) / Std.Error,
    nmeans = 5,
    df = dferror,
    lower.tail = FALSE
  )

  sigf <- signif_code(pvalue)

  det_results <- data.frame(
    Difference = diferenca,
    Std.Error = Std.Error,
    CI_lower = diferenca - dms,
    CI_upper = diferenca + dms,
    `q value` = diferenca / Std.Error,
    `Pr(>|q|)` = pvalue,
    ` `        = sigf,
    row.names = nomes,
    check.names = FALSE
  )


  # Simple results
  result <- ProcTest(groups)
  simple_results <- group.test(result)


  # Output
  output <- list(res1 = det_results,
                 res2 = simple_results)
  names(output) <- c(gettext("Details of results", domain = "R-MCP"),
                               gettext("Simple results", domain = "R-MCP"))

  return(output)
}
