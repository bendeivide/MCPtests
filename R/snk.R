snk <- function(y, trt, dferror, mserror, alpha) {
  n <- length(unique(trt))
  data.na    <- subset(data.frame(y, trt), is.na(y) == FALSE)
  Ybar <- tapply(y, trt, mean)
  Ybar <- sort(Ybar)
  rn   <- tapply(data.na[, 1], data.na[, 2], "length")
  rh   <- 1/mean(1/rn)
  std  <- sqrt(mserror / rh)

  # DIFERENCA 1: valor critico inicial depende de "range" (p = n)
  dms <- qtukey(1 - alpha, n, dferror) * std

  range <- n
  pos <- 1
  col <- 0
  qobs <- Ybar[n] - Ybar[1]
  aux <- rep(0, times = n)
  groups <- Ybar

  # Quando todas as medias sao estatisticamente iguais
  if ((qobs >= - dms) & (qobs <= dms)) aux[1:n] <- 1

  if (!(any(aux == 0))) {
    groups <- cbind(Ybar, aux)
  }

  if (any(aux == 0)) {
    continua1 = TRUE
  } else {
    continua1 = FALSE
  }

  if (continua1 == TRUE) {
    range <- range - 1
    pos <- 0
    ncomp <- n - range + 1
    ct <- 1
    if (range < 2) {
      continua2 <- FALSE
    } else {
      continua2 <- TRUE
    }

    if (continua2 == TRUE) {
      repeat {
        pos <- pos + 1
        qobs <- Ybar[pos + range - 1] - Ybar[pos]
        aux[1:n] <- 0

        # DIFERENCA 2: dms recalculado a cada passo com o "range" atual
        # No SNK, para uma amplitude com "range" medias, usa-se qtukey(1-alpha, range, dferror)
        dms <- qtukey(1 - alpha, range, dferror) * std

        if (((qobs >= - dms) & (qobs <= dms))) {
          aux[pos:(pos + range - 1)] <- 1
        }

        dentro <- FALSE
        if ((col > 0) & any(aux == 1)) {
          for (i in 1:col) {
            if (any(groups[aux == 1, i + 1] == 0)) {
              dentro <- FALSE
            } else {
              dentro <- TRUE
            }
            if (dentro == TRUE) break
          }
        }

        if ((dentro == FALSE) & any(aux == 1)) {
          groups <- cbind(groups, aux)
          col <- col + 1
        }

        ct <- ct + 1
        if (ct > ncomp) {
          range <- range - 1
          pos <- 0
          ncomp <- n - range + 1
          ct <- 1
        }

        if (range < 2) {
          continua2 <- FALSE
        } else {
          continua2 <- TRUE
        }
        if (continua2 == FALSE) break
      }
    }
  }

  # Detalhes dos resultados
  stdgeral <- sqrt(mserror / rh)

  # DIFERENCA 3: para comparacoes par-a-par no SNK, p = 2
  dms_pair <- qtukey(1 - alpha, 2, dferror) * stdgeral

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

  # DIFERENCA 4: nmeans = 2 (e nao 5, como estava hardcoded)
  pvalue <- ptukey(
    q = abs(diferenca) / Std.Error,
    nmeans = 2,
    df = dferror,
    lower.tail = FALSE
  )

  sigf <- signif_code(pvalue)

  det_results <- data.frame(
    Difference = diferenca,
    Std.Error  = Std.Error,
    CI_lower   = diferenca - dms_pair,
    CI_upper   = diferenca + dms_pair,
    `q value`  = diferenca / Std.Error,
    `Pr(>|q|)` = pvalue,
    ` `        = sigf,
    row.names  = nomes,
    check.names = FALSE
  )

  # Resultados simples (letras)
  result <- ProcTest(groups)
  simple_results <- group.test(result)

  output <- list(res1 = det_results,
                 res2 = simple_results)
  names(output) <- c(gettext("Details of results", domain = "R-MCP"),
                     gettext("Simple results",  domain = "R-MCP"))

  return(output)
}
