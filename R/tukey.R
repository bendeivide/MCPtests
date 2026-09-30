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
