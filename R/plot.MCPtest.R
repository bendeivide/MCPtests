plot.MCPtest <- function(result) {
  # 1. Recriar o data.frame (complete_results)
  # complete_results <- data.frame(
  #   Comparison = result$Comparison,
  #   Difference = ,
  #   Std.Error = rep(2.978028, 10),
  #   CI_lower = c(-3.1470460, -0.9393209, -10.7972777, -5.0594454, -14.9174022, -17.1251273, -15.0220423, -24.8799991, -27.0877243, -22.9675998),
  #   CI_upper = c(22.862960, 25.070685, 15.212728, 20.950560, 11.092604, 8.884878, 10.987963, 1.130007, -1.077719, 3.042406)
  # )
  complete_results <- result[[2]][[1]][[1]]

  # 2. Criar o vetor de cores
  # Se CI_lower for maior que 0 ou CI_upper for menor que 0, o intervalo nao contem zero, entao a cor e vermelha
  # Caso contrario, o intervalo contem zero e a cor e preta
  cores <- ifelse(complete_results$CI_lower > 0 | complete_results$CI_upper < 0,
                  "red", "black")

  # 3. Definir a posicao de cada linha no eixo Y
  # O primeiro item deve ficar no topo, entao invertemos a ordem
  n <- nrow(complete_results)
  y_pos <- n:1

  # 4. Configurar as margens para dar espaco aos rotulos
  par(mar = c(5, 6, 4, 2))

  # 5. Criar o grafico vazio com os limites corretos
  plot(x = complete_results$Difference,
       y = y_pos,
       type = "n",
       xlim = range(c(complete_results$CI_lower, complete_results$CI_upper)),
       ylim = c(0.5, n + 0.5),
       yaxt = "n",
       xlab = "Linear Function",
       ylab = "",
       main = "95% family-wise confidence level",
       cex.main = 1.3,
       font.main = 2)

  # 6. Adicionar a linha vertical tracejada no zero
  abline(v = 0, lty = 2, col = "black")

  # 7. Adicionar os intervalos de confianca como segmentos horizontais
  segments(x0 = complete_results$CI_lower,
           y0 = y_pos,
           x1 = complete_results$CI_upper,
           y1 = y_pos,
           col = cores,
           lwd = 1.5)

  # 8. Adicionar as barras verticais nas extremidades dos intervalos
  h <- 0.15
  segments(x0 = complete_results$CI_lower, y0 = y_pos - h,
           x1 = complete_results$CI_lower, y1 = y_pos + h,
           col = cores, lwd = 1.5)
  segments(x0 = complete_results$CI_upper, y0 = y_pos - h,
           x1 = complete_results$CI_upper, y1 = y_pos + h,
           col = cores, lwd = 1.5)

  # 9. Adicionar os pontos centrais com a diferenca estimada
  points(x = complete_results$Difference,
         y = y_pos,
         pch = 19,
         col = cores,
         cex = 1.2)

  # 10. Adicionar os rotulos do eixo Y com os nomes das comparacoes
  lab <- gsub(" == 0", "", rownames(complete_results))
  axis(side = 2, at = y_pos, labels = lab, las = 1)

  # 11. Adicionar linhas de grade horizontais leves
  abline(h = y_pos, col = "gray90", lty = 1)

  # 12. Replotar os elementos por cima das linhas de grade
  segments(x0 = complete_results$CI_lower, y0 = y_pos,
           x1 = complete_results$CI_upper, y1 = y_pos,
           col = cores, lwd = 1.5)
  points(x = complete_results$Difference, y = y_pos, pch = 19, col = cores, cex = 1.2)
  abline(v = 0, lty = 2, col = "black")

}


