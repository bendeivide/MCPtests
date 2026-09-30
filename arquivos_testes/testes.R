# Example
#########
# Response variable
y <- rv <- c(100.08, 105.66, 97.64, 100.11, 102.60, 121.29, 100.80,
             99.11, 104.43, 122.18, 119.49, 124.37, 123.19, 134.16,
             125.67, 128.88, 148.07, 134.27, 151.53, 127.31)
y <- rv <- rnorm(20, 100, 2)

# Treatments
trt <- treat <- factor(rep(LETTERS[1:5], each = 4))
#trt <- factor(rep(c("boi", "vaca", "bode", "cabra", "pato"), each = 4))

#ExpDes::crd(trt, y, quali = TRUE, mcomp = "sk")


# dados <- data.frame(trt, y)
# write.table(dados, "dados.csv", sep = ";")

# Anova
res     <- anova(aov(rv~treat))
dferror <- DFerror <- res$Df[2]
mserror <- MSerror <- res$`Mean Sq`[2]
replication <- n <- 4
alpha <- 0.05

#amostra <- SimulateData(5, 4, cenario = 2)
# Usar a base de dados de /data/
#load("/media/ben10/Backup/BEN_R/pkgs_published_22.09.2026/MCPtests_22.09.2026/MCPtests/data/dic.rda")
#save(amostra, file = "./data/dic.rda")
library(MCPtests)
data("dic")
y <- c(amostra$y)
trt <- as.factor(amostra$trat)
# MSerror
mserror <- summary(aov(amostra$y ~ amostra$trat))[[1]][[3]][[2]]

# SSerror
sserror <- summary(aov(amostra$y ~ amostra$trat))[[1]][[2]][[2]]

# DFerror
dferror <- summary(aov(amostra$y ~ amostra$trat))[[1]][[1]][[2]]
n <- length(unique(trt))
replication <- 4

library(MCPtests)
resultado <- MCPtest(y, trt, dferror, mserror, alpha, MCP = "tukey")
# Precisamos padronizar esta funcao
MCPtests:::plot.MCPtest(resultado)

resultado <- MCPtests:::tukey(y, trt, dferror, mserror, alpha)
snk(y, trt, replication, dferror, mserror, alpha)

# Testando com o ExpDes
ExpDes::tukey(y, trt, dferror, sserror, alpha = 0.05, group = TRUE,
  main = NULL)
ExpDes::crd(trt, y, quali = TRUE, mcomp = "snk", sigT = 0.05)

# Multiple comparison procedure: MGR test
MCPtest(aov(y ~trt), trt = "trt", alpha = 0.05,
        main = "Multiple Comparison Procedure: MGR test",
        MCP = c("tukey"))

# Comparando com o agricolae
library(agricolae)
HSD.test(y, trt, dferror, mserror, alpha = 0.05, console = TRUE)

# Usando o pacote multcomp
library(multcomp)
## Sumarizacao
plot(summary(glht(
  aov(y ~trt),
  linfct = mcp(trt = "Tukey")
)))
plot(confint(glht(
  aov(y ~trt),
  linfct = mcp(trt = "Tukey")
)))

