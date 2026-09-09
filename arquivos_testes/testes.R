# Example
#########
# Response variable
y <- rv <- c(100.08, 105.66, 97.64, 100.11, 102.60, 121.29, 100.80,
             99.11, 104.43, 122.18, 119.49, 124.37, 123.19, 134.16,
             125.67, 128.88, 148.07, 134.27, 151.53, 127.31)

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
replication <- 4
alpha <- 0.05

resultado <- MCPtests:::sktest(y, trt, dferror, mserror, replication, alpha,
                  parallel = FALSE)

resultadoMCPtest(y = rv,
        trt = treat,
        dferror = DFerror,
        mserror = MSerror,
        alpha = 0.05,
        main = "Multiple Comparison Procedure: MGM test",
        MCP = c("all"))

