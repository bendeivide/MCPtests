#' Multiple comparison procedures to the means.
#'
#' \code{MCPtest} applies several tests of multiple comparisons present
#'     in the literature. The tests chosen are based on the evaluation
#'     of the researcher to decide the choice of test for analysis
#'     in the experiment.
#' @param y Model (aov or lm), numeric vector containing the
#'     response variable or the mean of the treatments.
#' @param trt Constant (y = model) or a vector containing
#'     the treatments.
#' @param dferror Degrees of freedom of the Mean Square Error.
#' @param mserror Mean Square Error.
#' @param replication Number de repetitions of the treatments in the experiment.
#' For unbalanced data should be informed the harmonic mean of repetitions.
#' This argument should be informed only if \code{ismean = TRUE}.
#' @param alpha Significant level. The default is \code{alpha = 0.05}.
#' @param main Title of the analysis.
#' @param MCP Allows choosing the multiple comparison test;
#'     the \emph{defaut} is "all". This option will go perform all tests.
#'     However, the options are: the Skott-Knott midrange test ("MGM"),
#'     the Skott-Knott Range test ("MGR"), the Tukey midrange test ("TM"),
#'     the Scott-Knott's test ("SK").
#' @param ismean Logic. If \code{FALSE} (default), the \code{y} argument represents
#' a model (aov or lm) or a numeric vector containing the response variable.
#' If \code{TRUE} the \code{y} argument represents the mean of treatments.
#' @param parallel Logic. If \code{FALSE} (default), the \code{MCPtest} function
#' runs without parallelization. Otherwise, the code is executed with
#' parallelization. Note that the parallelization is not always more efficient.
#  @param console Logic. If \code{TRUE} (default), the \code{MCPtest} function
#  prints the results on the console even though the results are stored in an
#  object. If \code{FALSE}, otherside.
#' @return \code{MCPtest} returns the print of a list of results. First,
#'     the summary of \code{y}. Second, the statistics
#'     of the test chosen. And finally, the  mean group results for each test.
#'     If \code{MCPtest} function is stored
#'     in an object, the results will be printed and
#'     also stored in the object.
#' @details The \code{MCP} argument allows you to choose various tests
#'     of multiple comparisons at once. For example,
#'     \code{MCP = c("MGM", "MGR")}, and so on.
#' @references
#' BATISTA, Ben Deivide de Oliveira. Testes de comparacoes multiplas
#'     baseados na distribuicao da \emph{midrange} \emph{estudentizada}
#'     externamente. 2016. 194f. Tese (Doutorado em Estatistica e
#'     Experimentacao Agropecuaria) - Universidade Federal de Lavras, 2016.
#'     (Portuguese, Brazil)
#'
#' SCOTT, A. J.; KNOTT, M. A cluster analysis method for grouping means in
#'     the analysis of variance. Biometrics, International Biometric Society,
#'     v. 30, n. 3, p. 507-512, 1974.
#'
#' DUNCAN, D. B. Multiple range and multiple F Tests. Biometrics, International
#'     Biometric Society, v. 11, n. 1,p. 1-42, 1955.
#' @examples
#' # Simulated data (completely randomized design)
#'
#' # Response variable
#' rv <- c(100.08, 105.66, 97.64, 100.11, 102.60, 121.29, 100.80,
#'         99.11, 104.43, 122.18, 119.49, 124.37, 123.19, 134.16,
#'         125.67, 128.88, 148.07, 134.27, 151.53, 127.31)
#'
#' # Treatments
#' treat <- factor(rep(LETTERS[1:5], each = 4))
#'
#' # Anova
#' res     <- anova(aov(rv~treat))
#' DFerror <- res$Df[2]
#' MSerror <- res$`Mean Sq`[2]
#'
#' # Loading the MCPtests package
#' library(MCPtests)
#'
#' # applying the tests
#' results <- MCPtest(y = rv,
#'                   trt = treat,
#'                   dferror = DFerror,
#'                   mserror = MSerror,
#'                   alpha = 0.05,
#'                   main = "Multiple Comparison Procedure: MGM test",
#'                   MCP = c("MGM"))
#'
#' # Other option for the MCP argument is "all". All tests are used.
#'
#' results$Groups     # Results of the tests
#' results$Statistics # Main arguments of the tests
#' results$Summary    # Summary of the response variable
#'
#' # Using the y argument as aov or lm model
#' res  <- aov(rv~treat)
#'
#' MCPtest(y = res, trt = "treat", alpha = 0.05,
#'        main = "Multiple Comparison Procedure: MGM test",
#'        MCP = c("MGM"))
#'
#' # For unbalanced data: It will be used the harmonic mean of
#' #                       the number of experiment replicates
#'
#' # Using the previous example
#' rv <- rv[-1]
#' treat <- treat[-1]
#'
#' res  <- lm(rv~treat) # Linear model
#'
#' # Multiple comparison procedure: MGR test
#' MCPtest(y = res, trt = "treat", alpha = 0.05,
#'        main = "Multiple Comparison Procedure: MGR test",
#'        MCP = c("MGR"))
#'
#' # Assuming that the available data are the averages
#' #  of the treatments and the analysis of variance
#'
#' # Analysis of Variance Table
#'
#' # Response: rv
#' #            Df Sum Sq Mean Sq F value    Pr(>F)
#' # treat      4 4135.2 1033.80  14.669 4.562e-05 ***
#' # Residuals 15 1057.1   70.47
#'
#' mean.treat <- c(100.87, 105.95, 117.62, 127.97, 140.30)
#' treat <- factor(LETTERS[1:5])
#' DFerror <- 15
#' MSerror <- 70.47488
#' replic <- 4
#'
#' MCPtest(y = mean.treat,
#'        trt = treat,
#'        dferror = DFerror,
#'        mserror = MSerror,
#'        replication = replic,
#'        alpha = 0.05,
#'        main = "Multiple Comparison Procedure: MGM test",
#'        MCP = c("MGM"),
#'        ismean = TRUE)
#'
#' @import "utils" "graphics" "SMR" "foreach"
#' @importFrom "stats" "deviance" "df.residual" "qtukey" "aov" "pchisq" "sd"
#' @importFrom "parallel" "makeCluster" "detectCores" "stopCluster"
# @importFrom "doParallel" "registerDoParallel"
#' @export
MCPtest <- function(y, trt = NULL, dferror = NULL, mserror = NULL, replication = NULL, alpha = 0.05, main = NULL,
                   MCP = "all", ismean = FALSE, parallel = FALSE,
                   console = TRUE){
  #####################################################
  #Defensive programming
  if (ismean == TRUE) {
    if (is.null(replication)) {
      stop("The replication argument must be informed", call. = FALSE)
    }
  }
  if (is.numeric(y)) {
    if (!is.numeric(y)) {
      stop("The y argument must be numeric", call. = FALSE)
    }
    if (!is.factor(trt)) {
      stop("The trt argument must be factor", call. = FALSE)
    }
    if (length(y) != length(trt)) {
      stop("The y and trt arguments must have same length", call. = FALSE)
    }
    if (is.null(dferror)) {
      stop("The dferror argument must be informed", call. = FALSE)
    }
    if (is.null(mserror)) {
      stop("The mserror argument must be informed", call. = FALSE)
    }
  }
  #####################################################
  if (all(MCP == "all")){
    MCP = c("MGM", "MGR", "SNKM", "TM", "SK")  # FAZER ALTERACOES
  }
  #####################################################
  #Defensive programming
  if (is.null(trt)) {
    stop("The trt argument is required", call. = FALSE)
  }
  mcps <- c("MGM", "MGR", "SNKM", "TM", "SK") # FAZER ALTERACOES
  nas  <- pmatch(MCP, mcps)
  if (any(is.na(nas))) {
    stop("The options for the MCP argument are 'MGM', 'MGR', 'SNKM', 'TM' and 'SK'", call. = FALSE) # FAZER ALTERACOES
  }
  ################################################
  name.y   <- paste(deparse(substitute(y)))
  name.trt <- paste(deparse(substitute(trt)))
  if (is.null(main)) {
    main  <- paste(name.y, "~", name.trt)
  }
  #y: "aov" or "lm" class
  if ("aov" %in% class(y) | "lm" %in% class(y)) {
    if (is.null(main)) {
      main <- y$call
    }
    A <- y$model
    dferror <- df.residual(y)
    mserror <- deviance(y)/dferror
    y <- A[, 1]
    if (length(trt) > length(names(A)[-1]) + 1) {
      stop("The length of the trt argument is invalid", call. = FALSE)
    }
    ipch <- pmatch(trt, names(A))
    if (any(is.na(ipch))) {
      ncolnam <- length(colnames(A))
      optrt   <- colnames(A)[2]
      if (ncolnam > 2) {
        for (i in 3:ncolnam) {
          optrt <- paste(optrt, names(A)[i], sep = " ")
        }
      }
      stop("The options for the trt argument are\n Options: ", optrt, "\n Note: Observe which of the options have signified practice", call. = FALSE)
    }
    ipch  <- pmatch(trt, names(A))
    nipch <- length(ipch)
    ntrt  <- length(names(A))
    if (any(is.na(ipch))) {
      optrt <- names(A)[2]
      for (i in 3:ntrt) {
        optrt <- paste(optrt, names(A)[i], sep = " ")
      }
      stop("Any of the options of trt argument is wrong \n Options: ", optrt, call. = FALSE)
    }
    name.t <- names(A)[ipch][1]
    trt    <- A[, ipch]
    if (nipch > 1) {
      trt <- A[, ipch[1]]
      for (i in 2:nipch) {
        name.t <- paste(name.t, names(A)[ipch][i], sep = ":")
        trt <- paste(trt, A[, ipch[i]], sep = ":")
      }
    }
    name.y <- names(A)[1]
  }
  #if (console) cat(gettext("Multiple Comparison Procedures\n", domain = "R-MCP"))
  #if (console) cat(gettext("Study: ", domain = "R-MCP"), main, "\n\n")

  # Mean standard error
  ep <- function(x){
    epm <- sd(x)/sqrt(length(x))
    return(epm)
  }

  if (ismean == TRUE) {
    data.na    <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means      <- round(y, 2)
    n          <- length(y)
    rn         <- rep(replication, n)
    rh         <- 1/mean(1/rn)
    Mean       <- mean(y)
    CV         <- sqrt(mserror) * 100/Mean
  } else {
    data.na    <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means      <- round(tapply(data.na[, 1], data.na[, 2],  "mean"), 2)
    std        <- round(tapply(data.na[, 1], data.na[, 2], "ep"), 2)
    rn         <- tapply(data.na[, 1], data.na[, 2], "length")
    mi         <- round(tapply(data.na[, 1], data.na[, 2], "min"), 2)
    ma         <- round(tapply(data.na[, 1], data.na[, 2], "max"), 2)
    n          <- nrow(means)
    rh         <- 1/mean(1/rn)
    Mean       <- mean(data.na[, 1])
    CV         <- sqrt(mserror) * 100/Mean
  }


  ################################################
  #Observation for unbalanced data
  if (length(unique(rn)) != 1) {
    #if (console) cat(gettext("Unbalanced data: It will be used the harmonic mean of \n the number of experiment replicates \n", domain = "R-MCP"))
  }
  ################################################

  if (ismean == TRUE) {
    summarydata <- data.frame(Means  = means,
                              r      = rh)
    #if (console) cat(gettext("Summary:\n", domain = "R-MCP"))
    colnames(summarydata) <- c(gettext("Means", domain = "R-MCP"),
                               "r"
    )
    #if (console) cat(summarydata)
  } else {
    summarydata <- data.frame(Means  = means,
                              std    = std,
                              r      = rh,
                              Min    = mi,
                              Max    = ma)
    #if (console) cat(gettext("Summary:\n", domain = "R-MCP"))
    colnames(summarydata) <- c(gettext("Means", domain = "R-MCP"),
                               gettext("std", domain = "R-MCP"),
                               "r", "Min", "Max"
                               )
    #if (console) print(summarydata)
  }

  if (length(unique(rn)) != 1) {
    #if (console) cat(gettext("\n Harmonic mean of the number of experiment replicates",
    #            domain = "R-MCP"), rh, "\n")
  }
  #DMS midrange
  if (any(MCP == "SNKM")) {
    nnn      <- 2:n
    aaa      <- rep(1 - alpha / 2, times = (n - 1))
    MRq.snkm <- SMR::qSMR(aaa, nnn, dferror) # Studentized midrange of the SNKM test
    MRq      <- MRq.snkm[n - 1] # Studentized midrange of the MGM/TM test
    dms.int  <- MRq * sqrt(mserror / rh) # Internal DMS of the MGM/TM test
    dms      <- dms.int + sqrt(0.5) * sqrt(mserror / rh) / n^0.5 # add of mean square error/nmed^0.5
    dms      <- c(dms, dms.int) # External DMS of the MGM/TM test
    dmsint   <- MRq.snkm * sqrt(mserror / rh) # Internal DMS of the SNKM test
    dms.snkm <- dmsint + sqrt(0.5) * sqrt(mserror / rh) / n^0.5 # External DMS of the SNKM test
  } else{
    MRq     <- SMR::qSMR(1 - alpha / 2, n, dferror) # Studentized midrange of the MGM/TM test
    dms.int <- MRq * sqrt(mserror / rh) # Internal DMS of the MGM/TM test
    dms     <- dms.int + sqrt(0.5) * sqrt(mserror/rh) / n^0.5 # External DMS of the MGM/TM test
    dms     <- c(dms, dms.int) # DMS of the MGM/TM test
  }

  #DMS range:
  Rq   <- qtukey(1 - alpha, n, dferror) # Studentized range of the MGR test
  dms.range <- Rq * sqrt(mserror / rh) # DMS of the MGR test

  #Initial statistics:
  statistics.MGM  <- NA
  statistics.MGR  <- NA
  statistics.SNKM <- NA
  statistics.TM   <- NA
  statistics.SK   <- NA

  #Initial groups:
  test.MGM  <- NA
  test.MGR  <- NA
  test.SNKM <- NA
  test.TM   <- NA
  test.SK   <- NA

  #Defensive programming
  if (!any(parallel == c(TRUE,FALSE))) {
    stop("The options for the parallel argument are TRUE or FALSE!", call. = FALSE)
  }


  # Means Grouping Midrange Test (BATISTA, 2016)
  if (any(MCP == "MGM")){
    #if (console) cat(gettext("\nMeans Grouping Midrange Test (BATISTA, 2016)\n\n", domain = "R-MCP"))
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             DF      = dferror,
                             n       = n,
                             Stud.Midrange = MRq,
                             Ext.DMS = dms[1],
                             Int.DMS = dms[2])
    #if(console) cat(gettext("Statistics: \n", domain = "R-MCP"))
    rownames(statistics) <- " "
    colnames(statistics) <- c(gettext("Exp.Mean", domain = "R-MCP"),
                                "CV",
                                gettext("MSerror", domain = "R-MCP"),
                                gettext("DF", domain = "R-MCP"),
                                "n",
                                gettext("Stud.Midrange", domain = "R-MCP"),
                                gettext("Ext.DMS", domain = "R-MCP"),
                                gettext("Int.DMS", domain = "R-MCP")
                                )
    statistics.MGM <- statistics
    #if (console) print(statistics)
    test <- MGMtest(y, trt, n, dferror, mserror, alpha, dms)
    test[1] <- round(test[1], 2)
    test.MGM <- test
    #if(console) cat(gettext("\nGroups: \n", domain = "R-MCP"))
    #if (console) print(test)
  }

  # Means Grouping based on the Range Test (BATISTA, 2016)
  if (any(MCP == "MGR")){
    #if (console) cat(gettext("\nMean Grouping Range Test (BATISTA, 2016)\n\n", domain = "R-MCP"))
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             DF      = dferror,
                             n       = n,
                             Stud.Range = Rq,
                             DMS = dms.range)
    #if (console) cat(gettext("Statistics: \n", domain = "R-MCP"))
    rownames(statistics) <- " "
    colnames(statistics) <- c(gettext("Exp.Mean", domain = "R-MCP"),
                                "CV",
                                gettext("MSerror", domain = "R-MCP"),
                                gettext("DF", domain = "R-MCP"),
                                "n",
                                gettext("Stud.Range", domain = "R-MCP"),
                                "DMS"
    )
    statistics.MGR   <- statistics
    #if (console) print(statistics)
    test <- MGRtest(y, trt, n, dferror, mserror, alpha, dms.range)
    test[1] <- round(test[1], 2)
    test.MGR <- test
    #if (console) cat(gettext("\nGroups: \n", domain = "R-MCP"))
    #if (console) print(test)
  }



  # Student-Newman-Keuls (SNK) Midrange Test (BATISTA, 2016)
  if (any(MCP == "SNKM")){
    #if (console) cat(gettext("\nSNK Midrange Test (BATISTA, 2016)\n\n", domain = "R-MCP"))
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             DF      = dferror,
                             n       = nnn,
                             Stud.Midrange = MRq.snkm,
                             DMS = dms.snkm)
    statistics <- statistics[order(statistics[, 6], decreasing = FALSE),]
    cont <- seq(1, n - 1, by = 1)
    rownames(statistics) <- cont
    rownames(statistics) <- rownames(MRq.snkm, do.NULL = FALSE, prefix = "comp")
    #if (console) cat(gettext("Statistics: \n", domain = "R-MCP"))
    colnames(statistics) <- c(gettext("Exp.Mean", domain = "R-MCP"),
                                "CV",
                                gettext("MSerror", domain = "R-MCP"),
                                gettext("DF", domain = "R-MCP"),
                                "n",
                                gettext("Stud.Midrange", domain = "R-MCP"),
                                "DMS"
    )
    statistics.SNKM <- round(statistics, 4)
    #if (console) print(round(statistics, 4))
    test <- SNKMtest(y, trt, n, dferror, mserror, alpha, dms.snkm)
    test[1] <- round(test[1], 2)
    test.SNKM <- test
    #if (console) cat(gettext("\nGroups: \n", domain = "R-MCP"))
    #if (console) print(test)
  }

  # Tukey Midrange Test (BATISTA, 2016)
  if (any(MCP == "TM")) {
    #if (console) cat(gettext("\nTukey Midrange Test (BATISTA, 2016)\n\n", domain = "R-MCP"))
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             DF      = dferror,
                             n       = n,
                             Stud.Midrange = MRq,
                             Ext.DMS = dms[1],
                             Int.DMS = dms[2])
    #if (console) cat(gettext("Statistics: \n", domain = "R-MCP"))
    rownames(statistics) <- " "
    colnames(statistics) <- c(gettext("Exp.Mean", domain = "R-MCP"),
                                "CV",
                                gettext("MSerror", domain = "R-MCP"),
                                gettext("DF", domain = "R-MCP"),
                                "n",
                                gettext("Stud.Midrange", domain = "R-MCP"),
                                gettext("Ext.DMS", domain = "R-MCP"),
                                gettext("Int.DMS", domain = "R-MCP")
    )
    statistics.TM <- statistics
    #if (console) print(statistics)

    test <- TMtest(y, trt, n, dferror, mserror, alpha, dms)
    test[1] <- round(test[1], 2)
    test.TM <- test
    #if (console) cat(gettext("\nGroups: \n", domain = "R-MCP"))
    #if (console) print(test)
  }

  # Scott-Knott's test (SCOTT-KNOTT, 1974)
  if (any(MCP == "SK")) {
    #if (console) cat(gettext("\nScott-Knott's Test (SCOTT-KNOTT, 1974)\n\n", domain = "R-MCP"))
    statistics <- data.frame(Exp.Mean = Mean,
                             CV      = CV,
                             MSerror = mserror,
                             DF      = dferror,
                             n       = n
                             )
    #if (console) cat(gettext("Statistics: \n", domain = "R-MCP"))
    rownames(statistics) <- " "
    colnames(statistics) <- c(gettext("Exp.Mean", domain = "R-MCP"),
                                "CV",
                                gettext("MSerror", domain = "R-MCP"),
                                gettext("DF", domain = "R-MCP"),
                                "n"
    )
    statistics.SK <- statistics
    #if (console) print(statistics)

    test <- sktest(y, trt, dferror, mserror, rh, alpha, parallel)
    test.SK <- test
    #if (console) cat(gettext("\nGroups: \n", domain = "R-MCP"))
    #if (console) print(test)
  }

  #All statistics
  stat.tests <- list(Statistics.MGM  = statistics.MGM,
                     Statistics.MGR  = statistics.MGR,
                     Statistics.SNKM = statistics.SNKM,
                     Statistics.TM   = statistics.TM,
                     Statistics.SK   = statistics.SK)

  #All groups
  group.tests <- list(group.MGM  = test.MGM,
                      group.MGR  = test.MGR,
                      group.SNKM = test.SNKM,
                      group.TM   = test.TM,
                      group.SK   = test.SK)
  ################
  # Output results
  ################

  alltest    <- 1:5
  nas        <- nas[order(nas, na.last = NA)]
  ntest      <- alltest[-nas]

  if (length(nas) == 5) {
    statistics <- stat.tests
    grouptest  <- group.tests
  } else {
    statistics <- stat.tests[- ntest]
    grouptest  <- group.tests[- ntest]
  }

  # if(length(nas) == 4) {
  #   grouptest  <- group.tests
  # } else{
  #   grouptest  <- group.tests[- ntest]
  # }

  output <- list(Summary = summarydata,
                    Groups  = grouptest,
                 Statistics = statistics,
                     Tests  = MCP)
  #invisible(output)
  return(output)
}
