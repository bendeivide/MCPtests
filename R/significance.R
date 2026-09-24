signif_code <- function(p, alpha = 0.05) {

  ifelse(
    p <= 0.001, "***",
    ifelse(
      p <= 0.01, "**",
      ifelse(
        p <= 0.05, "*",
        ifelse(
          p <= 0.1, ".",
          ""
        )
      )
    )
  )
}
