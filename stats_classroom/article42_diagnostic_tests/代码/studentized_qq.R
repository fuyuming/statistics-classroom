# 输入lm模型；R侧实际调用car::qqPlot，关闭包络以统一三语比较。
studentized_qq <- function(fit, draw = TRUE) {
  z <- sort(as.numeric(rstudent(fit)))
  n <- length(z)
  degrees <- fit$df.residual-1
  if (degrees <= 0 || any(!is.finite(z))) stop("invalid studentized residuals")
  probability <- ppoints(n)
  theoretical <- qt(probability, df = degrees)
  sample_quartiles <- quantile(z, c(0.25,0.75), type = 7, names = FALSE)
  theory_quartiles <- qt(c(0.25,0.75), df = degrees)
  slope <- diff(sample_quartiles)/diff(theory_quartiles)
  intercept <- sample_quartiles[1]-slope*theory_quartiles[1]
  points <- data.frame(probability = probability, theoretical = theoretical,
                       studentized = z, reference_line = intercept+slope*theoretical)
  if (isTRUE(draw)) {
    car::qqPlot(fit, distribution = "t", line = "quartiles", envelope = FALSE,
                simulate = FALSE, id = FALSE, main = "Studentized residual Q-Q plot")
  }
  invisible(list(points = points, df = degrees, slope = slope, intercept = intercept))
}
