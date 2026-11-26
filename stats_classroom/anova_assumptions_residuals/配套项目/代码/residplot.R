# 按R in Action残差分布图思路改写；三语统一分箱和Gaussian核带宽。
residplot <- function(fit, nbreaks = 10L, draw = TRUE) {
  if (!inherits(fit, "lm") || !is.null(fit$weights)) {
    stop("fit must be an unweighted lm")
  }
  if (fit$rank != ncol(model.matrix(fit))) {
    stop("the model matrix must have full column rank")
  }
  if (length(nbreaks) != 1L || !is.finite(nbreaks) ||
      nbreaks < 2 || nbreaks != as.integer(nbreaks)) {
    stop("nbreaks must be an integer >= 2")
  }
  # rstudent采用删除当前观测后的误差方差，区别于rstandard。
  z <- as.numeric(rstudent(fit))
  n <- length(z)
  if (fit$df.residual <= 1 || any(!is.finite(z)) || diff(range(z)) == 0) {
    stop("finite, varying studentized residuals and df_resid > 1 are required")
  }
  sample_sd <- sd(z)
  quartiles <- quantile(z, c(0.25, 0.75), type = 7, names = FALSE)
  spread <- min(sample_sd, diff(quartiles)/1.34)
  if (spread <= 0) spread <- sample_sd
  bandwidth <- 0.9*spread*n^(-1/5)
  edges <- seq(min(z), max(z), length.out = nbreaks+1L)
  histogram <- hist(z, breaks = edges, plot = FALSE, right = FALSE, include.lowest = TRUE)
  grid <- seq(min(z)-3*bandwidth, max(z)+3*bandwidth, length.out = 512)
  normal <- dnorm(grid, mean(z), sample_sd)
  # 显式计算Gaussian核，和Python、MATLAB使用同一个定义。
  kde <- vapply(grid, function(point) {
    mean(exp(-0.5*((point-z)/bandwidth)^2))/(bandwidth*sqrt(2*pi))
  }, numeric(1))
  out <- list(
    studentized_residuals = z,
    histogram = data.frame(left = head(edges, -1), right = tail(edges, -1),
                           count = histogram$counts, density = histogram$density),
    curves = data.frame(x = grid, normal = normal, kde = kde),
    bandwidth = unname(bandwidth)
  )
  if (isTRUE(draw)) {
    ymax <- max(c(histogram$density, normal, kde))*1.08
    plot(histogram, freq = FALSE, col = "grey78", border = "white",
         xlim = range(grid), ylim = c(0, ymax), main = "Residual distribution",
         xlab = "Externally studentized residual", ylab = "Density")
    lines(grid, normal, col = "blue", lwd = 2)
    lines(grid, kde, col = "red", lty = 2, lwd = 2)
    rug(z, col = "brown")
    legend("topright", c("Fitted normal", "Gaussian KDE", "Residual positions"),
           col = c("blue", "red", "brown"), lty = c(1, 2, 1), bty = "n")
  }
  invisible(out)
}
