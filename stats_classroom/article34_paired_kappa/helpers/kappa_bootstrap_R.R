# 第34篇 Bootstrap helper：四格计数的整数式 Cohen κ及R boot端点。
kappa_from_counts <- function(v) {
  stopifnot(length(v) == 4L, all(is.finite(v)), all(v >= 0), sum(v) > 0)
  a <- v[1]; b <- v[2]; cc <- v[3]; d <- v[4]
  den <- (a+b)*(b+d) + (cc+d)*(a+cc)
  if (den == 0) return(NA_real_)
  2 * (a*d - b*cc) / den
}
bootstrap_summary_row <- function(stat, kappa_hat, jack, B, seed,
                                   percentile, bca, alpha_bca, shared_counts=NULL) {
  invalid <- sum(!is.finite(stat[, 1]))
  if (invalid != 0L) stop("发现不可定义的bootstrap样本；本例不允许静默删除。")
  if (!is.null(shared_counts) && !all(unname(stat[,2:5]) == unname(shared_counts)))
    stop("boot::boot生成的计数与共享bootstrap_counts.csv不一致。")
  # jack有n个delete-one值；影响值的倍数是n−1，而不是Bootstrap重复次数B。
  L <- (length(jack)-1L) * (mean(jack) - jack)
  z0 <- qnorm(mean(stat[,1] < kappa_hat))
  acc <- sum(L^3) / (6 * sum(L^2)^1.5)
  data.frame(B=B, seed=seed, invalid=invalid, kappa=kappa_hat,
    bootstrap_mean=mean(stat[,1]), bootstrap_bias=mean(stat[,1])-kappa_hat,
    bootstrap_variance=var(stat[,1]), bootstrap_se=sd(stat[,1]),
    bias_probability=mean(stat[,1] < kappa_hat), z0=z0, acceleration=acc,
    alpha_low=alpha_bca[1], alpha_high=alpha_bca[2],
    percentile_low=percentile[1], percentile_high=percentile[2],
    bca_low=bca[1], bca_high=bca[2])
}
