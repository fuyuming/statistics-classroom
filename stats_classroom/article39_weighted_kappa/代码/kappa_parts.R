# 输入是4×4联合计数表和同顺序4×4一致权重；输出依次为Po、Pe、κ、一般delta SE、Wald下限、上限。
# 原生点估计/区间由DescTools/vcd另外核对；这里的梯度表示联合格概率改变对κ的一阶影响。
kappa_parts <- function(counts, weight) {
  if (!identical(dim(counts), c(4L, 4L))) stop("counts必须是4×4计数表。")
  if (!identical(dim(weight), c(4L, 4L))) stop("weight必须是4×4一致权重矩阵。")
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) stop("counts必须为有限的非负整数。")
  if (any(!is.finite(weight)) || any(weight < 0) || any(weight > 1)) stop("weight必须为0到1之间的有限数。")
  if (any(abs(diag(weight) - 1) > 1e-12)) stop("weight对角线必须全为1。")
  n <- sum(counts)
  if (!is.finite(n) || n <= 0) stop("配对表总人数必须为正数。")
  # 每格人数/n是联合比例；行列求和分别得到两位评定者的边际比例。
  p <- counts / n
  row <- rowSums(p)
  col <- colSums(p)
  po <- sum(weight * p)
  # outer把行比例与列比例逐一相乘，构造机会参照表。
  pe <- sum(weight * outer(row, col))
  if (1 - pe <= 1e-14) {
    return(c(Po=po, Pe=pe, kappa=NA, SE=NA, CI_low=NA, CI_high=NA))
  }
  kappa <- (po - pe) / (1 - pe)
  # 行、列边际同时随联合概率改变；这个一般delta梯度保留配对关联。
  row_score <- as.vector(weight %*% col) # 对每个行等级汇总机会一致分。
  col_score <- as.vector(row %*% weight) # 对每个列等级汇总机会一致分。
  gradient <- (weight - (1 - kappa) * outer(row_score, col_score, "+")) / (1 - pe)
  variance <- (sum(p * gradient^2) - sum(p * gradient)^2) / n
  se <- sqrt(max(0, variance))
  z <- qnorm(0.975) # 双侧95%区间，在标准正态左侧累计概率0.975处取分位点。
  c(Po=po, Pe=pe, kappa=kappa, SE=se, CI_low=kappa-z*se, CI_high=kappa+z*se)
}
