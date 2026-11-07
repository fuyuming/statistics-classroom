# 样本率与率差：确定性二项分布与率差卷积（UTF-8，base R）
here <- getwd()
if (basename(here) != "代码") here <- file.path(here, "代码")
root <- normalizePath(file.path(here, ".."))
out <- file.path(root, "运行结果", "R")
fig <- file.path(out, "配图")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
data <- read.csv(file.path(root, "数据", "共享教学数据.csv"), check.names = FALSE)
cell <- function(k, m, s, w = 1) pnorm(k + w / 2, m, s) - pnorm(k - w / 2, m, s)
font <- if ("STHeiti" %in% names(pdfFonts())) "STHeiti" else "sans"
one <- function(n, p, name, title, hi = n) {
  k <- 0:n
  e <- dbinom(k, n, p)
  s <- sqrt(n * p * (1 - p))
  a <- cell(k, n * p, s)
  a[1] <- pnorm(.5, n * p, s)
  a[n + 1] <- pnorm(n - .5, n * p, s, lower.tail = FALSE)
  write.csv(data.frame(k = k, rate = k / n, exact_pmf = e, normal_cell = a), file.path(out, paste0(name, ".csv")), row.names = FALSE)
  sh <- k <= hi
  fn <- if (name == "one_p05_100") "01_样本率精确与近似.png" else if (name == "one_p01_100") "01_样本率精确与近似_低概率.png" else paste0("01_样本率精确与近似_", name, ".png")
  png(file.path(fig, fn), 1450, 850, res = 160)
  par(family = font)
  plot(k[sh], e[sh], type = "h", lwd = 8, col = "#4C78A8", xlab = "计数 k（X）", ylab = "概率", main = title)
  points(k[sh], a[sh], pch = 16, col = "#E45756")
  legend("topright", c("精确概率", "连续性校正正态单格概率"), col = c("#4C78A8", "#E45756"), pch = c(15, 16), bty = "n")
  dev.off()
  cat(name, " rate_variance(pq/n)=", p * (1 - p) / n, " P(X=0)=", (1 - p)^n, "\n")
}
for (i in which(data$kind == "one")) {
  r <- data[i, ]
  n <- as.integer(r$n)
  p <- as.numeric(r$p)
  hi <- if (p < .1) 8 else n
  one(n, p, as.character(r$scenario), paste0("样本率：n=", n, "，p=", p, ""), hi)
}
diff <- function(n, p1, p2) {
  d <- (-n):n
  a <- dbinom(0:n, n, p1)
  b <- dbinom(0:n, n, p2)
  e <- sapply(d, function(v) {
    j <- max(0, -v):min(n, n - v)
    sum(a[v + j + 1] * b[j + 1])
  })
  s <- sqrt(p1 * (1 - p1) / n + p2 * (1 - p2) / n)
  m <- p1 - p2
  ac <- cell(d / n, m, s, 1 / n)
  ac[1] <- pnorm(-1 + .5 / n, m, s)
  ac[length(ac)] <- pnorm(1 - .5 / n, m, s, lower.tail = FALSE)
  data.frame(difference_count = d, D = d / n, exact_pmf = e, normal_cell = ac)
}
dr <- data[data$kind == "difference", ]
zs <- lapply(seq_len(nrow(dr)), function(i) diff(as.integer(dr$n[i]), as.numeric(dr$p1[i]), as.numeric(dr$p2[i])))
for (i in seq_along(zs)) write.csv(zs[[i]], file.path(out, paste0(dr$scenario[i], ".csv")), row.names = FALSE)
png(file.path(fig, "02_率差精确与近似.png"), 1500, 1500, res = 160)
par(mfrow = c(2, 1), family = font, mar = c(4, 4, 3, 1))
for (i in seq_along(zs)) {
  z <- zs[[i]]
  n <- as.integer(dr$n[i])
  sh <- abs(z$difference_count) <= ifelse(n >= 1000, 8, 20)
  plot(z$difference_count[sh], z$exact_pmf[sh], type = "h", lwd = 6, col = "#59A14F", xlab = "差值计数 k=X1−X2（D=k/n）", ylab = "概率", main = paste0("率差：n=", n, "，p1=", dr$p1[i], "，p2=", dr$p2[i]))
  points(z$difference_count[sh], z$normal_cell[sh], pch = 16, col = "#F28E2B")
}
dev.off()
x1 <- 60
n1 <- 300
x2 <- 20
n2 <- 200
r1 <- x1 / n1
r2 <- x2 / n2
pp <- (x1 + x2) / (n1 + n2)
se0 <- sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
zz <- (r1 - r2) / se0
se <- sqrt(r1 * (1 - r1) / n1 + r2 * (1 - r2) / n2)
write.csv(data.frame(pooled_p = pp, z = zz, two_sided_p = 2 * pnorm(abs(zz), lower.tail = FALSE), pearson_Q = zz^2, unpooled_se = se, wald95_low = r1 - r2 - 1.96 * se, wald95_high = r1 - r2 + 1.96 * se), file.path(out, "observed_comparison.csv"), row.names = FALSE)
# 同一数据对照连续性校正开关。FALSE用于与本文未校正Z平方核对。
uncorrected_test <- prop.test(c(x1, x2), c(n1, n2), correct = FALSE)
corrected_test <- prop.test(c(x1, x2), c(n1, n2), correct = TRUE)
print(uncorrected_test)
print(corrected_test)
# 本次两组且未校正时，返回区间采用未合并标准误与精确正态分位数。
expected_interval <- r1 - r2 + c(-1, 1) * qnorm(0.975) * se
stopifnot(isTRUE(all.equal(as.numeric(uncorrected_test$conf.int), expected_interval)))
comparison <- data.frame(
  correction = c(FALSE, TRUE),
  pearson_Q = c(unname(uncorrected_test$statistic), unname(corrected_test$statistic)),
  df = c(unname(uncorrected_test$parameter), unname(corrected_test$parameter)),
  two_sided_p = c(uncorrected_test$p.value, corrected_test$p.value)
)
write.csv(comparison, file.path(out, "continuity_correction_comparison.csv"), row.names = FALSE)
print(binom.test(0, 20))
write.csv(data.frame(x = 0, n = 20, wald_low = 0, wald_high = 0, clopper_pearson_low = 0, clopper_pearson_high = qbeta(.975, 1, 20)), file.path(out, "boundary_zero_20.csv"), row.names = FALSE)

# Cohen's h：用样本率计算反正弦变换后的效应量（方向为第一组减第二组）。
# asin() 接受0到1的数值，返回弧度角；h的大小分类只作经验参考，医学意义仍需结合结局判断。
cohens_h <- function(p1, p2) {
  2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
}
conventional_band <- function(h) {
  ah <- abs(h)
  if (ah < .2) {
    "below_small_reference"
  } else if (ah < .5) {
    "small"
  } else if (ah < .8) {
    "medium"
  } else {
    "large"
  }
}
effect <- data.frame(
  p1 = c(r1, .55),
  p2 = c(r2, .45)
)
effect$risk_difference <- effect$p1 - effect$p2
effect$difference_percentage_points <- 100 * effect$risk_difference
effect$cohens_h <- cohens_h(effect$p1, effect$p2)
effect$abs_h <- abs(effect$cohens_h)
effect$conventional_band <- vapply(effect$cohens_h, conventional_band, character(1))
write.csv(effect, file.path(out, "effect_size_comparison.csv"), row.names = FALSE)
cat("\nCohen's h（反正弦变换后的效应量；方向为第一组减第二组）\n")
cat("h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2))；医学意义需结合结局、临床阈值和研究设计判断。\n")
cat(sprintf("主例：第一组率=%.6f，第二组率=%.6f，率差=%.1f百分点，h=%.13f，惯例分类=%s\n", r1, r2, 100 * (r1-r2), cohens_h(r1, r2), conventional_band(cohens_h(r1, r2))))
cat("惯例参考：|h|<.2低于小效应参考值；.2<=|h|<.5小；.5<=|h|<.8中；>=.8大；不代表医学重要性。\n")
print(effect)


# 图3的单根柱形怎样算 -----------------------------------------
# 教学模型给定两组总体率均为0.30，率差0.02对应相差1名阳性。
n <- 50
p1 <- 0.30
p2 <- 0.30
x2 <- 0:(n - 1)
x1 <- x2 + 1
# 每对计数同时出现的概率，由两组独立的二项概率相乘得到。
prob_each_pair <- dbinom(x1, size = n, prob = p1) *
  dbinom(x2, size = n, prob = p2)
# 50种互斥组合都产生同一率差，概率相加得到约0.08488994。
exact_probability <- sum(prob_each_pair)
cat("图3中D=0.02的精确概率\n")
print(exact_probability)
write.csv(
  data.frame(x1 = x1, x2 = x2, rate1 = x1 / n, rate2 = x2 / n,
             difference = (x1 - x2) / n, probability = prob_each_pair),
  file.path(out, "图3单柱概率核对.csv"), row.names = FALSE
)
