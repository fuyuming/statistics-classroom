# 统计课堂24：8名受试者×4种声音条件
# 打开根目录的统计课堂24.Rproj，再打开本文件，点击Source。
# 也可在项目根目录执行 Rscript 代码/repeated_anova.R。
# ① 读取原表：每行一人，每列一个条件；数值单位为%。 ----
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
candidates <- c(getwd(), dirname(getwd()))
if (length(script)) candidates <- c(dirname(dirname(normalizePath(script))), candidates)
root <- candidates[vapply(candidates, function(p) file.exists(file.path(p, "数据", "四种声音条件反应率.csv")), logical(1))][1]
if (is.na(root)) stop("请先打开项目根目录的统计课堂24.Rproj。")
for (pkg in c("rstatix", "effectsize")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("请先安装扩展包：", pkg))
}
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(file.path(root, "数据", "四种声音条件反应率.csv"))
conds <- c("A", "B", "C", "D")
x0 <- as.matrix(dat[, conds])

# 下面封装一次完整分析，最后对原表和练习表分别运行。
run_one <- function(x, label) {
  n <- nrow(x)
  k <- ncol(x)
  colnames(x) <- conds
  save_csv <- function(z, name) write.csv(z, file.path(out, paste0(label, "_", name, ".csv")), row.names = FALSE)

  # ② 长表与标准函数：同一个人的4行使用同一个subject编号。 ----
  long <- data.frame(
    subject = factor(rep(dat$subject, each = k)),
    condition = factor(rep(conds, n), levels = conds),
    response = as.vector(t(x))
  )
  save_csv(long, "long_data")
  fit <- aov(response ~ condition + Error(subject/condition), data = long)
  capture.output(summary(fit), file = file.path(out, paste0(label, "aov.txt")))
  fit_table <- summary(fit)[["Error: subject:condition"]][[1]]
  rf <- rstatix::anova_test(data = long, dv = response, wid = subject,
                          within = condition, effect.size = "pes")
  save_csv(as.data.frame(rf$`Mauchly's Test for Sphericity`), "mauchly")
  save_csv(as.data.frame(rstatix::get_anova_table(rf, correction = "GG")), "rstatix_GG")

  # ③ 手动计算：全体均值、各条件均值、各人的平均值。 ----
  grand <- mean(x)
  cm <- colMeans(x)
  sm <- rowMeans(x)
  ss_total <- sum((x - grand)^2)
  ss_condition <- n * sum((cm - grand)^2)  # 每列有n个人
  ss_subject <- k * sum((sm - grand)^2)    # 每人有k次测量
  ss_error <- ss_total - ss_condition - ss_subject
  # sweep的第二个参数：1按行，2按列。避免R向量循环造成错误减法。
  centered <- sweep(x, 1, sm, "-")
  remain <- sweep(centered, 2, cm, "-") + grand
  stopifnot(abs(sum(remain^2) - ss_error) < 1e-9)
  save_csv(as.data.frame(centered), "person_centered")
  save_csv(as.data.frame(remain), "centered_residual")
  save_csv(data.frame(condition = conds, condition_mean = cm), "condition_means")
  save_csv(data.frame(subject = dat$subject, subject_mean = sm), "subject_means")

  df1 <- k - 1
  df_subject <- n - 1
  df2 <- df1 * df_subject
  ms_c <- ss_condition / df1
  ms_e <- ss_error / df2
  fval <- ms_c / ms_e
  pval <- pf(fval, df1, df2, lower.tail = FALSE)
  # 与标准aov的实际输出比较，而非仅运行函数。
  stopifnot(abs(fval - fit_table$`F value`[1]) < 1e-9,
            abs(pval - fit_table$`Pr(>F)`[1]) < 1e-12)

  # ④ GG：两个自由度都乘epsilon；推导将在第25篇展开。 ----
  # D在每个人内部做中心化；S是四列的样本协方差矩阵。
  D <- diag(k) - matrix(1/k, k, k)
  S <- cov(x)
  DS <- D %*% S
  eps <- sum(diag(DS))^2 / ((k - 1) * sum(diag(DS %*% DS)))
  gg1 <- eps * df1
  gg2 <- eps * df2
  pgg <- pf(fval, gg1, gg2, lower.tail = FALSE)

  # ⑤ 偏eta平方及经典双侧95%区间：明确使用未GG调整的F模型。 ----
  eta <- ss_condition / (ss_condition + ss_error)
  # 对非中心参数求根，使观察F左侧概率分别等于.975和.025。
  ncp <- sapply(c(.975, .025), function(target) {
    fn <- function(z) suppressWarnings(pf(fval, df1, df2, ncp = z) - target)
    if (fn(0) <= 0) return(0)
    hi <- 1
    while (fn(hi) > 0) hi <- hi * 2
    uniroot(fn, c(0, hi), tol = 1e-10)$root
  })
  etaci <- ncp / (ncp + df2)
  eta_function <- as.data.frame(effectsize::eta_squared(
    fit, partial = TRUE, ci = .95, alternative = "two.sided"
  ))
  save_csv(eta_function, "effectsize")
  stopifnot(abs(eta_function$Eta2_partial[1] - eta) < 1e-9,
            max(abs(c(eta_function$CI_low[1], eta_function$CI_high[1]) - etaci)) < 1e-4)

  rows <- list(
    data.frame(scenario=label, measure="repeated_rm", term="condition", SS=ss_condition, df=df1, MS=ms_c, F=fval, P=pval, partial_eta2=eta, CI_low=etaci[1], CI_high=etaci[2]),
    data.frame(scenario=label, measure="repeated_rm", term="subject", SS=ss_subject, df=df_subject, MS=ss_subject/df_subject),
    data.frame(scenario=label, measure="repeated_rm", term="error", SS=ss_error, df=df2, MS=ms_e),
    data.frame(scenario=label, measure="repeated_rm", term="total", SS=ss_total, df=n*k-1),
    data.frame(scenario=label, measure="GG", term="condition", epsilon_GG=eps, df1_GG=gg1, df2_GG=gg2, F=fval, P_GG=pgg)
  )
  # 错误示范：忽略subject，假装四列来自独立组；只能用于比较算法。
  wrong <- summary(aov(response ~ condition, data=long))[[1]]
  rows[[6]] <- data.frame(scenario=label, measure="independent_wrong_demo", term="condition",
                         F=wrong$`F value`[1], P=wrong$`Pr(>F)`[1], df1=k-1, df2=n*k-k)

  # ⑥ 六个配对t检验：每对用自己的差值SD；方向为后列−前列。 ----
  pairs <- combn(k, 2, simplify = FALSE)
  m <- length(pairs)
  paired <- lapply(pairs, function(ij) {
    i <- ij[1]; j <- ij[2]
    d <- x[, j] - x[, i]
    tt <- t.test(x[, j], x[, i], paired = TRUE, conf.level = 1-.05/m)
    se <- sd(d) / sqrt(n)
    half <- qt(1-.05/(2*m), n-1) * se
    stopifnot(max(abs(tt$conf.int - c(mean(d)-half, mean(d)+half))) < 1e-9)
    data.frame(scenario=label, measure="paired_t", term=paste0(conds[j], "-", conds[i]),
               mean_diff_percentage_points=mean(d), SD=sd(d), SE=se, df=n-1,
               t=unname(tt$statistic), rawP=tt$p.value, bonferroniP=min(m*tt$p.value,1),
               CI_level=1-.05/m, CI_low=tt$conf.int[1], CI_high=tt$conf.int[2],
               Shapiro_P=shapiro.test(d)$p.value)
  })
  paired <- do.call(rbind, paired)
  stopifnot(max(abs(p.adjust(paired$rawP, "bonferroni") - paired$bonferroniP)) < 1e-12)
  save_csv(paired, "paired_t")
  # 六张差值QQ图；每个面板与正文的对应差值相同。
  pdf(file.path(out, paste0(label, "_difference_QQ.pdf")), width=9, height=6)
  par(mfrow=c(2,3))
  for (ij in pairs) {
    d <- x[, ij[2]] - x[, ij[1]]
    qqnorm(d, main=paste0(conds[ij[2]], "-", conds[ij[1]]))
    qqline(d, col="#C7664D")
  }
  dev.off()
  rows <- c(rows, split(paired, seq_len(nrow(paired))))
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(z) {
    z[setdiff(all_names, names(z))] <- NA
    z[all_names]
  })
  save_csv(do.call(rbind, rows), "summary")
  c(F=fval, P=pval, epsilon=eps, P_GG=pgg, eta=eta)
}

original <- run_one(x0, "original")
# ⑦ 练习：给同一个人的四个数加同一个数，保留所有条件差值。 ----
offsets <- seq(0, 14, by=2)
changed <- run_one(sweep(x0, 1, offsets, "+"), "offsets_0_14")
stopifnot(max(abs(original - changed)) < 1e-10)
print(original)
cat("原表与练习表计算完毕，标准函数核对通过。\n")
