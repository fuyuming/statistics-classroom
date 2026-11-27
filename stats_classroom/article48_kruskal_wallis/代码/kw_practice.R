# 统计课堂48：H构造、六瓶精确枚举、兔例Monte Carlo。UTF-8。
# 24只独立兔，四组各6；每兔两眼形成一个改善差（分）。
# 来源Rosner8e表12.17；六瓶为教学假设资料。完整说明见README。
# RStudio打开根目录Rproj，按节运行；也可source("代码/kw_practice.R")。
# 统计计算使用base R，第5节绘图需要ggplot2。图默认显示，保存默认关。

# 1 准备：项目根目录、参数与共享数据 ----
root <- getwd()
if (!dir.exists(file.path(root, "数据"))) {
  stop("请先打开项目根目录的Rproj，或将工作目录设为项目根。")
}
save_figures <- Sys.getenv("SAVE_FIGURES", "0") == "1"
show_figures <- Sys.getenv("KW_HEADLESS", "0") != "1"
alpha <- 0.05
mc_draws <- 199999L
seed <- 20261130L
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rabbit <- read.csv(file.path(root, "数据", "rabbit.csv"))
cat("每行一只兔，drug为药物，score为改善差（分）：\n")
print(head(rabbit))
stopifnot(!anyNA(rabbit))

# 2 合并秩、B、秩方差、C及H ----
# rank的average参数对并列值取平均秩。全部24只共同编秩。
rabbit$rank <- rank(rabbit$score, ties.method = "average")
order <- unique(rabbit$drug)
N <- nrow(rabbit)
k <- length(order)
center <- (N + 1) / 2
sizes <- numeric(k)
rank_sums <- numeric(k)
mean_ranks <- numeric(k)
for (i in seq_len(k)) {
  selected <- rabbit$drug == order[i]
  sizes[i] <- sum(selected)
  rank_sums[i] <- sum(rabbit$rank[selected])
  mean_ranks[i] <- mean(rabbit$rank[selected])
}
rank_table <- data.frame(drug = order, n = sizes,
                         rank_sum = rank_sums, mean_rank = mean_ranks)
print(rank_table)
B <- sum(sizes * (mean_ranks - center)^2)
rank_total_ss <- sum((rabbit$rank - center)^2)
rank_variance <- rank_total_ss / (N - 1)
counts <- as.numeric(table(rabbit$score))
C <- 1 - sum(counts^3 - counts) / (N^3 - N)
if (rank_variance <= 0) stop("全部评分相同，秩变异为0。")
h_raw <- 12 * B / (N * (N + 1))
H <- B / rank_variance
# pchisq的lower.tail=FALSE返回χ²右尾；自由度为组数减1。
p_chisq <- pchisq(H, df = k - 1, lower.tail = FALSE)
native <- kruskal.test(score ~ drug, data = rabbit)
cat(sprintf("B=%.9f SS_T=%.9f s_r2=%.9f\n", B, rank_total_ss, rank_variance))
cat(sprintf("C=%.12f H_raw=%.9f H=%.9f chi2_P=%.12f\n", C, h_raw, H, p_chisq))
Vbase <- N * (N + 1) / 12
Tbase <- (N - 1) * Vbase
tie_values <- as.numeric(table(rabbit$score))
tie_losses <- (tie_values^3 - tie_values) / 12
Tactual <- Tbase - sum(tie_losses)
Vactual <- rank_variance
cat(sprintf("Vbase=%.12f Tbase=%.12f\n", Vbase, Tbase))
cat("tie frequencies=", paste(tie_values, collapse = ","),
    "losses=", paste(sprintf("%.12f", tie_losses), collapse = ","), "\n")
cat(sprintf("Tactual=%.12f Vactual=%.12f C=Vactual/Vbase=%.12f\n",
            Tactual, Vactual, Vactual / Vbase))
group_var <- Vactual * (1 / sizes - 1 / N)
group_cov <- -Vactual / N
cat("group mean-rank variances=", paste(sprintf("%.12f", group_var), collapse = ","), "\n")
cat(sprintf("covariance between any two distinct groups=%.12f\n", group_cov))
if (length(unique(sizes)) == 1 && k == 4) {
  n <- sizes[1]
  # 逐行命名四个平均秩，便于新手把每个对比代回公式。
  m1 <- mean_ranks[1]
  m2 <- mean_ranks[2]
  m3 <- mean_ranks[3]
  m4 <- mean_ranks[4]
  differences <- c(m1-m2, (m1+m2)/2-m3, (m1+m2+m3)/3-m4)
  SE <- sqrt(c(2*Vactual/n, 3*Vactual/(2*n), 4*Vactual/(3*n)))
  Z <- differences / SE
  Z2 <- Z^2
  contrast_components <- data.frame(
    contrast = c("m1-m2", "(m1+m2)/2-m3", "(m1+m2+m3)/3-m4"),
    difference = differences, SE = SE, Z = Z,
    Z_squared_contribution = Z2)
  write.csv(contrast_components, file.path(out, "contrast_components.csv"), row.names = FALSE)
  print(contrast_components)
  cat(sprintf("contrast Z^2 sum=%.12f; H=%.12f\n", sum(Z2), H))
  contrast_sum_z2 <- sum(Z2)
} else {
  contrast_sum_z2 <- NA_real_
  cat("contrast components skipped: requires four equal-sized groups.\n")
}
print(native) # 返回校正后统计量、近似自由度与P值。
cat("χ²近似判断：", if (p_chisq <= alpha) "拒绝同分布H0。\n" else "不足以拒绝同分布H0。\n")
write.csv(rank_table, file.path(out, "rabbit_ranks.csv"), row.names = FALSE)

# 3 六瓶：观察H来自当前组别，完整枚举全部身份划分 ----
bottles <- read.csv(file.path(root, "数据", "six_bottles.csv"))
print(bottles)
small_ranks <- rank(bottles$concentration_mg_L, ties.method = "average")
small_center <- (nrow(bottles) + 1) / 2
small_var <- var(small_ranks) # var使用N-1分母。
small_sizes <- table(bottles$condition)
stopifnot(nrow(bottles) == 6, length(small_sizes) == 3,
          all(small_sizes == 2), small_var > 0)
small_means <- tapply(small_ranks, bottles$condition, mean)
small_B <- sum(as.numeric(small_sizes) * (small_means - small_center)^2)
small_H <- small_B / small_var
all_H <- numeric(0)
# combn列出选2瓶的全部组合；剩余两瓶自动归C，组名不去重。
for (a in combn(1:6, 2, simplify = FALSE)) {
  remaining <- setdiff(1:6, a)
  for (b in combn(remaining, 2, simplify = FALSE)) {
    cc <- setdiff(remaining, b)
    allocation_means <- c(mean(small_ranks[a]), mean(small_ranks[b]),
                           mean(small_ranks[cc]))
    allocation_B <- sum(2 * (allocation_means - small_center)^2)
    all_H <- c(all_H, allocation_B / small_var)
  }
}
tail_count <- sum(all_H >= small_H - 1e-10) # 包含相等，容纳浮点误差。
exact_P <- tail_count / length(all_H)
small_chisq_P <- pchisq(small_H, df = 2, lower.tail = FALSE)
values <- sort(unique(round(all_H, 12)))
frequencies <- sapply(values, function(x) sum(abs(all_H - x) < 1e-10))
distribution <- data.frame(H = values, frequency = frequencies,
                            probability = frequencies / length(all_H))
print(distribution)
cat(sprintf("辅助 H=%.9f total=%d tail=%d exact_P=%.12f chi2_P=%.12f\n",
            small_H, length(all_H), tail_count, exact_P, small_chisq_P))
cat("精确枚举判断：", if (exact_P <= alpha) "拒绝同分布H0。\n" else "不足以拒绝同分布H0。\n")
write.csv(distribution, file.path(out, "aux_exact_distribution.csv"), row.names = FALSE)

# 4 兔例：固定组大小，Monte Carlo抽取合法划分 ----
set.seed(seed)
ends <- cumsum(sizes)
starts <- c(1, head(ends, -1) + 1)
exceedances <- 0L
for (draw in seq_len(mc_draws)) {
  permuted <- sample(rabbit$rank, replace = FALSE) # 每轮均匀不放回重排。
  perm_B <- 0
  for (i in seq_len(k)) {
    perm_mean <- mean(permuted[starts[i]:ends[i]])
    perm_B <- perm_B + sizes[i] * (perm_mean - center)^2
  }
  perm_H <- perm_B / rank_variance
  exceedances <- exceedances + (perm_H >= H - 1e-10)
}
mc_P <- (exceedances + 1) / (mc_draws + 1)
mc_se <- sqrt(mc_P * (1 - mc_P) / mc_draws)
cat(sprintf("MC B=%d seed=%d exceedances=%d P_plus1=%.9f MCSE=%.9f\n",
            mc_draws, seed, exceedances, mc_P, mc_se))
cat("条件置换判断：", if (mc_P <= alpha) "拒绝同分布H0。\n" else "不足以拒绝同分布H0。\n")
cat("这是整体分布比较，具体药物组对需要另作多重比较。\n")
write.csv(data.frame(H, C, chi2_P = p_chisq, mc_draws,
                     tail_count = exceedances, mc_P, mc_se,
                     small_H, small_exact_P = exact_P,
                     Vbase, Tbase, Tactual, Vactual,
                     contrast_Z2_sum = contrast_sum_z2),
          file.path(out, "summary.csv"), row.names = FALSE)

# 5 绘图：R + ggplot2，也是正文主图的来源 ----
# ggplot建立图对象；geom_step(direction="hv")画右连续的精确CDF阶梯。
# 本节使用第2、3节结果，不读取作者的旧输出。配色与坐标在此统一。
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop('绘图需要ggplot2：请先运行install.packages("ggplot2")。')
}
library(ggplot2)
rank_plot_data <- data.frame(
  drug = factor(order, levels = rev(order)),
  mean_rank = mean_ranks
)
rank_plot <- ggplot(rank_plot_data, aes(x = mean_rank, y = drug)) +
  geom_col(fill = "#2463B5", width = 0.62) +
  geom_vline(xintercept = center, linetype = "dashed", color = "#D66B24") +
  geom_text(aes(label = sprintf("%.2f", mean_rank)), hjust = -0.15, size = 4) +
  scale_x_continuous(limits = c(0, max(mean_ranks) * 1.15), expand = c(0, 0)) +
  labs(x = "Mean rank", y = NULL,
       subtitle = paste("Dashed line: pooled mean =", center)) +
  theme_classic(base_size = 14) +
  theme(plot.subtitle = element_text(size = 11, margin = margin(b = 12)))

exact_curve <- data.frame(
  H = c(-0.1, values, 5.5),
  cumulative = c(0, cumsum(frequencies) / length(all_H), 1)
)
curve_grid <- seq(0, 5.5, length.out = 500)
approx_curve <- data.frame(H = curve_grid, cumulative = pchisq(curve_grid, df = 2))
cdf_plot <- ggplot() +
  geom_step(data = exact_curve, aes(H, cumulative, color = "Exact"),
            direction = "hv", linewidth = 0.9) +
  geom_line(data = approx_curve, aes(H, cumulative, color = "Chi-square(2)"),
            linewidth = 0.9) +
  geom_vline(xintercept = small_H, linetype = "dashed", color = "#777777", linewidth = 0.6) +
  annotate("label", x = ifelse(small_H < 3, 0.2, 5.3),
           y = ifelse(small_H < 3, 0.95, 0.16),
           label = sprintf("Observed H = %.3f", small_H), angle = 0,
           hjust = ifelse(small_H < 3, 0, 1), size = 4.1,
           fill = "white", color = "#555555") +
  scale_color_manual(values = c("Exact" = "#2463B5", "Chi-square(2)" = "#D66B24"),
                     breaks = c("Exact", "Chi-square(2)")) +
  scale_x_continuous(limits = c(-0.1, 5.5), breaks = 0:5, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = "H", y = "Cumulative probability", color = NULL) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top", legend.text = element_text(size = 14),
        legend.key.width = grid::unit(24, "pt"), legend.spacing.x = grid::unit(8, "pt"),
        plot.margin = margin(6, 8, 6, 8))
if (show_figures) {
  print(rank_plot)
  print(cdf_plot)
}
if (save_figures) {
  ggsave(file.path(root, "文章配图", "R_mean_ranks.png"), rank_plot,
         width = 6.5, height = 3.6, dpi = 300, bg = "white")
  ggsave(file.path(root, "文章配图", "R_aux_cdf.png"), cdf_plot,
         width = 4.2, height = 4.2, dpi = 300, bg = "white")
}

# 6 改分组练习 ----
# 保留6个浓度与每组2瓶，将condition按浓度1到6改成A,B,C,C,B,A。
# 重跑第3节：三组平均秩都为3.5，预期H=0、精确P=1；再运行第5节。
