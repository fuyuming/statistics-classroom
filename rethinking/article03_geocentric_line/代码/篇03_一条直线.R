# ============================================================
# 精读 03｜12 个时间点，就敢画一条生长曲线吗
# 对应：McElreath 2023 第 03 讲 Geocentric Models；教材第 4 章
# 用法： Rscript 代码/篇03_一条直线.R
#
# 约定：
#   · 主线案例＝教学用模拟数据（生物膜厚度随培养时间，12 个时间点），残差在脚本里固定给出
#   · 对照案例＝教材 Howell1 成年人（体重~身高）；R 先导出 CSV，Python／MATLAB 读同一份
#   · 全部确定性计算：网格后验、解析最小二乘对照、Halton 序列代替随机抽样
#     → 三语言输出逐位一致（含 CSV 文件本身）
#   · 主线后验用「矩形范围内的平坦先验」；图 1 的细线也从这个先验抽（前后一致）
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(rethinking)
})

# 目录自动定位：脚本放在 <篇目录>/代码/ 下，所以篇目录是它的上一级（读者换机器也不用改路径）
.args <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .args, value = TRUE)
root <- if (length(.f) > 0) dirname(dirname(normalizePath(sub("^--file=", "", .f[1])))) else getwd()
OUT  <- file.path(root, "文章配图")
RES  <- file.path(root, "运行结果")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

A_LIM <- c(10, 26); B_LIM <- c(2.0, 6.0); N_GRID <- 200      # 参数框与网格密度（图 1 的细线也在这个框里抽）

save_fig <- function(name, plot, w = 7.4, h = 4.2) {
  ggsave(file.path(OUT, name), plot, width = w, height = h, dpi = 150,
         device = ragg::agg_png, bg = "white")
}
num <- function(x, d = 6) {
  s <- formatC(round(as.numeric(x), d), format = "f", digits = d)
  s <- sub("0+$", "", s); s <- sub("\\.$", "", s)
  ifelse(s == "" | s == "-0", "0", s)
}
write_csv_fixed <- function(df, path) {
  hdr <- paste0('"', names(df), '"', collapse = ",")
  cols <- lapply(df, function(cl) if (is.numeric(cl)) num(cl) else as.character(cl))
  writeLines(c(hdr, do.call(paste, c(cols, sep = ","))), path, useBytes = TRUE)
}
halton <- function(n, base) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    f <- 1; r <- 0; k <- i
    while (k > 0) { f <- f / base; r <- r + f * (k %% base); k <- k %/% base }
    out[i] <- r
  }
  out
}

# ---------- 0) 教材对照数据 ----------
data(Howell1)
d <- Howell1[Howell1$age >= 18, c("height", "weight")]
d <- d[order(d$height), ]
d$height_c <- d$height - mean(d$height)
write_csv_fixed(data.frame(height = d$height, height_c = d$height_c, weight = d$weight),
                file.path(RES, "00_教材数据_Howell1成年人.csv"))

# ---------- 1) 主线案例（教学模拟数据）----------
t_hours   <- 1:12
resid_fix <- c(0.6, -0.9, 1.2, -0.4, 0.8, -1.3, 0.3, 1.1, -0.7, 0.2, -0.5, 1.4)
thickness <- 18 + 3.4 * t_hours + resid_fix
bio <- tibble(时间 = t_hours, 厚度 = thickness)
write_csv_fixed(bio, file.path(RES, "01_我们的案例_生物膜厚度.csv"))

# ---------- 2) 二维网格后验 ----------
grid_post <- function(x, y, n_grid = N_GRID, sigma, a_lim = A_LIM, b_lim = B_LIM,
                      prior_a = NULL, prior_b = NULL) {
  step_a <- (a_lim[2] - a_lim[1]) / (n_grid - 1)
  step_b <- (b_lim[2] - b_lim[1]) / (n_grid - 1)
  a_seq <- a_lim[1] + (0:(n_grid - 1)) * step_a
  b_seq <- b_lim[1] + (0:(n_grid - 1)) * step_b
  g <- expand.grid(a = a_seq, b = b_seq)          # expand.grid：a 变得最快
  M <- outer(x, b_seq, "*")
  b_idx <- (seq_len(nrow(g)) - 1) %/% n_grid + 1
  lp <- numeric(nrow(g))
  for (i in seq_len(nrow(g))) {
    e <- (y - (g$a[i] + M[, b_idx[i]])) / sigma
    lp[i] <- -0.5 * sum(e^2) - length(y) * (log(sigma) + 0.5 * log(2 * pi))
    if (!is.null(prior_a)) lp[i] <- lp[i] - 0.5 * ((g$a[i] - prior_a[1]) / prior_a[2])^2 -
      log(prior_a[2]) - 0.5 * log(2 * pi)
    if (!is.null(prior_b)) lp[i] <- lp[i] - 0.5 * ((g$b[i] - prior_b[1]) / prior_b[2])^2 -
      log(prior_b[2]) - 0.5 * log(2 * pi)
  }
  g$post <- exp(lp - max(lp)); g$post <- g$post / sum(g$post)
  g
}
ols <- function(x, y) {
  b <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
  a <- mean(y) - b * mean(x)
  res <- y - (a + b * x)
  s2 <- sum(res^2) / (length(y) - 2)
  list(a = a, b = b, sigma_hat = sqrt(s2), se_b = sqrt(s2 / sum((x - mean(x))^2)))
}
marg <- function(g, v) {
  s <- tapply(g$post, g[[v]], sum)
  tibble(v = as.numeric(names(s)), p = as.numeric(s) / sum(as.numeric(s)))
}
interval89 <- function(m) {
  c_ <- cumsum(m$p)
  c(m$v[which(c_ >= 0.055)[1]], m$v[which(c_ >= 0.945)[1]])
}
joint_draw <- function(g, n, base = 5) {
  u <- halton(n, base); c_ <- cumsum(g$post)
  g[vapply(u, function(uu) which(c_ >= uu)[1], integer(1)), ]
}
wquantile <- function(v, w, p) {                 # 加权分位数（线性插值），三语言同法
  o <- order(v); v <- v[o]; cw <- cumsum(w[o]) / sum(w)
  k <- which(cw >= p)[1]
  if (k == 1) return(v[1])
  cw0 <- cw[k - 1]; v0 <- v[k - 1]
  v0 + (v[k] - v0) * (p - cw0) / (cw[k] - cw0)
}

o_bio <- ols(bio$时间, bio$厚度)
g_bio <- grid_post(bio$时间, bio$厚度, sigma = o_bio$sigma_hat)
o_hw  <- ols(d$height_c, d$weight)
g_hw  <- grid_post(d$height_c, d$weight, sigma = o_hw$sigma_hat, a_lim = c(30, 60), b_lim = c(0.4, 0.9))

iv_bio <- interval89(marg(g_bio, "b")); iv_hw <- interval89(marg(g_hw, "b"))
mean_a_bio <- weighted.mean(g_bio$a, g_bio$post); mean_b_bio <- weighted.mean(g_bio$b, g_bio$post)
mean_a_hw  <- weighted.mean(g_hw$a, g_hw$post);  mean_b_hw  <- weighted.mean(g_hw$b, g_hw$post)

summ <- tibble(
  案例 = c("我们的案例：厚度~时间", "教材对照：体重~身高"),
  最小二乘斜率 = round(c(o_bio$b, o_hw$b), 4),
  后验均值斜率 = round(c(mean_b_bio, mean_b_hw), 4),
  后验均值截距 = round(c(mean_a_bio, mean_a_hw), 4),
  斜率89下 = round(c(iv_bio[1], iv_hw[1]), 4),
  斜率89上 = round(c(iv_bio[2], iv_hw[2]), 4),
  残差标准差 = round(c(o_bio$sigma_hat, o_hw$sigma_hat), 4)
)
write_csv_fixed(summ, file.path(RES, "02_后验摘要.csv"))
print(summ)

# ---------- 3) 小波动相加：二项的概率质量 vs 正态的"整点区间概率" ----------
n_bin <- 10; ks <- 0:n_bin
sd_bin <- sqrt(n_bin * 0.25)
binom_df <- tibble(成功数 = ks, 概率 = dbinom(ks, n_bin, 0.5))
norm_df <- tibble(成功数 = ks,
                  区间概率 = pnorm(ks + 0.5, n_bin / 2, sd_bin) - pnorm(ks - 0.5, n_bin / 2, sd_bin))
write_csv_fixed(binom_df, file.path(RES, "03_二项概率.csv"))
write_csv_fixed(norm_df, file.path(RES, "03_正态的整点区间概率.csv"))

# ---------- 4) 先验预测：从"分析里真正用的那个平坦先验"抽线 ----------
n_prior <- 60
prior_lines <- tibble(线号 = 1:n_prior,
                      截距 = A_LIM[1] + (A_LIM[2] - A_LIM[1]) * halton(n_prior, 2),
                      斜率 = B_LIM[1] + (B_LIM[2] - B_LIM[1]) * halton(n_prior, 3))
write_csv_fixed(data.frame(线号 = prior_lines$线号,
                           截距 = round(prior_lines$截距, 6), 斜率 = round(prior_lines$斜率, 6)),
                file.path(RES, "04_先验预测的直线.csv"))

# ---------- 5) 后验预测（时间点 1–12；均值带与单次观测带分开）----------
n_draw <- 2000
dr <- joint_draw(g_bio, n_draw, 5)
eps <- o_bio$sigma_hat * qnorm(halton(n_draw, 11))
pred_summary <- tibble(时间 = t_hours) %>%
  mutate(预测均值 = sapply(t_hours, function(tt) mean(dr$a + dr$b * tt)),
         均值89下 = sapply(t_hours, function(tt) wquantile(dr$a + dr$b * tt, rep(1 / n_draw, n_draw), 0.055)),
         均值89上 = sapply(t_hours, function(tt) wquantile(dr$a + dr$b * tt, rep(1 / n_draw, n_draw), 0.945)),
         单次89下 = sapply(t_hours, function(tt) wquantile(dr$a + dr$b * tt + eps, rep(1 / n_draw, n_draw), 0.055)),
         单次89上 = sapply(t_hours, function(tt) wquantile(dr$a + dr$b * tt + eps, rep(1 / n_draw, n_draw), 0.945)))
write_csv_fixed(data.frame(时间 = pred_summary$时间,
                           预测均值 = round(pred_summary$预测均值, 3),
                           均值89下 = round(pred_summary$均值89下, 3),
                           均值89上 = round(pred_summary$均值89上, 3),
                           单次89下 = round(pred_summary$单次89下, 3),
                           单次89上 = round(pred_summary$单次89上, 3)),
                file.path(RES, "05_后验预测.csv"))

# ---------- 6) 先验敏感性（只动斜率先验，截距保持平坦；参数框不变）----------
sens <- list(list("平坦先验", NULL), list("斜率先验 N(3.4, 0.2)", c(3.4, 0.2)),
             list("斜率先验 N(5.0, 1.2)", c(5.0, 1.2)))
sens_df <- do.call(rbind, lapply(sens, function(s) {
  gg <- grid_post(bio$时间, bio$厚度, sigma = o_bio$sigma_hat, prior_b = s[[2]])
  iv <- interval89(marg(gg, "b"))
  data.frame(斜率先验 = s[[1]], 斜率均值 = round(weighted.mean(gg$b, gg$post), 4),
             斜率89下 = round(iv[1], 4), 斜率89上 = round(iv[2], 4),
             先验权重 = round(o_bio$se_b^2 / (o_bio$se_b^2 + if (is.null(s[[2]])) Inf else s[[2]][2]^2), 4))
}))
write_csv_fixed(sens_df, file.path(RES, "06_先验敏感性.csv"))
print(sens_df)

# ---------- 7) 网格分辨率：同一份后验，换网格密度 ----------
res_grid <- c(20, 50, 100, 200)
res_df <- do.call(rbind, lapply(res_grid, function(ng) {
  gg <- grid_post(bio$时间, bio$厚度, n_grid = ng, sigma = o_bio$sigma_hat)
  data.frame(网格点数 = ng, 后验均值斜率 = round(weighted.mean(gg$b, gg$post), 6),
             后验密度最高点斜率 = round(gg$b[which.max(gg$post)], 4))
}))
res_df$后验均值斜率 <- round(res_df$后验均值斜率, 4)
res_df$最小二乘参照 <- round(o_bio$b, 4)
write_csv_fixed(res_df, file.path(RES, "07_网格分辨率.csv"))
print(res_df)

# ---------- 图 ----------
p1 <- ggplot(bio, aes(时间, 厚度)) +
  geom_abline(data = prior_lines, aes(intercept = 截距, slope = 斜率),
              colour = "#167d80", alpha = 0.20, linewidth = 0.4) +
  geom_point(size = 2.6, colour = "#c1462c") +
  geom_abline(intercept = mean_a_bio, slope = mean_b_bio, colour = "#163d48", linewidth = 1.1) +
  labs(title = "先验给出一个「直线框」，数据把它收成一条窄带",
       subtitle = "点：教学模拟的 12 个厚度值；细线：从分析所用先验里抽的 60 条线；粗线：后验均值",
       x = "培养时间（小时）", y = "生物膜厚度（微米）")
save_fig("01-我们的案例-先验与后验直线.png", p1)

p2 <- ggplot(d, aes(height, weight)) +
  geom_point(size = 1.1, colour = "#9aa8a8") +
  geom_abline(intercept = mean_a_hw - mean_b_hw * mean(d$height), slope = mean_b_hw,
              colour = "#c1462c", linewidth = 1.1) +
  labs(title = "教材对照：成年人身高与体重",
       subtitle = paste0("Howell1 成年人 ", nrow(d), " 人；红线是后验均值直线（斜率 ",
                         round(mean_b_hw, 3), "，身高已中心化）"),
       x = "身高（厘米）", y = "体重（公斤）")
save_fig("02-教材对照-身高与体重.png", p2)

p3 <- ggplot() +
  geom_col(data = binom_df, aes(成功数, 概率), fill = "#167d80", width = 0.55, alpha = 0.85) +
  geom_line(data = norm_df, aes(成功数, 区间概率), colour = "#c1462c", linewidth = 1.1) +
  geom_point(data = norm_df, aes(成功数, 区间概率), colour = "#c1462c", size = 2.2) +
  labs(title = "十个小波动加起来，就近似钟形了",
       subtitle = "柱：10 次伯努利求和（p=0.5）的概率；红线：正态分布在每个整点区间的概率（同一刻度）",
       x = "成功次数（0 到 10）", y = "概率")
save_fig("03-小波动相加.png", p3)

p4 <- ggplot(res_df, aes(网格点数, 后验均值斜率)) +
  geom_hline(aes(yintercept = 最小二乘参照), colour = "#c1462c", linetype = "dashed", linewidth = 0.8) +
  geom_line(colour = "#167d80", linewidth = 0.9) +
  geom_point(colour = "#163d48", size = 3) +
  geom_text(aes(label = sprintf("%.4f", 后验均值斜率)), vjust = -0.9, size = 3.2) +
  scale_x_continuous(breaks = res_grid) +
  labs(title = "网格越细，后验均值越稳（虚线：最小二乘解 3.4224）",
       subtitle = "同一个后验、同一份数据，只换网格密度：20 → 200 点时估计如何收敛",
       x = "每个参数的网格点数", y = "斜率后验均值（微米/小时）")
save_fig("04-网格分辨率.png", p4)

p5 <- ggplot(bio, aes(时间, 厚度)) +
  geom_ribbon(data = pred_summary, aes(x = 时间, ymin = 单次89下, ymax = 单次89上),
              inherit.aes = FALSE, fill = "#167d80", alpha = 0.22) +
  geom_ribbon(data = pred_summary, aes(x = 时间, ymin = 均值89下, ymax = 均值89上),
              inherit.aes = FALSE, fill = "#163d48", alpha = 0.35) +
  geom_point(size = 2.4, colour = "#c1462c") +
  geom_line(data = pred_summary, aes(时间, 预测均值), inherit.aes = FALSE,
            colour = "#163d48", linewidth = 1) +
  labs(title = "同一个时间点，「平均厚度」与「再测一次的读数」",
       subtitle = "深线：预测均值；深带：平均厚度的 89% 区间；浅带：单次新观测的 89% 预测区间",
       x = "培养时间（小时）", y = "生物膜厚度（微米）")
save_fig("05-后验预测.png", p5)

cat("\n完成：5 张图 + 8 份 CSV 已写入\n")
