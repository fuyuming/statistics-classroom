# ============================================================
# 精读 07｜过拟合与正则化：样本内好不等于样本外好
# 对应：McElreath 2023 第 07 讲 Overfitting；教材第 7、8 章
# 用法： Rscript 代码/篇07_过拟合.R
#
# 约定：与 05、06 篇同一套工程约定——显式最小二乘（解正规方程）、
#       留一交叉验证用帽子矩阵的闭式解（确定性）、厚尾似然自己按 t 密度公式算。
# ============================================================
suppressPackageStartupMessages({ library(tidyverse) })

.args <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .args, value = TRUE)
root <- if (length(.f) > 0) dirname(dirname(normalizePath(sub("^--file=", "", .f[1])))) else getwd()
OUT <- file.path(root, "文章配图"); RES <- file.path(root, "运行结果")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(name, plot, w = 7.4, h = 4.4) {
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
ols_ls <- function(X, y) {
  b <- solve(t(X) %*% X, t(X) %*% y)
  e <- as.numeric(y - X %*% b)
  list(coef = as.numeric(b), resid = e,
       sigma_hat = sqrt(sum(e^2) / (length(y) - ncol(X))),
       sse = sum(e^2), r2 = 1 - sum(e^2) / sum((y - mean(y))^2))
}
# 留一交叉验证（最小二乘闭式解）：把 e_i 除以 (1 - h_ii)
loo_error <- function(X, y) {
  XtX_inv <- solve(t(X) %*% X)
  h <- rowSums((X %*% XtX_inv) * X)
  fit <- ols_ls(X, y)
  loo_i <- fit$resid / (1 - h)
  sqrt(mean(loo_i^2))
}
orth_resid <- function(v, X) v - X %*% solve(t(X) %*% X, t(X) %*% v)
r4 <- function(x) round(as.numeric(x), 4)

# ---------- 1) 数据：12 个观测，生成均值为一个三次曲线 ----------
n <- 12
t_x <- 1:n
X0 <- cbind(1, t_x)
# 多项式基用标准化后的时间轴（否则 t^8 之类的量级会让正规方程数值奇异）
t_s <- (t_x - mean(t_x)) / sd(t_x)
ts_fine_of <- function(tv) (tv - mean(t_x)) / sd(t_x)
# 生成均值：3 阶形状（真值），残差与 [1, t, t², t³] 正交
mu_true <- 20 + 4 * t_x - 0.35 * t_x^2 + 0.010 * t_x^3
resid_g <- orth_resid(c(0.8, -1.2, 0.6, -0.9, 1.1, -0.5, 0.7, -1.0, 0.4, 0.9, -0.6, 1.0),
                      cbind(1, t_x, t_x^2, t_x^3))
y_obs <- mu_true + resid_g
write_csv_fixed(data.frame(时间 = t_x, 观测 = r4(y_obs), 生成均值 = r4(mu_true)),
                file.path(RES, "01_数据.csv"))

# ---------- 2) 多项式阶数 1–8：样本内 vs 留一交叉验证 ----------
poly_X <- function(k) {
  out <- cbind(1, t_s)
  if (k >= 2) for (j in 2:k) out <- cbind(out, t_s^j)
  out
}
rows <- lapply(1:7, function(k) {
  X <- poly_X(k)
  f <- ols_ls(X, y_obs)
  data.frame(阶数 = k, 系数个数 = k + 1,
             样本内残差标准差 = r4(f$sigma_hat),
             样本内平方误差和 = r4(f$sse),
             留一交叉验证误差 = r4(loo_error(X, y_obs)))
})
df_res <- do.call(rbind, rows)
write_csv_fixed(df_res, file.path(RES, "02_阶数与误差.csv"))
print(df_res)
best_k <- df_res$阶数[which.min(df_res$留一交叉验证误差)]
cat("留一交叉验证最优阶数 =", best_k, "\n")

# ---------- 3) 过拟合：3 阶与 8 阶的拟合曲线（含外推） ----------
t_fine <- seq(1, 16, length.out = 200)
truth_fine <- 20 + 4 * t_fine - 0.35 * t_fine^2 + 0.010 * t_fine^3
fit_curve <- function(k, y = y_obs) {
  X <- poly_X(k)
  b <- ols_ls(X, y)$coef
  tsf <- ts_fine_of(t_fine)
  Xf <- cbind(1, tsf)
  if (k >= 2) for (j in 2:k) Xf <- cbind(Xf, tsf^j)
  as.numeric(Xf %*% b)
}
curve_df <- tibble(时间 = t_fine, `3 阶` = fit_curve(3), `7 阶` = fit_curve(7),
                   `真值形状` = truth_fine) %>%
  pivot_longer(-时间, names_to = "模型", values_to = "拟合")
curve_df$模型 <- factor(curve_df$模型, levels = c("真值形状", "3 阶", "7 阶"))
write_csv_fixed(data.frame(时间 = r4(t_fine), 真值形状 = r4(truth_fine),
                           三阶 = r4(fit_curve(3)), 七阶 = r4(fit_curve(7))),
                file.path(RES, "03_两条拟合曲线.csv"))

# ---------- 4) 正则化＝先验宽度：同一模型三种先验 ----------
# 用一阶模型（截距＋斜率），把斜率先验从窄到宽改三次；σ 固定为残差标准差
X1 <- cbind(1, t_x)
s1 <- ols_ls(X1, y_obs)$sigma_hat
grid_post1 <- function(prior_b) {
  a_seq <- 15 + (0:199) * (30 / 199); b_seq <- 0 + (0:199) * (8 / 199)
  g <- expand.grid(a = a_seq, b = b_seq)
  lp <- numeric(nrow(g))
  for (i in seq_len(nrow(g))) {
    mu <- g$a[i] + g$b[i] * t_x
    lp[i] <- -0.5 * sum(((y_obs - mu) / s1)^2) - n * (log(s1) + 0.5 * log(2 * pi))
    if (!is.null(prior_b)) lp[i] <- lp[i] - 0.5 * ((g$b[i] - prior_b[1]) / prior_b[2])^2
  }
  g$post <- exp(lp - max(lp)); g$post <- g$post / sum(g$post)
  list(g = g, a_seq = a_seq, b_seq = b_seq)
}
b_mean <- function(res) weighted.mean(res$g$b, res$g$post)
b_sd   <- function(res) sqrt(weighted.mean((res$g$b - weighted.mean(res$g$b, res$g$post))^2, res$g$post))
prior_list <- list(c(3.0, 0.2), c(3.0, 1.0), NULL)
prior_names <- c("窄先验 N(3, 0.2)", "中等先验 N(3, 1.0)", "平坦先验")
reg_rows <- do.call(rbind, lapply(seq_along(prior_list), function(i) {
  res <- grid_post1(prior_list[[i]])
  data.frame(先验 = prior_names[i], 后验均值斜率 = r4(b_mean(res)), 后验标准差 = r4(b_sd(res)),
             预测第16期 = r4(weighted.mean(res$g$a + res$g$b * 16, res$g$post)))
}))
write_csv_fixed(reg_rows, file.path(RES, "04_先验宽度.csv"))
print(reg_rows)

# ---------- 5) 厚尾：一个异常点，正态 vs t(df=4) ----------
y_out2 <- y_obs
y_out2[9] <- y_out2[9] + 22                       # 第 9 个观测被污染
X3 <- poly_X(3)
f_norm <- ols_ls(X3, y_out2)
# t 似然的权重（IRLS 一步近似，确定性）：w_i = (df+1)/(df + z_i^2)，z_i 用正态残差标准化
z <- f_norm$resid / f_norm$sigma_hat
w_t <- (4 + 1) / (4 + z^2)
f_t <- ols_ls(X3 * sqrt(w_t), y_out2 * sqrt(w_t))
robust_df <- data.frame(
  观测 = 1:12,
  标准化残差 = r4(z),
  厚尾权重 = r4(w_t)
)
write_csv_fixed(robust_df, file.path(RES, "05_厚尾权重.csv"))
coef_df <- data.frame(
  模型 = c("正态最小二乘", "厚尾（t, df=4，IRLS 一步）"),
  三次项系数 = r4(c(f_norm$coef[4], f_t$coef[4])),
  残差标准差 = r4(c(f_norm$sigma_hat, f_t$sigma_hat))
)
write_csv_fixed(coef_df, file.path(RES, "06_正态与厚尾.csv"))
print(coef_df)

# ---------- 图 ----------
p1 <- ggplot(df_res, aes(阶数)) +
  geom_line(aes(y = 样本内平方误差和, colour = "样本内平方误差和（越小越好骗）"), linewidth = 1) +
  geom_point(aes(y = 样本内平方误差和, colour = "样本内平方误差和（越小越好骗）"), size = 2.6) +
  geom_line(aes(y = 留一交叉验证误差^2 * 12, colour = "留一交叉验证误差（放大到同一尺度）"), linewidth = 1) +
  geom_point(aes(y = 留一交叉验证误差^2 * 12, colour = "留一交叉验证误差（放大到同一尺度）"), size = 2.6) +
  geom_vline(xintercept = best_k, linetype = "dashed", colour = "#9aa8a8") +
  labs(title = "样本内误差一路变小，样本外误差却先降后升",
       subtitle = paste0("虚线是最优阶数 k = ", best_k, "（按留一交叉验证）"),
       x = "多项式阶数", y = "误差", colour = NULL)
save_fig("01_阶数与误差.png", p1, 7.6, 4.4)

p2 <- ggplot() +
  geom_point(data = data.frame(t_x = t_x, y = y_obs), aes(t_x, y), size = 2.6, colour = "#c1462c") +
  geom_line(data = curve_df, aes(时间, 拟合, colour = 模型), linewidth = 1) +
  geom_vline(xintercept = 12, linetype = "dashed", colour = "#9aa8a8") +
  labs(title = "3 阶与 7 阶：样本内都贴得住，出了数据范围就分道扬镳",
       subtitle = "虚线右侧没有观测；高阶曲线在那里迅速偏离生成它的形状",
       x = "时间", y = "观测值", colour = NULL)
save_fig("02_两阶对照.png", p2, 7.6, 4.4)

dens_df <- do.call(rbind, lapply(seq_along(prior_list), function(i) {
  # 用后验标准差画示意密度（横轴围绕各自均值）
  res <- grid_post1(prior_list[[i]])
  m <- b_mean(res); s <- b_sd(res)
  x <- m + seq(-4, 4, length.out = 300) * s
  data.frame(斜率 = x, 密度 = dnorm(x, m, s), 先验 = prior_names[i])
}))
dens_df$先验 <- factor(dens_df$先验, levels = prior_names)
p3 <- ggplot(dens_df, aes(斜率, 密度, colour = 先验)) +
  geom_line(linewidth = 1.1) +
  labs(title = "正则化就是先验宽度：窄先验把斜率拉回来",
       subtitle = "同一条直线、同一份数据，只改斜率先验的宽度",
       x = "斜率", y = "后验密度", colour = NULL)
save_fig("03_先验宽度.png", p3, 7.6, 4.4)

p4 <- ggplot(robust_df, aes(factor(观测), 厚尾权重)) +
  geom_col(fill = "#167d80", width = 0.6) +
  geom_text(aes(label = ifelse(厚尾权重 < 0.5, "被压低的异常点", "")), vjust = -0.5, size = 3.4, colour = "#c1462c") +
  labs(title = "厚尾似然自动压低异常点的权重",
       subtitle = "第 9 个观测被污染：它的厚尾权重明显低于其他观测（正态模型里它的权重始终是 1）",
       x = "观测编号", y = "厚尾权重") +
  ylim(0, 1.35)
save_fig("04_厚尾权重.png", p4, 7.6, 4.4)

p5 <- ggplot() +
  geom_point(data = data.frame(t_x = t_x, y = y_out2), aes(t_x, y), size = 2.6, colour = "#c1462c") +
  geom_line(data = data.frame(t_fine = t_fine,
                              干净数据拟合 = fit_curve(3, y_obs),
                              含异常点拟合 = fit_curve(3, y_out2)) %>%
              pivot_longer(-t_fine, names_to = "模型", values_to = "拟合"),
            aes(t_fine, 拟合, colour = 模型), linewidth = 1) +
  labs(title = "一个异常点就把三次曲线拽走",
       subtitle = "第 9 个观测被抬高 22：含它的拟合（实线之一）与不含它的拟合明显分开",
       x = "时间", y = "观测值", colour = NULL)
save_fig("05_异常点影响.png", p5, 7.6, 4.4)

cat("\n完成：5 张图 + 6 份 CSV 已写入\n")
