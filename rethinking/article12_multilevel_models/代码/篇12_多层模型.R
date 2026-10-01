# ============================================================
# 精读 12｜多层模型：为什么"分组各算"和"全部混算"都不对
# 对应：McElreath 2023 第 12 讲 Multilevel Models；教材第 13 章 Models With Memory
# 用法： Rscript 代码/篇12_多层模型.R
#
# 模型：y_ij = α_j + b·x_ij + ε_ij，α_j ~ Normal(μ, σ_α)（随机截距）
# 计算：对 (b, σ_α, σ) 铺三维网格；α_j 与 μ 解析积分掉，得到闭式边际似然：
#   瓶间部分（组均值）→ 提供 σ_α 与 σ 的信息；瓶内部分（组内偏离）→ 提供斜率 b 的信息
# 这样网格只有三维、完全确定性，三语言可以逐字节一致。
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
r4 <- function(x) round(as.numeric(x), 4)
wq <- function(v, w, p) {
  o <- order(v); v <- v[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  k <- which(cw >= p)[1]
  if (k == 1) return(v[1])
  v[k - 1] + (v[k] - v[k - 1]) * (p - cw[k - 1]) / (cw[k] - cw[k - 1])
}

# ---------- 1) 数据：12 个培养瓶 × 每瓶 3 个孔 ----------
x <- rep(c(1, 2, 3), 12)
# 目标：瓶间差异 σ_α ≈ 0.7、瓶内噪声 σ ≈ 0.55。这样 σ_α² 远大于 σ²/n（≈0.10），
# 部分合并的收缩看得见，且"新瓶预测"比"同瓶预测"更宽（这是本篇要讲的第二件事）。
alpha_true <- c(4.60, 4.90, 5.15, 5.35, 5.50, 5.60, 5.70, 5.85, 6.05, 6.30, 6.60, 7.00)
b_true <- 1.0
# 手写残差（不用随机数）：每瓶三个孔的和为 0，四种模式轮换，使残差与瓶内剂量梯度尽量不相关
resid <- c(-0.66, 0.11, 0.55,   0.55, -0.66, 0.11,   -0.26, 0.78, -0.52,
            0.39, -0.65, 0.26,  -0.66, 0.11, 0.55,    0.55, -0.66, 0.11,
           -0.26, 0.78, -0.52,   0.39, -0.65, 0.26,   -0.66, 0.11, 0.55,
            0.55, -0.66, 0.11,  -0.26, 0.78, -0.52,    0.39, -0.65, 0.26)
grid_id <- rep(1:12, each = 3)
y <- alpha_true[grid_id] + b_true * x + resid
write_csv_fixed(data.frame(瓶 = grid_id, 剂量 = x, 长势 = y), file.path(RES, "01_数据.csv"))

n_j <- 3; J <- 12; N <- length(y)
gmean <- as.numeric(tapply(y, grid_id, mean))
gxbar <- mean(x); gbar <- mean(y)
within_y <- y - ave(y, grid_id)          # 瓶内偏离（不含 b）
within_x <- x - ave(x, grid_id)

# ---------- 2) 三种算法 ----------
Xp <- cbind(1, x)
bp <- solve(t(Xp) %*% Xp, t(Xp) %*% y)
b_fe <- sum(within_y * within_x) / sum(within_x^2)     # 瓶内（固定效应）斜率
pooled_level <- mean(y - bp[2] * x)
unpooled <- gmean - b_fe * gxbar                       # 各瓶自己的水平
pooled <- rep(mean(unpooled), J)                       # 完全合并且忽略瓶

# 边际似然（μ 与各 α_j 解析积分掉）
marg_ll <- function(b, sa, s) {
  r <- gmean - b * gxbar
  v <- rep(s^2 / n_j + sa^2, J)     # 每瓶一个 v（本例每瓶孔数相同，但必须是长度为 J 的向量）
  S <- sum(1 / v)
  mu_hat <- sum(r / v) / S
  q_between <- sum((r - mu_hat)^2 / v)                 # 只与瓶间差异有关
  w <- within_y - b * within_x
  -0.5 * (log(S) + sum(log(v)) + q_between) - 0.5 * sum(w^2) / s^2 - (N - J) * log(s)
}
b_seq  <- seq(0, 2, length.out = 60)
sa_seq <- seq(0.005, 2, length.out = 60)
s_seq  <- seq(0.05, 2, length.out = 60)
gg <- expand.grid(b = b_seq, sa = sa_seq, s = s_seq)
lp <- numeric(nrow(gg))
for (i in seq_len(nrow(gg))) lp[i] <- marg_ll(gg$b[i], gg$sa[i], gg$s[i])
gg$post <- exp(lp - max(lp)); gg$post <- gg$post / sum(gg$post)

shrink_f <- function(b, sa, s) sa^2 / (sa^2 + s^2 / n_j)
mu_of <- function(b) mean(gmean) - b * gxbar      # μ̂ 随 b 线性变化（这里是向量化写法）
mu_draw <- mu_of(gg$b)
alpha_draw <- sapply(seq_len(nrow(gg)), function(i)
  mu_draw[i] + (gmean - gg$b[i] * gxbar - mu_draw[i]) * shrink_f(gg$b[i], gg$sa[i], gg$s[i]))
partial <- as.numeric(alpha_draw %*% gg$post)
b_post <- weighted.mean(gg$b, gg$post)
mu_post <- weighted.mean(mu_draw, gg$post)
sa_post <- weighted.mean(gg$sa, gg$post)
s_post <- weighted.mean(gg$s, gg$post)
shrink_post <- weighted.mean(sapply(seq_len(nrow(gg)), function(i) shrink_f(gg$b[i], gg$sa[i], gg$s[i])), gg$post)

cmp <- data.frame(瓶 = 1:12, 瓶内均值 = r4(unpooled), 完全合并 = r4(pooled),
                  完全不合并 = r4(unpooled), 部分合并 = r4(partial),
                  收缩量 = r4(unpooled - partial))
write_csv_fixed(cmp, file.path(RES, "02_三种算法.csv"))
print(cmp)

# 用字符列写 "NA"，避免 R 的 formatC 把 NA 补成带空格的 "     NA"（三语言要不一致）
post_tab <- data.frame(
  参数 = c("共同斜率 b", "总体水平 μ", "组间标准差 σ_α", "残差标准差 σ", "收缩因子（每瓶都是 3 个孔）"),
  后验均值 = num(r4(c(b_post, mu_post, sa_post, s_post, shrink_post))),
  下界89 = c(num(r4(wq(gg$b, gg$post, 0.055))), "NA", num(r4(wq(gg$sa, gg$post, 0.055))),
             num(r4(wq(gg$s, gg$post, 0.055))), "NA"),
  上界89 = c(num(r4(wq(gg$b, gg$post, 0.945))), "NA", num(r4(wq(gg$sa, gg$post, 0.945))),
             num(r4(wq(gg$s, gg$post, 0.945))), "NA")
)
write_csv_fixed(post_tab, file.path(RES, "03_超参数后验.csv"))
print(post_tab)
cat("网格：60³ =", nrow(gg), "个格点\n")

# ---------- 3) 两种不确定性的来源 ----------
se_existing <- sqrt(s_post^2 + s_post^2 / n_j)
se_new <- sqrt(s_post^2 + sa_post^2)
pred_tab <- data.frame(
  情形 = c("同一个瓶再测一个孔", "换一个全新的瓶再测一个孔"),
  预测标准差 = r4(c(se_existing, se_new)),
  区间半宽89 = r4(c(1.598 * se_existing, 1.598 * se_new)),
  说明 = c("只多了观测误差与瓶均值自身的不确定性", "还要加上瓶与瓶之间的真实差异 σ_α")
)
write_csv_fixed(pred_tab, file.path(RES, "04_两种不确定性.csv"))
print(pred_tab)

write_csv_fixed(data.frame(
  量 = c("完全合并的水平", "完全不合并的极差", "部分合并的极差",
         "收缩量与偏离的相关", "瓶间真实差异（生成时设定）", "瓶内噪声（生成时设定）"),
  值 = r4(c(mean(unpooled), diff(range(unpooled)), diff(range(partial)),
             cor(unpooled - mean(unpooled), unpooled - partial),
             sd(alpha_true), sqrt(mean(resid^2))))
), file.path(RES, "05_汇总.csv"))

# ---------- 图 ----------
p1 <- ggplot(data.frame(瓶 = factor(grid_id), 长势 = y), aes(瓶, 长势)) +
  geom_jitter(width = 0.12, size = 2.2, colour = "#167d80", alpha = 0.85) +
  geom_point(data = data.frame(瓶 = factor(1:12), 长势 = gmean), colour = "#c1462c", size = 3) +
  geom_hline(yintercept = gbar, linetype = "dashed", colour = "#9aa8a8") +
  labs(title = "12 个培养瓶，每瓶 3 个孔",
       subtitle = "红点＝瓶内均值；虚线＝整体均值。瓶与瓶确实有差异，但差异里也混着噪声",
       x = "培养瓶", y = "长势")
save_fig("01_数据.png", p1, 7.4, 4.4)

p2 <- ggplot(cmp, aes(完全不合并, 部分合并)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9aa8a8") +
  geom_hline(yintercept = mean(unpooled), linetype = "dotted", colour = "#163d48") +
  geom_segment(aes(xend = 完全不合并, yend = 完全不合并), colour = "#c1462c",
               arrow = arrow(length = unit(0.12, "cm"))) +
  geom_point(size = 3, colour = "#167d80") +
  labs(title = "部分合并＝把每个瓶的估计往整体拉一点",
       subtitle = paste0("虚线＝不收缩；红箭头＝收缩量；点线＝所有瓶共同的水平 ", r4(mean(unpooled))),
       x = "完全不合并（瓶内均值）", y = "部分合并（多层模型）")
save_fig("02_收缩.png", p2, 7.0, 5.0)

shrink_df <- data.frame(偏离 = unpooled - mean(unpooled), 收缩量 = unpooled - partial)
p3 <- ggplot(shrink_df, aes(偏离, 收缩量)) +
  geom_point(size = 3, colour = "#167d80") +
  geom_smooth(method = "lm", se = FALSE, colour = "#c1462c", linewidth = 1) +
  labs(title = "偏离整体越远，被拉回的越多",
       subtitle = "收缩量与偏离成正比——这就是部分合并的机制（每瓶都是 3 个孔，收缩比例相同）",
       x = "瓶内均值 − 整体水平", y = "收缩量")
save_fig("03_收缩规律.png", p3, 7.4, 4.4)

dens <- rbind(
  data.frame(参数 = "共同斜率 b", 值 = gg$b, 权重 = gg$post),
  data.frame(参数 = "组间标准差 σ_α", 值 = gg$sa, 权重 = gg$post),
  data.frame(参数 = "残差标准差 σ", 值 = gg$s, 权重 = gg$post)
)
dens_sum <- dens %>% group_by(参数) %>%
  summarise(均值 = weighted.mean(值, 权重), 下 = wq(值, 权重, 0.055), 上 = wq(值, 权重, 0.945), .groups = "drop")
p4 <- ggplot(dens_sum, aes(参数, 均值)) +
  geom_point(size = 3, colour = "#167d80") +
  geom_errorbar(aes(ymin = 下, ymax = 上), width = 0.15, colour = "#163d48") +
  geom_text(aes(label = sprintf("%.3f", 均值)), vjust = -0.8, size = 3.8) +
  labs(title = "三个参数一起估：斜率与两个标准差",
       subtitle = "σ_α 是瓶与瓶之间的真实差异，σ 是同一瓶内的观测噪声——两者要靠数据分开",
       x = NULL, y = "后验均值（89% 区间）")
save_fig("04_超参数.png", p4, 7.4, 4.2)

p5 <- ggplot(pred_tab, aes(情形, 区间半宽89, fill = 情形)) +
  geom_col(width = 0.45) +
  geom_text(aes(label = sprintf("%.3f", 区间半宽89)), vjust = -0.6, size = 4) +
  scale_fill_manual(values = c("#e8b4a8", "#167d80")) +
  labs(title = "预测全新瓶，比预测同一瓶的新孔不确定得多",
       subtitle = paste0("只多了一层 σ_α = ", r4(sa_post), "，预测标准差就从 ",
                         r4(se_existing), " 涨到 ", r4(se_new)),
       x = NULL, y = "89% 区间半宽") +
  theme(legend.position = "none")
save_fig("05_预测区间.png", p5, 7.4, 4.2)

cat("\n完成：5 张图 + 5 份 CSV 已写入\n")
