# ============================================================
# 精读 09｜事件数据：0/1 结果、计数、率与有序类别
# 对应：McElreath 2023 第 09 讲 Modeling Events；教材第 10、11 章
# 用法： Rscript 代码/篇09_事件数据.R
#
# 约定：与前几篇同一套工程约定——确定性网格后验／显式公式，三语言逐位一致。
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
# 网格后验工具：log 似然函数 + 两参数网格
grid2 <- function(a_seq, b_seq, loglik) {
  g <- expand.grid(a = a_seq, b = b_seq)
  lp <- numeric(nrow(g))
  for (i in seq_len(nrow(g))) lp[i] <- loglik(g$a[i], g$b[i])
  p <- exp(lp - max(lp)); g$post <- p / sum(p)
  g
}

# ---------- 1) 0/1 结果：定植成功与否 ----------
dose <- c(1, 2, 3, 4, 5, 6, 7, 8)
n_units <- rep(10, 8)
n_succ <- c(1, 2, 4, 5, 7, 8, 9, 10)
binom_df <- data.frame(剂量 = dose, 单元数 = n_units, 成功数 = n_succ,
                       成功比例 = r4(n_succ / n_units))
write_csv_fixed(binom_df, file.path(RES, "01_定植数据.csv"))

loglik_binom <- function(a, b) {
  p <- 1 / (1 + exp(-(a + b * dose)))
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  sum(n_succ * log(p) + (n_units - n_succ) * log(1 - p))
}
# 斜率网格要包含负值（否则"平坦先验"其实带着一个非负斜率的先验限制）
g_binom <- grid2(seq(-8, 3, length.out = 300), seq(-1.5, 2.5, length.out = 300), loglik_binom)
binom_sum <- data.frame(
  参数 = c("截距 a", "斜率 b"),
  后验均值 = r4(c(weighted.mean(g_binom$a, g_binom$post), weighted.mean(g_binom$b, g_binom$post))),
  下界89 = r4(c(wq(g_binom$a, g_binom$post, 0.055), wq(g_binom$b, g_binom$post, 0.055))),
  上界89 = r4(c(wq(g_binom$a, g_binom$post, 0.945), wq(g_binom$b, g_binom$post, 0.945)))
)
write_csv_fixed(binom_sum, file.path(RES, "02_二项后验.csv"))
print(binom_sum)

# 预测曲线：每个剂量点的成功概率后验
dose_fine <- seq(0, 9, length.out = 60)
pred_curve <- sapply(dose_fine, function(d) {
  p <- 1 / (1 + exp(-(g_binom$a + g_binom$b * d)))
  c(mean = weighted.mean(p, g_binom$post),
     lo = wq(p, g_binom$post, 0.055), hi = wq(p, g_binom$post, 0.945))
})
# 概率曲线按 3 位小数写出：三语言的浮点求和顺序不同，4 位小数会在进位边界上出现末位差
pred_df <- data.frame(剂量 = round(dose_fine, 3), 均值 = round(pred_curve["mean", ], 3),
                      下 = round(pred_curve["lo", ], 3), 上 = round(pred_curve["hi", ], 3))
write_csv_fixed(pred_df, file.path(RES, "03_定植概率曲线.csv"))
p50 <- weighted.mean(1 / (1 + exp(-(g_binom$a + g_binom$b * 5))), g_binom$post) %>% r4()
cat("剂量 5 的成功概率后验均值 =", p50, "\n")

# ---------- 2) 计数：菌落数（泊松）----------
# 每个剂量做 3 个重复：这样才能算"同一剂量内"的方差/均值，做真正的过离散检查
treat <- rep(0:7, each = 3)
counts <- c(1, 2, 3,      # 剂量 0，均值 2
            1, 3, 5,      # 剂量 1，均值 3
            2, 5, 8,      # 剂量 2，均值 5
            2, 6, 10,     # 剂量 3，均值 6
            4, 9, 14,     # 剂量 4，均值 9
            4, 11, 18,    # 剂量 5，均值 11
            5, 14, 23,    # 剂量 6，均值 14
            6, 18, 30)    # 剂量 7，均值 18
cnt_df <- data.frame(处理水平 = treat, 菌落数 = counts)
write_csv_fixed(cnt_df, file.path(RES, "04_计数数据.csv"))
mean_cnt <- mean(counts)
# 同一剂量内的方差与均值 → 合并的组内离散度指标（描述性，不是通用诊断参数）
within <- sapply(0:7, function(d) c(var(counts[treat == d]), mean(counts[treat == d])))
within_disp <- sum(within[1, ]) / sum(within[2, ])

loglik_pois <- function(a, b) {
  lam <- exp(a + b * treat)
  sum(counts * log(lam) - lam - lgamma(counts + 1))
}
g_pois <- grid2(seq(-1, 3, length.out = 300), seq(-1, 1, length.out = 300), loglik_pois)
pois_sum <- data.frame(
  参数 = c("截距 a", "斜率 b"),
  后验均值 = r4(c(weighted.mean(g_pois$a, g_pois$post), weighted.mean(g_pois$b, g_pois$post))),
  下界89 = r4(c(wq(g_pois$a, g_pois$post, 0.055), wq(g_pois$b, g_pois$post, 0.055))),
  上界89 = r4(c(wq(g_pois$a, g_pois$post, 0.945), wq(g_pois$b, g_pois$post, 0.945)))
)
write_csv_fixed(pois_sum, file.path(RES, "05_泊松后验.csv"))
print(pois_sum)
# 交叉核对：拟合泊松回归后的 Pearson 离散度（与上面的组内指标不是同一个量）
lam_fit <- exp(weighted.mean(g_pois$a, g_pois$post) + weighted.mean(g_pois$b, g_pois$post) * treat)
pearson_disp <- sum((counts - lam_fit)^2 / lam_fit) / (length(counts) - 2)
# 率比：处理水平 7 相对 0
rr <- exp(g_pois$b * 7)
rr_tab <- data.frame(量 = c("每皿平均菌落数（观测）", "方差/均值（同一剂量内合并，过离散检查）",
                            "率比 exp(7b) 后验均值", "率比 89% 区间下", "率比 89% 区间上"),
                     值 = r4(c(mean_cnt, within_disp, weighted.mean(rr, g_pois$post),
                               wq(rr, g_pois$post, 0.055), wq(rr, g_pois$post, 0.945))))
write_csv_fixed(rr_tab, file.path(RES, "06_率比与过离散.csv"))
print(rr_tab)
write_csv_fixed(data.frame(量 = c("同一剂量内合并的离散度（描述性）", "拟合泊松后的 Pearson 离散度"),
                           值 = r4(c(within_disp, pearson_disp))),
                file.path(RES, "07_两种离散度.csv"))

# ---------- 3) 率：暴露量不同（offset）----------
plate <- c(10, 10, 10, 10)          # 接种体积（微升）
obs_small <- c(3, 5, 4, 4)          # 小体积皿的菌落数
plate_big <- c(40, 40, 40, 40)      # 大体积皿
obs_big <- c(12, 18, 15, 19)
naive <- mean(obs_big) / mean(obs_small)                       # 直接比计数（错）
offset_correct <- (sum(obs_big) / sum(plate_big)) / (sum(obs_small) / sum(plate_small <- plate))
off_tab <- data.frame(
  比较 = c("直接比计数（忽略体积）", "按体积折算成率再比"),
  比值 = r4(c(naive, offset_correct)),
  说明 = c("大体积皿当然菌落更多，这是体积的功劳", "把体积放进模型（对数偏移），比的是每微升的密度")
)
write_csv_fixed(off_tab, file.path(RES, "07_offset对照.csv"))
print(off_tab)

# ---------- 4) 有序类别：轻/中/重 ----------
ord_x <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
ord_y <- c(0, 0, 0, 1, 0, 1, 1, 1, 1, 2, 2, 2)     # 0=轻 1=中 2=重
ord_df <- data.frame(剂量 = ord_x, 等级 = ord_y)
write_csv_fixed(ord_df, file.path(RES, "08_有序数据.csv"))
# 累积 logit：logit(F_k) = c_k − b·x，切点必须 c0 < c1
# （初版把切点写反了：F(≤轻) > F(≤中)，导致 P(中) 为负；已修正并去掉掩盖负值的 pmax）
# 本例暂把两个切点固定为 0 和 1.5，只估斜率——所以区间只反映"切点已知"时的不确定性
c0_fix <- 0; c1_fix <- 1.5
loglik_ord1 <- function(b) {
  F0 <- 1 / (1 + exp(-(c0_fix - b * ord_x)))     # P(Y ≤ 轻)
  F1 <- 1 / (1 + exp(-(c1_fix - b * ord_x)))     # P(Y ≤ 中)
  ll <- 0
  for (i in seq_along(ord_y)) {
    if (ord_y[i] == 0) ll <- ll + log(F0[i])
    else if (ord_y[i] == 1) ll <- ll + log(F1[i] - F0[i])
    else ll <- ll + log(1 - F1[i])
  }
  ll
}
b_seq <- seq(-2, 3, length.out = 400)
lp_ord <- sapply(b_seq, loglik_ord1)
w_ord <- exp(lp_ord - max(lp_ord)); w_ord <- w_ord / sum(w_ord)
ord_sum <- data.frame(
  参数 = "斜率 b",
  后验均值 = r4(weighted.mean(b_seq, w_ord)),
  下界89 = r4(wq(b_seq, w_ord, 0.055)),
  上界89 = r4(wq(b_seq, w_ord, 0.945))
)
write_csv_fixed(ord_sum, file.path(RES, "09_有序后验.csv"))
print(ord_sum)
cat("切点固定为 c0 =", c0_fix, "（轻/中之间的切点）、c1 =", c1_fix, "（中/重之间的切点）；只估斜率\n")
# 把等级当连续变量 vs 累积 logit：同一份数据的两种读法
cont_fit <- coef(lm(ord_y ~ ord_x))
# 概率类输出按 3 位小数写：三语言浮点求和顺序不同，4 位小数会在进位边界造成末位差
r3 <- function(x) round(as.numeric(x), 3)
ord_curve <- data.frame(
  剂量 = 1:6,
  轻的累积概率 = r3(sapply(1:6, function(d) weighted.mean(1 / (1 + exp(-(c0_fix - b_seq * d))), w_ord))),
  中的累积概率 = r3(sapply(1:6, function(d) weighted.mean(1 / (1 + exp(-(c1_fix - b_seq * d))), w_ord))),
  当连续变量_平均编码得分 = r3(cont_fit[1] + cont_fit[2] * (1:6))
)
ord_curve$重的概率 <- r3(1 - ord_curve$中的累积概率)
ord_curve$轻的类别概率 <- r3(ord_curve$轻的累积概率)
ord_curve$中的类别概率 <- r3(ord_curve$中的累积概率 - ord_curve$轻的累积概率)
write_csv_fixed(ord_curve, file.path(RES, "10_连续vs有序.csv"))
print(ord_curve)

# ---------- 图 ----------
p1 <- ggplot() +
  geom_point(data = binom_df, aes(剂量, 成功比例), size = 3, colour = "#c1462c") +
  geom_line(data = pred_df, aes(剂量, 均值), colour = "#167d80", linewidth = 1.1) +
  geom_ribbon(data = pred_df, aes(剂量, ymin = 下, ymax = 上), fill = "#167d80", alpha = 0.18) +
  labs(title = "0/1 结果要的是概率曲线，不是直线",
       subtitle = "点＝各剂量的成功比例；带＝成功概率的 89% 后验区间（logistic 联系）",
       x = "剂量", y = "成功概率") +
  coord_cartesian(ylim = c(0, 1))
save_fig("01_logit.png", p1, 7.4, 4.4)

p2 <- ggplot(cnt_df, aes(处理水平, 菌落数)) +
  geom_point(size = 3, colour = "#c1462c") +
  geom_line(data = data.frame(处理水平 = seq(0, 7, length.out = 60),
                              lam = exp(weighted.mean(g_pois$a, g_pois$post) +
                                          weighted.mean(g_pois$b, g_pois$post) * seq(0, 7, length.out = 60))),
            aes(处理水平, lam), colour = "#167d80", linewidth = 1.1) +
  labs(title = "计数数据：把线性预测放进指数里（泊松的 log 联系）",
       subtitle = paste0("曲线是后验均值下的预测均值；同一剂量内合并的离散度 = ", r4(within_disp),
                         "，明显大于 1 ——重复间的波动超过泊松预期"),
       x = "处理水平", y = "每皿菌落数")
save_fig("02_泊松.png", p2, 7.4, 4.4)

p3 <- ggplot(off_tab, aes(比较, 比值)) +
  geom_col(fill = c("#e8b4a8", "#167d80"), width = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#9aa8a8") +
  geom_text(aes(label = sprintf("%.4f", 比值)), vjust = -0.5, size = 4) +
  labs(title = "体积不一样，计数不能直接比",
       subtitle = "直接比计数得 4.75；按体积折算成密度后，比值降到 1.1875 附近——差的全是体积",
       x = NULL, y = "比值") +
  ylim(0, 5.6)
save_fig("03_offset.png", p3, 7.4, 4.2)

p4 <- ggplot(data.frame(等级 = factor(c("轻", "中", "重"), levels = c("轻", "中", "重")),
                        个数 = c(sum(ord_y == 0), sum(ord_y == 1), sum(ord_y == 2))),
             aes(等级, 个数)) +
  geom_col(fill = "#167d80", width = 0.55) +
  geom_text(aes(label = 个数), vjust = -0.5, size = 4) +
  labs(title = "有序类别：等级之间有顺序，但不是等距的刻度",
       subtitle = "把轻/中/重当 0/1/2 直接回归，等于假设三个等级等距——累积 logit 不假设这一点",
       x = "等级", y = "观测个数") +
  ylim(0, max(table(ord_y)) * 1.3)
save_fig("04_有序类别.png", p4, 7.4, 4.2)

# 图 5 改成两个面板：左边「平均编码得分」，右边「三个等级的预测概率」——两者不是同一个量
fig5_df <- rbind(
  data.frame(面板 = "左：平均编码得分（把等级当 0/1/2）",
             剂量 = c(1:6, ord_x),
             值 = c(ord_curve$当连续变量_平均编码得分, ord_y),
             系列 = c(rep("模型预测的平均得分", 6), rep("原始观测等级", 12))),
  data.frame(面板 = "右：三个等级的预测概率（累积 logit）",
             剂量 = rep(1:6, 3),
             值 = c(ord_curve$轻的类别概率, ord_curve$中的类别概率, ord_curve$重的概率),
             系列 = rep(c("轻", "中", "重"), each = 6))
)
p5 <- ggplot(fig5_df, aes(剂量, 值, colour = 系列)) +
  geom_line(data = subset(fig5_df, 系列 != "原始观测等级"), linewidth = 1.1) +
  geom_point(data = subset(fig5_df, 系列 == "原始观测等级"), size = 2) +
  facet_wrap(~面板, scales = "free_y") +
  labs(title = "两种模型回答的问题不同",
       subtitle = "左：预测的是这套 0/1/2 编码下的平均得分；右：预测的是三个等级各自的概率",
       x = "剂量", y = "值", colour = NULL)
save_fig("05_连续vs有序.png", p5, 7.4, 4.4)

cat("\n完成：5 张图 + 10 份 CSV 已写入\n")
