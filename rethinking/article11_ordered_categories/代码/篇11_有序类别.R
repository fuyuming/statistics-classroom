# ============================================================
# 精读 11｜有序类别：切点自由估计，以及"合并成 0/1"会丢什么
# 对应：McElreath 2023 第 11 讲 Ordered Categories；教材第 11 章
# 用法： Rscript 代码/篇11_有序类别.R
#
# 与第 09 篇的关系：09 篇为了讲清累积概率，把切点固定；这一篇把切点放开，
# 用"有序约束"参数化（c2 = c1 + exp(v)）让两个切点在网格上自动保持次序。
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

# ---------- 1) 数据：30 个样本，剂量 1–10，三档等级 ----------
x  <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10)
y  <- c(0,0,1, 0,0,1, 0,1,1, 0,1,1, 1,1,1, 1,1,2, 1,1,2, 1,2,2, 2,2,2, 2,2,2)
d  <- data.frame(剂量 = x, 等级 = y)
write_csv_fixed(d, file.path(RES, "01_数据.csv"))

# ---------- 2) 累积 logit（切点自由）：logit(F_k) = c_k − b·x，k = 0,1 ----------
# 参数化：c1 = u，c2 = u + exp(v)，保证 c1 < c2
loglik <- function(b, u, v) {
  c1 <- u; c2 <- u + exp(v)
  F1 <- 1 / (1 + exp(-(c1 - b * x)))
  F2 <- 1 / (1 + exp(-(c2 - b * x)))
  P <- cbind(F1, F2 - F1, 1 - F2)                       # 每行一个观测的三档概率
  sum(log(pmax(P[cbind(seq_along(y), y + 1)], 1e-300)))
}
b_seq <- seq(-2, 3, length.out = 70)
u_seq <- seq(-4, 4, length.out = 70)
v_seq <- seq(log(0.05), log(5), length.out = 70)
g <- expand.grid(b = b_seq, u = u_seq, v = v_seq)
lp <- numeric(nrow(g))
for (i in seq_len(nrow(g))) lp[i] <- loglik(g$b[i], g$u[i], g$v[i])
g$post <- exp(lp - max(lp)); g$post <- g$post / sum(g$post)
g$c1 <- g$u; g$c2 <- g$u + exp(g$v)

post_tab <- data.frame(
  参数 = c("斜率 b", "切点 c1（轻/中）", "切点 c2（中/重）"),
  后验均值 = r4(c(weighted.mean(g$b, g$post), weighted.mean(g$c1, g$post), weighted.mean(g$c2, g$post))),
  下界89 = r4(c(wq(g$b, g$post, 0.055), wq(g$c1, g$post, 0.055), wq(g$c2, g$post, 0.055))),
  上界89 = r4(c(wq(g$b, g$post, 0.945), wq(g$c1, g$post, 0.945), wq(g$c2, g$post, 0.945)))
)
write_csv_fixed(post_tab, file.path(RES, "02_后验.csv"))
print(post_tab)
cat("网格：70×70×70 =", nrow(g), "个格点（三参数网格，切点次序由参数化保证）\n")

# ---------- 3) 累积概率曲线与各等级概率 ----------
dose_fine <- 1:10
curve_rows <- list()
for (dd in dose_fine) {
  F1 <- 1 / (1 + exp(-(g$c1 - g$b * dd)))
  F2 <- 1 / (1 + exp(-(g$c2 - g$b * dd)))
  curve_rows[[length(curve_rows) + 1]] <- data.frame(
    剂量 = dd,
    累积_轻 = r4(weighted.mean(F1, g$post)),
    累积_中 = r4(weighted.mean(F2, g$post)),
    概率_轻 = r4(weighted.mean(F1, g$post)),
    概率_中 = r4(weighted.mean(F2 - F1, g$post)),
    概率_重 = r4(weighted.mean(1 - F2, g$post))
  )
}
curve_df <- do.call(rbind, curve_rows)
# 每一行的三个概率之和应当为 1（数值检查）
curve_df$概率合计 <- r4(curve_df$概率_轻 + curve_df$概率_中 + curve_df$概率_重)
write_csv_fixed(curve_df, file.path(RES, "03_累积概率与等级概率.csv"))
print(curve_df)

# ---------- 4) 后验预测检查：预测的等级比例 vs 观测 ----------
obs_prop <- as.numeric(table(factor(y, levels = 0:2))) / length(y)
pred_rows <- sapply(dose_fine, function(dd) {
  F1 <- 1 / (1 + exp(-(g$c1 - g$b * dd)))
  F2 <- 1 / (1 + exp(-(g$c2 - g$b * dd)))
  c(weighted.mean(F1, g$post), weighted.mean(F2 - F1, g$post), weighted.mean(1 - F2, g$post))
})
pred_prop <- rowMeans(pred_rows)
ppc <- data.frame(
  等级 = c("轻", "中", "重"),
  观测比例 = r4(obs_prop),
  预测比例 = r4(pred_prop),
  差 = r4(pred_prop - obs_prop)
)
write_csv_fixed(ppc, file.path(RES, "04_后验预测检查.csv"))
print(ppc)

# ---------- 5) 把等级合并成 0/1 会丢什么 ----------
# 二分类：等级 ≥ 2 记 1（重 vs 其余），用同一个 logistic 模型
y_bin <- as.integer(y >= 2)
loglik_bin <- function(a, b) {
  p <- 1 / (1 + exp(-(a + b * x))); p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  sum(y_bin * log(p) + (1 - y_bin) * log(1 - p))
}
ga <- expand.grid(a = seq(-8, 3, length.out = 300), b = seq(-1.5, 2.5, length.out = 300))
lpa <- numeric(nrow(ga))
for (i in seq_len(nrow(ga))) lpa[i] <- loglik_bin(ga$a[i], ga$b[i])
ga$post <- exp(lpa - max(lpa)); ga$post <- ga$post / sum(ga$post)
se_bin <- sqrt(weighted.mean((ga$b - weighted.mean(ga$b, ga$post))^2, ga$post))
se_ord <- sqrt(weighted.mean((g$b - weighted.mean(g$b, g$post))^2, g$post))
cmp <- data.frame(
  写法 = c("有序模型（三档全用）", "合并成二分类（重 vs 其余）"),
  斜率后验均值 = r4(c(weighted.mean(g$b, g$post), weighted.mean(ga$b, ga$post))),
  斜率后验标准差 = r4(c(se_ord, se_bin)),
  说明 = c("斜率是「累积几率的共同位移」", "斜率是「重 vs 其余」的几率变化——不是同一个量")
)
write_csv_fixed(cmp, file.path(RES, "05_有序vs二分类.csv"))
print(cmp)
# 信息损失的具体形态：二分类模型说不出「轻 vs 中」的区别
loss_rows <- list()
for (dd in c(3, 5, 7)) {
  F1 <- 1 / (1 + exp(-(g$c1 - g$b * dd)))
  F2 <- 1 / (1 + exp(-(g$c2 - g$b * dd)))
  p_bin <- 1 / (1 + exp(-(weighted.mean(ga$a, ga$post) + weighted.mean(ga$b, ga$post) * dd)))
  loss_rows[[length(loss_rows) + 1]] <- data.frame(
    剂量 = dd,
    有序_轻 = r4(weighted.mean(F1, g$post)),
    有序_中 = r4(weighted.mean(F2 - F1, g$post)),
    有序_重 = r4(weighted.mean(1 - F2, g$post)),
    二分类_重 = r4(p_bin),
    二分类_其余 = r4(1 - p_bin)
  )
}
loss_df <- do.call(rbind, loss_rows)
write_csv_fixed(loss_df, file.path(RES, "06_信息损失的形态.csv"))
print(loss_df)

# ---------- 图 ----------
stack_df <- curve_df %>%
  select(剂量, 轻 = 概率_轻, 中 = 概率_中, 重 = 概率_重) %>%
  pivot_longer(-剂量, names_to = "等级", values_to = "概率")
stack_df$等级 <- factor(stack_df$等级, levels = c("重", "中", "轻"))
p1 <- ggplot(stack_df, aes(factor(剂量), 概率, fill = 等级)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(重 = "#c1462c", 中 = "#e8b4a8", 轻 = "#167d80")) +
  labs(title = "三档等级的概率怎么随剂量变化",
       subtitle = "每一根柱子的三段加起来都是 1；切点是自由估出来的（本例未固定）",
       x = "剂量", y = "概率", fill = NULL)
save_fig("01_等级概率.png", p1, 7.4, 4.4)

dens <- rbind(
  data.frame(参数 = "斜率 b", 值 = g$b, 权重 = g$post),
  data.frame(参数 = "切点 c1", 值 = g$c1, 权重 = g$post),
  data.frame(参数 = "切点 c2", 值 = g$c2, 权重 = g$post)
)
p2 <- dens %>%
  group_by(参数) %>%
  summarise(均值 = weighted.mean(值, 权重),
            下 = wq(值, 权重, 0.055), 上 = wq(值, 权重, 0.945), .groups = "drop") %>%
  ggplot(aes(参数, 均值)) +
  geom_point(size = 3, colour = "#167d80") +
  geom_errorbar(aes(ymin = 下, ymax = 上), width = 0.15, colour = "#163d48") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#9aa8a8") +
  geom_text(aes(label = sprintf("%.3f", 均值)), vjust = -0.8, size = 3.8) +
  labs(title = "三个参数一起估：斜率与两个切点",
       subtitle = "c2 明显大于 c1 —— 次序由参数化 c2 = c1 + exp(v) 保证",
       x = NULL, y = "后验均值（89% 区间）")
save_fig("02_参数.png", p2, 7.4, 4.2)

curve_long <- curve_df %>%
  select(剂量, 轻的累积 = 累积_轻, 中的累积 = 累积_中) %>%
  pivot_longer(-剂量, names_to = "曲线", values_to = "概率")
p3 <- ggplot(curve_long, aes(剂量, 概率, colour = 曲线)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2) +
  labs(title = "累积概率：这才是模型直接预测的东西",
       subtitle = "两条累积曲线之差就是中间那一档的概率；次序由切点保证",
       x = "剂量", y = "累积概率", colour = NULL)
save_fig("03_累积概率.png", p3, 7.4, 4.4)

ppc_long <- ppc %>%
  select(等级, 观测比例, 预测比例) %>%
  pivot_longer(-等级, names_to = "来源", values_to = "比例")
p4 <- ggplot(ppc_long, aes(等级, 比例, fill = 来源)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", 比例)), position = position_dodge(width = 0.7), vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = c(观测比例 = "#167d80", 预测比例 = "#e8b4a8")) +
  labs(title = "后验预测检查：预测的等级比例 vs 观测",
       subtitle = "差得不远说明模型抓住了大方向；差异大的档位要回去看数据",
       x = NULL, y = "比例", fill = NULL)
save_fig("04_预测检查.png", p4, 7.4, 4.2)

# 图 5：两种模型能回答的问题根本不同——二分类说不出「轻 vs 中」
loss_long <- rbind(
  data.frame(剂量 = loss_df$剂量, 模型 = "有序模型", 类别 = "轻", 概率 = loss_df$有序_轻),
  data.frame(剂量 = loss_df$剂量, 模型 = "有序模型", 类别 = "中", 概率 = loss_df$有序_中),
  data.frame(剂量 = loss_df$剂量, 模型 = "有序模型", 类别 = "重", 概率 = loss_df$有序_重),
  data.frame(剂量 = loss_df$剂量, 模型 = "合并成二分类", 类别 = "重", 概率 = loss_df$二分类_重),
  data.frame(剂量 = loss_df$剂量, 模型 = "合并成二分类", 类别 = "其余（轻+中，无法再分）", 概率 = loss_df$二分类_其余)
)
loss_long$类别 <- factor(loss_long$类别, levels = c("重", "中", "轻", "其余（轻+中，无法再分）"))
p5 <- ggplot(loss_long, aes(factor(剂量), 概率, fill = 类别)) +
  geom_col(width = 0.7) +
  facet_wrap(~模型) +
  scale_fill_manual(values = c(重 = "#c1462c", 中 = "#e8b4a8", 轻 = "#167d80",
                               "其余（轻+中，无法再分）" = "#9aa8a8")) +
  labs(title = "合并成 0/1 之后，有一类问题就再也答不了",
       subtitle = "左：三档概率都能算；右：只剩「重 vs 其余」，轻与中的区别被永久合并",
       x = "剂量", y = "概率", fill = NULL)
save_fig("05_有序vs二分类.png", p5, 7.4, 4.2)

cat("\n完成：5 张图 + 5 份 CSV 已写入\n")
