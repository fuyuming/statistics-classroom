# ============================================================
# 精读 02｜数路径：贝叶斯推断为什么可以先不背公式
# 对应：McElreath 公开课 2023 第 02 讲 Garden of Forking Data；教材第 2 章
# 用法： Rscript 代码/篇02_数路径.R
# 约定：确定性部分（路径计数、网格、解析 Beta、beta-二项预测）三语言逐位一致；
#       随机模拟部分只求形状一致（本篇不出现随机数，图 1–5 全部是确定性计算）。
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
})

OUT <- "文章配图"; RES <- "运行结果"
dir.create(OUT, showWarnings = FALSE); dir.create(RES, showWarnings = FALSE)

# 中文字体：macOS 用 PingFang SC（ragg 设备直接吃系统字体，不用 showtext）
save_fig <- function(name, p, w, h) {
  ragg::agg_png(file.path(OUT, name), width = w, height = h, units = "in", res = 200)
  print(p); dev.off()
  cat("已保存图：", name, "\n")
}
theme_set(theme_minimal(base_size = 12, base_family = "PingFang SC") +
            theme(plot.title = element_text(face = "bold"),
                  panel.grid.minor = element_blank()))

# ------------------------------------------------------------
# 1) 四面地球仪：把"路径"数出来
#    候选水面比例 p = 0、0.25、0.5、0.75、1（分别有 0、1、2、3、4 个水面）
#    观测序列 W L W：每抛一次，路径数乘以"这一面有多少种取法"
# ------------------------------------------------------------
n_faces <- 4                            # 练习 1：把这一行改成 6，就是六面地球仪（三语言同名）
face_water <- seq(0, n_faces)           # 水面数
face_land  <- n_faces - face_water      # 陆面数
ways <- face_water * face_land * face_water   # W、L、W 三抛各自的取法相乘

globe <- tibble(
  p        = face_water / n_faces,        # 比例由"水面数 / 面数"推出，改 n_faces 时不用动这里
  水面数   = face_water,
  陆面数   = face_land,
  W1取法   = face_water,
  L1取法   = face_land,
  W2取法   = face_water,
  路径数   = ways,
  后验     = ways / sum(ways)
)
print(globe)
write.csv(globe, file.path(RES, "01_四面地球仪路径计数.csv"), row.names = FALSE)
cat("路径总数 =", sum(ways), "（= ways 之和）\n")

# ------------------------------------------------------------
# 2) 20 点网格：把"数路径"推广到连续比例
#    观测 9 次里 6 次是水（与 01 篇一致，便于对照）
# ------------------------------------------------------------
W <- 6; L <- 3; N <- W + L
p_grid <- seq(0, 1, length.out = 20)

grid_post <- function(prior) {
  like <- dbinom(W, size = N, prob = p_grid)
  post <- like * prior
  tibble(p = p_grid, prior = prior / sum(prior),
         like = like / sum(like), post = post / sum(post))
}
prior_flat <- rep(1, 20)
prior_exp  <- exp(-5 * abs(p_grid - 0.5))
prior_step <- ifelse(p_grid < 0.5, 0, 1)

g_flat <- grid_post(prior_flat); g_exp <- grid_post(prior_exp); g_step <- grid_post(prior_step)

summ <- function(g, name) {
  lo <- g$p[which(cumsum(g$post) >= 0.055)[1]]
  hi <- g$p[which(cumsum(g$post) >= 0.945)[1]]
  # 统一精度：三语言写出的 CSV 要逐位一致，所以这里把网格点也钉到 6 位小数
  tibble(先验 = name,
         后验均值 = round(sum(g$p * g$post), 6),
         最大点 = round(g$p[which.max(g$post)], 6),
         区间下 = round(lo, 4), 区间上 = round(hi, 4))
}
tab_grid <- bind_rows(
  summ(g_flat, "平坦先验"),
  summ(g_exp,  "exp(-5|p-0.5|)"),
  summ(g_step, "p<0.5 时为 0")
)
print(tab_grid)
write.csv(tab_grid, file.path(RES, "02_三种先验下的后验汇总.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 3) 解析解：Beta(W+1, L+1) 与网格结果的对照
# ------------------------------------------------------------
a <- W + 1; b <- L + 1
an <- tibble(
  指标 = c("后验均值", "后验众数", "后验标准差", "89% 区间下", "89% 区间上"),
  解析值 = c(a / (a + b), (a - 1) / (a + b - 2),
             sqrt(a * b / ((a + b)^2 * (a + b + 1))),
             qbeta(0.055, a, b), qbeta(0.945, a, b))
)
an$解析值 <- round(an$解析值, 6)
cat("解析后验 Beta(", a, ",", b, ")\n", sep = ""); print(an)
write.csv(an, file.path(RES, "03_解析后验对照.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 4) 从后验到预测：下一次取点、以及未来两次取点
#    beta-二项分布是确定性的，所以三语言可以逐位一致
# ------------------------------------------------------------
pred_next <- a / (a + b)                                  # 下一次取到水的概率
k <- 0:2
# P(k 次水 in 2 次) = C(2,k) * B(a+k, b+2-k) / B(a,b)
logB <- function(x, y) lgamma(x) + lgamma(y) - lgamma(x + y)
# 对照：把后验均值当成固定 p 去算（会丢掉后验的宽度）
p_hat <- pred_next
pred2 <- tibble(
  未来两次里的水数 = k,
  `概率（整条后验）` = round(exp(log(choose(2, k)) + logB(a + k, b + 2 - k) - logB(a, b)), 6),
  `概率（固定均值）` = round(choose(2, k) * p_hat^k * (1 - p_hat)^(2 - k), 6)
)
cat("下一次取到水的概率 =", round(pred_next, 6), "\n"); print(pred2)
write.csv(pred2, file.path(RES, "04_预测分布.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 5) 测量误差：判定会看错时，"数路径"怎么重算
#    真样本 3 个水、1 个陆；每次判定有 1/3 的概率看错（3 种可能里 2 种是对的）
#    观测到"水"的路径数 = 真的水：3×2 + 真的陆：1×1
# ------------------------------------------------------------
true_water <- 3; true_land <- 1
# 练习 2：判定取法改成 10 种（判对 9、判错 1）→ judge_correct <- 9; judge_wrong <- 1
judge_correct <- 2; judge_wrong <- 1
mis <- tibble(
  来源       = c("真实是水", "真实是水", "真实是陆", "真实是陆"),
  判定结果   = c("记录为水", "记录为陆", "记录为水", "记录为陆"),
  真样本数   = c(true_water, true_water, true_land, true_land),
  # 每次判定有 3 种取法：判对的 2 种、判错的 1 种
  每次判定的取法 = c(judge_correct, judge_wrong, judge_wrong, judge_correct)
) %>% mutate(路径数 = 真样本数 * 每次判定的取法)
print(mis)
ways_obs_water <- sum(mis$路径数[mis$判定结果 == "记录为水"])
ways_obs_land  <- sum(mis$路径数[mis$判定结果 == "记录为陆"])
cat("观测到水的路径数 = 6 + 1 =", ways_obs_water,
    "；观测到陆的路径数 = 3 + 2 =", ways_obs_land,
    "；合计 =", ways_obs_water + ways_obs_land, "\n")
cat("记录为水的概率 =", round(ways_obs_water / (ways_obs_water + ways_obs_land), 6),
    "（= 7/12）；记录为陆的概率 =", round(ways_obs_land / (ways_obs_water + ways_obs_land), 6), "\n")
write.csv(mis, file.path(RES, "05_误分类路径计数.csv"), row.names = FALSE)

# ============================================================
# 图 1｜分岔路径计数（确定性）
# ============================================================
tree <- globe %>%
  select(p, W1取法, L1取法, W2取法, 路径数) %>%
  pivot_longer(-c(p, 路径数), names_to = "步", values_to = "取法") %>%
  mutate(步 = factor(步, levels = c("W1取法", "L1取法", "W2取法"),
                     labels = c("第 1 抛 W", "第 2 抛 L", "第 3 抛 W")))

p1a <- ggplot(tree, aes(步, factor(p), label = 取法)) +
  geom_tile(fill = "#eef4f3", colour = "white", linewidth = 1.2) +
  geom_text(colour = "#163d48", size = 4.6) +
  labs(title = "每一抛，把这个比例下的取法乘一次",
       subtitle = "四面地球仪：水面 0/1/2/3/4 面，观测序列 W L W",
       x = NULL, y = "候选水面比例 p") +
  theme(panel.grid = element_blank())

p1b <- ggplot(globe, aes(factor(p), 路径数, fill = 路径数 > 0)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0(路径数, "/", sum(ways))), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#c9c9c9", "#167d80"), guide = "none") +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0, 10, 2)) +
  labs(title = "路径数 ÷ 总数 = 后验概率",
       subtitle = paste0("总数 = ", sum(ways), " 条路径；比例 0 和 1 都被排除"),
       x = "候选水面比例 p", y = "路径数")

p1 <- p1a | p1b
save_fig("01-分岔路径计数.png", p1, 8.6, 4.2)

# ============================================================
# 图 2｜20 点网格：先验、似然、后验（确定性）
# ============================================================
curves <- bind_rows(
  g_flat %>% mutate(线 = "先验（平坦）", 值 = prior),
  g_flat %>% mutate(线 = "似然（归一化）", 值 = like),
  g_flat %>% mutate(线 = "后验", 值 = post)
) %>% mutate(线 = factor(线, levels = c("先验（平坦）", "似然（归一化）", "后验")))

p2 <- ggplot(curves, aes(p, 值, colour = 线, linetype = 线)) +
  geom_line(linewidth = 1.15) +
  scale_colour_manual(values = c("#b9b9b9", "#c1462c", "#167d80")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid")) +
  labs(title = "把每个候选比例算一遍，就得到后验",
       subtitle = "观测 9 次里 6 次落在水上；平坦先验下，归一化后的似然与后验重合",
       x = "水面比例 p", y = "归一化权重", colour = NULL, linetype = NULL) +
  theme(legend.position = "top")
save_fig("02-先验似然后验.png", p2, 7.4, 4.2)

# ============================================================
# 图 3｜换先验，后验怎么动（确定性）
# ============================================================
pri_tab <- bind_rows(
  g_flat %>% mutate(情形 = "平坦先验（每个候选一视同仁）"),
  g_exp  %>% mutate(情形 = "exp(-5|p-0.5|)：压在 0.5 附近"),
  g_step %>% mutate(情形 = "p<0.5 时为 0：只认 ≥0.5")
) %>% mutate(情形 = factor(情形, levels = c("平坦先验（每个候选一视同仁）",
                                            "exp(-5|p-0.5|)：压在 0.5 附近",
                                            "p<0.5 时为 0：只认 ≥0.5")))

p3 <- ggplot(pri_tab, aes(p, post, colour = 情形)) +
  geom_line(linewidth = 1.1) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("#b9b9b9", "#167d80", "#c1462c")) +
  labs(title = "同一份数据，换一个先验，后验会被拉多远？",
       subtitle = "观测不变（9 次里 6 次是水）；先验若把某一片区域设成零权重，数据也救不回来",
       x = "水面比例 p", y = "后验概率", colour = NULL) +
  theme(legend.position = "top")
save_fig("03-换先验后的后验.png", p3, 7.6, 4.3)

# ============================================================
# 图 4｜从后验到预测（确定性：beta-二项）
# ============================================================
post_curve <- tibble(p = seq(0, 1, length.out = 200),
                     密度 = dbeta(p, a, b),
                     下一次取到水的概率 = a / (a + b))
p4a <- ggplot(post_curve, aes(p, 密度)) +
  geom_area(fill = "#167d80", alpha = 0.25) +
  geom_line(colour = "#167d80", linewidth = 1.15) +
  geom_vline(xintercept = a / (a + b), linetype = "dashed", colour = "#c1462c") +
  annotate("label", x = a / (a + b), y = max(post_curve$密度) * 0.9,
           label = paste0("下一次取到水的概率\n", round(a / (a + b), 3)),
           family = "PingFang SC", size = 3.4, colour = "#c1462c", hjust = -0.05) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "后验：Beta(7, 4)", x = "水面比例 p", y = "后验密度") +
  theme(plot.title = element_text(size = 12))

pred_long <- pred2 %>%
  pivot_longer(-未来两次里的水数, names_to = "算法", values_to = "概率")
p4b <- ggplot(pred_long, aes(factor(未来两次里的水数), 概率, fill = 算法)) +
  geom_col(position = position_dodge(0.62), width = 0.58) +
  geom_text(aes(label = sprintf("%.3f", 概率)), position = position_dodge(0.62),
            vjust = -0.5, size = 3.6) +
  scale_fill_manual(values = c("#c1462c", "#e8b4a8")) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(title = "预测：未来两次取点里有几次是水",
       subtitle = "深色用整条后验；浅色把后验均值当成固定 p",
       x = "未来两次里的水数", y = "概率", fill = NULL) +
  theme(legend.position = "top")

p4 <- p4a | p4b
save_fig("04-从后验到预测.png", p4, 8.8, 4.0)

# ============================================================
# 图 5｜测量误差：判定会看错时的路径计数（确定性）
# ============================================================
p5 <- ggplot(mis, aes(来源, 路径数, fill = 判定结果)) +
  geom_col(width = 0.55, position = position_dodge(0.62)) +
  geom_text(aes(label = paste0(真样本数, "×", 每次判定的取法, "=", 路径数)),
            position = position_dodge(0.62), vjust = -0.6, size = 3.8) +
  scale_fill_manual(values = c("#167d80", "#c9c9c9")) +
  scale_y_continuous(limits = c(0, 7.6)) +
  labs(title = "判定会看错时，观测结果有多少条路？",
       subtitle = paste0("真样本 3 水 1 陆；判定 3 种取法里判对 2 种、判错 1 种 → ",
                         "观测到水 ", ways_obs_water, " 条路、观测到陆 ", ways_obs_land, " 条路，合计 ",
                         ways_obs_water + ways_obs_land, " 条"),
       x = NULL, y = "路径数", fill = NULL) +
  theme(legend.position = "top")
save_fig("05-误分类路径计数.png", p5, 7.4, 4.0)

cat("\nR 版本:", R.version.string, "\n")
cat("全部完成。\n")
