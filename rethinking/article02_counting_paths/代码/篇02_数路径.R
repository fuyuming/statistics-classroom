# ============================================================
# 精读 02｜数路径：贝叶斯推断为什么可以先不背公式
# 对应：McElreath 公开课 2023 第 02 讲 Garden of Forking Data；教材第 2 章
# 用法： Rscript 代码/篇02_数路径.R
# 约定：本篇全部为确定性教学计算，无随机模拟；三语言按数值容差核对。
# 新增 Beta 密度在六位小数舍入边界允许末位差 0.000001，详见 VALIDATION.md。
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
})

# 命令行与 RStudio Source 自动定位脚本；逐行运行请先打开本课 Rproj。
source_files <- unlist(lapply(sys.frames(), function(x) x$ofile))
cli <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
script_file <- if (length(source_files)) tail(source_files, 1) else cli
root <- if (length(script_file)) dirname(dirname(normalizePath(script_file))) else getwd()
OUT <- file.path(root, "文章配图_v3")
RES <- file.path(root, "运行结果")
dir.create(OUT, showWarnings = FALSE); dir.create(RES, showWarnings = FALSE)

# 中文字体：macOS 用 PingFang SC（ragg 设备直接吃系统字体，不用 showtext）
save_fig <- function(name, p, w, h) {
  ragg::agg_png(file.path(OUT, name), width = w, height = h, units = "in", res = 200)
  print(p); dev.off()
  cat("已保存图：", name, "\n")
}
font_family <- if (Sys.info()[["sysname"]] == "Darwin") "PingFang SC" else if (.Platform$OS.type == "windows") "Microsoft YaHei" else "Noto Sans CJK SC"
theme_set(theme_minimal(base_size = 12, base_family = font_family) +
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
write.csv(mutate(globe, across(everything(), ~round(.x, 6))),
          file.path(RES, "01_四面地球仪路径计数.csv"), row.names = FALSE)
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
  # 统一精度：三语言采用相同输出精度，这里把网格点保留到 6 位小数
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
# 平坦先验是 Beta(1,1)。一次水乘一个 p，一次陆乘一个 (1-p)。
# Beta(a,b) 的幂次是 a-1、b-1，因此 6 水 3 陆对应 Beta(7,4)。
a <- W + 1; b <- L + 1
cat("为什么是 Beta？平坦先验 Beta(1,1)，加上水和陆的次数，参数成为", a, b, "。\n")
cat("dbeta 返回密度高度；qbeta 返回累计面积达到指定比例时的横轴位置。\n")
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
#    这里采用确定性预测计算，三语言的结果可以按数值容差核对
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
cat("记录为水有", ways_obs_water, "条路；记录为陆有", ways_obs_land,
    "条路；总数", ways_obs_water + ways_obs_land, "\n")
cat("记录为水的概率 =", ways_obs_water / (ways_obs_water + ways_obs_land), "\n")
write.csv(mis, file.path(RES, "05_误分类路径计数.csv"), row.names = FALSE)

# 6) 同一组 W-L-W：把测量误差一路算进后验
# 各次取点与误判独立，误判对称，候选先验等权。
error_rate <- judge_wrong / (judge_correct + judge_wrong)
obs_water_ways <- face_water * judge_correct + face_land * judge_wrong
obs_land_ways <- face_water * judge_wrong + face_land * judge_correct
q_record_water <- obs_water_ways / (n_faces * (judge_correct + judge_wrong))
# q 是记录为水的概率；q*(1-q)*q 是指定序列的似然。
mis_likelihood <- q_record_water * (1 - q_record_water) * q_record_water
mis_posterior <- mis_likelihood / sum(mis_likelihood)
mis_paths <- obs_water_ways * obs_land_ways * obs_water_ways
mis_post <- tibble(p = globe$p, 误判率 = error_rate,
                   记录为水的概率 = q_record_water, 相容路径数 = mis_paths,
                   无误判后验 = globe$后验, 含误判后验 = mis_posterior)
write.csv(mutate(mis_post, across(everything(), ~round(.x, 6))),
          file.path(RES, "06_误判后验对照.csv"), row.names = FALSE)
print(mis_post)
cat("端点重新获得支持，是因为含误判模型的似然不再为零；此处先验一直等权。\n")

# 7) Beta 读图：先画出三个累计阶段，不用分布名字代替理解
# 总面积均为 1，图上某一点的高度不是该点的概率。
beta_cases <- tibble(阶段 = c("平坦先验", "1水1陆", "6水3陆"),
                     a = c(1, 2, a), b = c(1, 2, b))
beta_read <- bind_rows(lapply(seq_len(nrow(beta_cases)), function(i) {
  p <- seq(0, 1, length.out = 201)
  tibble(阶段 = beta_cases$阶段[i], a = beta_cases$a[i], b = beta_cases$b[i],
         p = p, 密度 = dbeta(p, beta_cases$a[i], beta_cases$b[i]))
}))
write.csv(mutate(beta_read, across(where(is.numeric), ~round(.x, 6))),
          file.path(RES, "07_Beta读图.csv"), row.names = FALSE)

# 8) 两候选手算：这是独立教学设定，不是 6 水 3 陆的后验
candidate_p <- c(0.25, 0.75)
candidate_weight <- c(0.5, 0.5)
toy_mean <- sum(candidate_p * candidate_weight)
toy_prediction <- tibble(
  算法 = c("逐个候选预测再平均", "先平均再预测"),
  下一次是水 = c(sum(candidate_weight * candidate_p), toy_mean),
  两次都是水 = c(sum(candidate_weight * candidate_p^2), toy_mean^2)
)
write.csv(toy_prediction, file.path(RES, "08_两候选预测.csv"), row.names = FALSE)
print(toy_prediction)
cat("预测两次时，0.3125 与 0.25 不同：先平方再平均，不等于先平均再平方。\n")

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
       subtitle = paste0(n_faces, " 面地球仪；各面等可能；观测 W L W"),
       x = NULL, y = "候选水面比例 p") +
  theme(panel.grid = element_blank())

p1b <- ggplot(globe, aes(factor(p), 路径数, fill = 路径数 > 0)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0(路径数, "/", sum(ways))), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("#c9c9c9", "#167d80"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "路径数 ÷ 总数 = 后验概率",
       subtitle = paste0("相容路径合计 ", sum(ways), " 条\n先验等权、独立取点、无误判"),
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
save_fig("03-换先验后的后验.png", p3, 9, 4.8)

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
           family = font_family, size = 3.4, colour = "#c1462c", hjust = -0.05) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "6 水 3 陆后的水面比例分布", x = "水面比例 p", y = "后验密度") +
  theme(plot.title = element_text(size = 12))

pred_long <- pred2 %>%
  pivot_longer(-未来两次里的水数, names_to = "算法", values_to = "概率") %>%
  # 钉死因子顺序，否则中文列名会按字典序重排，图例颜色就跟说明对不上
  mutate(算法 = factor(算法, levels = c("概率（整条后验）", "概率（固定均值）")))
p4b <- ggplot(pred_long, aes(factor(未来两次里的水数), 概率, fill = 算法)) +
  geom_col(position = position_dodge(0.8), width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", 概率)), position = position_dodge(0.8),
            vjust = -0.45, size = 3.1) +
  scale_fill_manual(values = c("#c1462c", "#e8b4a8")) +
  scale_y_continuous(limits = c(0, 0.5)) +
  labs(title = "预测：未来两次取点里有几次是水",
       subtitle = "深色：用整条后验；浅色：把后验均值当成固定的 p",
       x = "未来两次里的水数", y = "概率", fill = NULL) +
  theme(legend.position = "top")

p4 <- p4a | p4b
save_fig("04-从后验到预测.png", p4, 10.8, 4.8)

# ============================================================
# 图 5｜测量误差：判定会看错时的路径计数（确定性）
# ============================================================
p5 <- ggplot(mis, aes(来源, 路径数, fill = 判定结果)) +
  geom_col(width = 0.55, position = position_dodge(0.62)) +
  geom_text(aes(label = paste0(真样本数, "×", 每次判定的取法, "=", 路径数)),
            position = position_dodge(0.62), vjust = -0.6, size = 3.8) +
  scale_fill_manual(values = c("#167d80", "#c9c9c9")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "判定会看错时，观测结果有多少条路？",
       subtitle = paste0("3 水 1 陆；判对分支 ", judge_correct, "，判错分支 ", judge_wrong, " → ",
                         "观测到水 ", ways_obs_water, " 条路、观测到陆 ", ways_obs_land, " 条路，合计 ",
                         ways_obs_water + ways_obs_land, " 条"),
       x = NULL, y = "路径数", fill = NULL) +
  theme(legend.position = "top")
save_fig("05-误分类路径计数.png", p5, 8.6, 4.5)

cat("\nR 版本:", R.version.string, "\n")
cat("全部完成。\n")

# 图3（文件编号08）：Beta 来历与面积；使用相同坐标尺度逐层阅读。
stage_labels <- c("起点：还没有看数据\nBeta(1, 1)",
                  "累计看到 1 水 1 陆\nBeta(2, 2)",
                  paste0("累计看到 ", W, " 水 ", L, " 陆\nBeta(", a, ", ", b, ")"))
beta_read$阶段图 <- factor(beta_read$阶段, levels = beta_cases$阶段, labels = stage_labels)
ci_bounds <- qbeta(c(0.055, 0.945), a, b)
shade_p <- seq(ci_bounds[1], ci_bounds[2], length.out = 201)
shade <- tibble(p = shade_p, 密度 = dbeta(shade_p, a, b),
                阶段图 = factor(stage_labels[3], levels = stage_labels))
p_beta <- ggplot(beta_read, aes(p, 密度)) +
  geom_area(data = shade, fill = "#167d80", alpha = 0.22) +
  geom_line(colour = "#167d80", linewidth = 1.1) +
  facet_wrap(~阶段图, ncol = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 3.1), breaks = 0:3) +
  labs(title = "数完路径，为什么得到一条 Beta 曲线？",
       subtitle = "同样宽的范围，面积越大，后验概率越大；每幅图总面积均为 1",
       x = "候选水面比例 p", y = "密度：概率在附近有多集中",
       caption = "最下图阴影：89% 后验概率；左右各留下 5.5%\n平坦先验、独立均匀取点、水陆无误判") +
  theme(strip.text = element_text(size = 12, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 10),
        panel.spacing = grid::unit(0.7, "lines"))
save_fig("08-Beta从哪里来.png", p_beta, 7.5, 9)
