# ============================================================
# 《Statistical Rethinking》精读 第 01 篇
# 科学先于统计：泥人、DAG 与贝叶斯工作流
# 对应课程：2023 Lecture 01 "Science Before Statistics"（书第 1、2 章）
# 三语言配套之一（R 为主），Python 与 MATLAB 版见同目录
# ============================================================

# 中文字体：ggplot2 装了 ragg 就自动用 ragg 出 png，字体给系统字体名即可。
# macOS 用 "PingFang SC"；Windows 换 "Microsoft YaHei"；Linux 换 "Noto Sans CJK SC"。
library(ggplot2)
library(ragg)
library(patchwork)
FONT <- "PingFang SC"
stopifnot(nrow(systemfonts::match_fonts(FONT)) == 1)

OUT <- "文章配图"; RES <- "运行结果"
dir.create(OUT, showWarnings = FALSE); dir.create(RES, showWarnings = FALSE)

base_theme <- theme_bw(base_size = 13, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey30", size = 11))
theme_set(base_theme)
void_theme <- theme_void(base_family = FONT, base_size = 13) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5),
        plot.subtitle = element_text(colour = "grey30", size = 10, hjust = 0.5))
save_fig <- function(file, plot, w, h) ggsave(file.path(OUT, file), plot,
                                              width = w, height = h, dpi = 240,
                                              bg = "white", device = ragg::agg_png)

# ------------------------------------------------------------
# 第 1 步：贝叶斯更新 —— 地球表面有多少水
# 观测：抛 N 次小球，W 次落在水上。网格近似（grid approximation）
# 89% 等尾区间直接按后验累积和取网格点：不用随机抽样，三语言结果才能逐位一致
# ------------------------------------------------------------
grid_approx <- function(W, N, n_grid = 20, prior_type = "flat") {
  p_grid <- seq(0, 1, length.out = n_grid)
  # 先验：默认平坦（每个候选一视同仁）；练习 2 用 "triangular"（偏向 0.5 的三角形先验）
  prior <- if (prior_type == "flat") rep(1, n_grid) else
    ifelse(p_grid < 0.5, p_grid / 0.5, (1 - p_grid) / 0.5)
  likelihood <- dbinom(W, size = N, prob = p_grid)
  posterior <- likelihood * prior
  posterior <- posterior / sum(posterior)
  data.frame(p = p_grid, prior = prior / sum(prior),
             likelihood = likelihood / sum(likelihood), posterior = posterior)
}
# 练习 2 一行启用（把上一行的默认换成下面这句即可）：
# g_tri <- grid_approx(6, 9, prior_type = "triangular")

cases <- list(c(2, 3), c(6, 9), c(20, 30))

# 用 tidyverse 写法把三组结果拼成一张表（purrr::map_dfr 逐个算完再合并）
library(dplyr); library(purrr); library(tibble)

summarise_one <- function(cs) {
  W <- cs[1]; N <- cs[2]
  g <- grid_approx(W, N)
  cum <- cumsum(g$posterior)
  tibble(
    样本量N = N, 水上次数W = W, 观测比例 = round(W / N, 4),
    网格后验均值 = round(sum(g$p * g$posterior), 6),
    网格最大值点 = round(g$p[which.max(g$posterior)], 6),
    区间下89 = round(g$p[which(cum >= 0.055)[1]], 4),
    区间上89 = round(g$p[which(cum >= 0.945)[1]], 4),
    后验密度峰值 = round(max(g$posterior), 6)
  )
}
res <- map_dfr(cases, summarise_one)
gs <- map(cases, function(cs) {
  g <- grid_approx(cs[1], cs[2])
  g$组 <- sprintf("%d 次里 %d 次是水", cs[2], cs[1])
  g
})
write.csv(res, file.path(RES, "01_网格近似的后验汇总.csv"), row.names = FALSE)
print(as.data.frame(res))

# ------------------------------------------------------------
# 第 2 步：同一组数据，null 模型并不唯一
# 泊松（每个个体期望相同）与负二项（期望本身在变）能产生几乎一样的计数分布
# ------------------------------------------------------------
set.seed(1)
n <- 500
pois <- rpois(n, lambda = 3)
nb <- rnbinom(n, mu = 3, size = 20)          # size 大 = 离散度小 = 贴近泊松
d <- rbind(data.frame(计数 = pois, 过程 = "过程一：泊松（每人期望相同）"),
           data.frame(计数 = nb,  过程 = "过程二：负二项（期望本身在变）"))
d$过程 <- factor(d$过程, levels = c("过程一：泊松（每人期望相同）",
                                    "过程二：负二项（期望本身在变）"))
tab <- table(d$过程, d$计数)
write.csv(as.data.frame.matrix(tab), file.path(RES, "02_两种过程的计数分布.csv"))

p_null <- ggplot(d, aes(x = 计数, fill = 过程)) +
  geom_histogram(position = "identity", alpha = 0.5, binwidth = 1,
                 colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("#c1462c", "#167d80")) +
  scale_x_continuous(breaks = seq(0, 15, 3)) +
  labs(title = "不同生成机制，也可能产生相近的计数分布",
       subtitle = "各模拟 500 个个体，平均事件数都设为 3；相近不等于相同，两者方差并不一样",
       x = "某个体在一段时间内记录到的事件数", y = "个体数", fill = NULL) +
  theme(legend.position = "top")
save_fig("01-两种过程相近的分布.png", p_null, 7.2, 4.4)

# 这两种过程到底有多像：直接比理论概率（确定性计算，三语言可逐位复核）
k <- 0:40
pm_pois <- dpois(k, 3)
pm_nb <- dnbinom(k, size = 20, mu = 3)
tv <- 0.5 * sum(abs(pm_pois - pm_nb))          # 全变差距离
md <- max(abs(pm_pois - pm_nb))                # 最大单点概率差
write.csv(data.frame(指标 = c("全变差距离", "最大单点概率差"), 数值 = round(c(tv, md), 6)),
          file.path(RES, "02_两种过程的近似程度.csv"), row.names = FALSE)
cat(sprintf("泊松 vs 负二项(size=20) 全变差距离 = %.6f，最大单点概率差 = %.6f\n", tv, md))

# ------------------------------------------------------------
# 第 3 步：DAG —— 同一个问题，不同的图对应不同的模型
# 自己用 ggplot 画（不依赖 igraph）：混杂与对撞两种基本结构
# ------------------------------------------------------------
dag_plot <- function(nodes, edges, ttl, note) {
  edges$x    <- nodes$x[match(edges$from, nodes$name)]
  edges$y    <- nodes$y[match(edges$from, nodes$name)]
  edges$xend <- nodes$x[match(edges$to,   nodes$name)]
  edges$yend <- nodes$y[match(edges$to,   nodes$name)]
  dx <- edges$xend - edges$x; dy <- edges$yend - edges$y
  len <- sqrt(dx^2 + dy^2)
  edges$x    <- edges$x    + dx / len * 0.30   # 箭头端点收进节点外沿
  edges$y    <- edges$y    + dy / len * 0.30
  edges$xend <- edges$xend - dx / len * 0.34
  edges$yend <- edges$yend - dy / len * 0.34
  ggplot() +
    geom_segment(data = edges, aes(x, y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
                 linewidth = 1.15, colour = "#163d48") +
    geom_label(data = nodes, aes(x, y, label = name), size = 5.6,
               family = FONT, label.padding = unit(0.45, "lines"),
               linewidth = 0.8, colour = "#163d48", fill = "#e8f3f3",
               label.r = unit(0.5, "lines")) +
    annotate("text", x = 0, y = -1.15, label = note, family = FONT,
             size = 3.5, colour = "grey25") +
    coord_equal(xlim = c(-1.6, 1.6), ylim = c(-1.5, 1.6)) +
    labs(title = ttl) + void_theme
}
nA <- data.frame(name = c("Z", "X", "Y"), x = c(0, -1.1, 1.1), y = c(1.1, -0.2, -0.2))
eA <- data.frame(from = c("Z", "Z"), to = c("X", "Y"))
nB <- data.frame(name = c("X", "Y", "Z"), x = c(-1.1, 1.1, 0), y = c(1.1, 1.1, -0.2))
eB <- data.frame(from = c("X", "Y"), to = c("Z", "Z"))
p_dag <- dag_plot(nA, eA, "甲：混杂（共同原因）", "共同原因 Z 打开了一条混杂路径；估因果效应时需要阻断它") |
  dag_plot(nB, eB, "乙：对撞（共同结果）", "按共同结果 Z 筛选或调整，通常会引入条件关联")
save_fig("02-DAG两种结构.png", p_dag, 8.6, 4.1)

# ------------------------------------------------------------
# 第 4 步：贝叶斯工作流（猫头鹰的五步）
# ------------------------------------------------------------
steps <- data.frame(
  id = 1:5,
  x  = 0,
  y  = c(0, -1.05, -2.10, -3.15, -4.20),
  lab = c("① 说清目标量：我们究竟想知道什么？",
          "② 写科学模型（含因果关系）",
          "③ 由 ①② 推出统计模型",
          "④ 从 ② 生成数据，检验 ③ 能否还原 ①",
          "⑤ 分析真实数据")
)
arrows_df <- data.frame(x = 0, y = steps$y[-5] - 0.34, yend = steps$y[-1] + 0.34)
p_flow <- ggplot() +
  geom_rect(data = steps, aes(xmin = -1.52, xmax = 1.52,
                              ymin = y - 0.34, ymax = y + 0.34),
            fill = "#e8f3f3", colour = "#167d80", linewidth = 1) +
  geom_text(data = steps, aes(0, y, label = lab), family = FONT,
            size = 3.9, lineheight = 1.15, colour = "#163d48") +
  geom_segment(data = arrows_df, aes(x, y, xend = x, yend = yend),
               arrow = arrow(length = unit(0.26, "cm"), type = "closed"),
               colour = "#167d80", linewidth = 0.95) +
  annotate("text", x = 0, y = 1.15,
           label = "统计模型只是这个流程的第三步；跳过 ①② 直接拟合，就是在用没长脑子的泥人",
           family = FONT, size = 3.6, colour = "grey30") +
  coord_cartesian(xlim = c(-1.75, 1.75), ylim = c(-4.9, 1.45)) +
  labs(title = "贝叶斯工作流：先有科学模型，再有统计模型") +
  theme_void(base_family = FONT, base_size = 13) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
save_fig("03-贝叶斯工作流五步.png", p_flow, 7.4, 5.0)

# ------------------------------------------------------------
# 第 5 步：网格近似的三张曲线（先验 / 似然 / 后验）
# ------------------------------------------------------------
g2 <- gs[[2]]
gl <- rbind(transform(g2, 成分 = "先验（还没看数据）", 密度 = g2$prior),
            transform(g2, 成分 = "似然（数据说什么）", 密度 = g2$likelihood),
            transform(g2, 成分 = "后验（看完之后）", 密度 = g2$posterior))
gl$成分 <- factor(gl$成分, levels = c("先验（还没看数据）", "似然（数据说什么）", "后验（看完之后）"))
p_curve <- ggplot(gl, aes(p, 密度, colour = 成分, linetype = 成分)) +
  geom_line(linewidth = 1.15) + geom_point(size = 1.7) +
  # 平坦先验下后验与似然形状完全重合，把似然的虚线单独再画一层压在上面，否则看不见
  geom_line(data = subset(gl, 成分 == "似然（数据说什么）"),
            aes(p, 密度), inherit.aes = FALSE, colour = "#c1462c",
            linetype = "dashed", linewidth = 1.1) +
  scale_colour_manual(values = c("#8a8a8a", "#c1462c", "#167d80")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid")) +
  labs(title = "把 20 个候选比例各自算一遍，就得到后验",
       subtitle = "观测 9 次里 6 次落在水上；先验平坦时，后验与似然形状完全重合",
       x = "水面比例 p", y = "归一化权重", colour = NULL, linetype = NULL) +
  theme(legend.position = "top")
save_fig("04-先验似然后验.png", p_curve, 7.2, 4.2)

# ------------------------------------------------------------
# 第 6 步：数据越多，后验越窄
# ------------------------------------------------------------
wide <- do.call(rbind, gs)
wide$组 <- factor(wide$组, levels = c("3 次里 2 次是水", "9 次里 6 次是水", "30 次里 20 次是水"))
p_narrow <- ggplot(wide, aes(p, posterior, colour = 组)) +
  geom_line(linewidth = 1.15) +
  scale_colour_manual(values = c("#b9b9b9", "#198c8b", "#163d48")) +
  labs(title = "观测比例都是 2/3，样本越多，后验越集中",
       subtitle = "每一个点都是「这一步之后你还能接受的候选比例」",
       x = "水面比例 p", y = "网格点后验概率", colour = NULL) +
  theme(legend.position = "top")
save_fig("05-样本量越大后验越窄.png", p_narrow, 7.2, 4.2)

cat("\n出图完成：\n"); print(list.files(OUT))
cat("\nR 版本:", R.version.string, "\n")
