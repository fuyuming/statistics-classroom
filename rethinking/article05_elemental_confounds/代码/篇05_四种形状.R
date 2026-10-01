# ============================================================
# 精读 05｜混杂的四种基本形状：叉、管、对撞、后代
# 对应：McElreath 2023 第 05 讲 Elemental Confounds；教材第 5、6 章
# 用法： Rscript 代码/篇05_四种形状.R
#
# 约定：
#   · 主线案例＝我们自己的场景（批次/接种/菌群密度/筛选）
#   · 对照案例＝教材 Howell1 与 WaffleDivorce（华夫饼店与离婚率那类荒谬相关）
#   · 全部确定性：显式最小二乘（解正规方程）、显式皮尔逊相关 → 三语言逐位一致
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(rethinking)
})

.args <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .args, value = TRUE)
root <- if (length(.f) > 0) dirname(dirname(normalizePath(sub("^--file=", "", .f[1])))) else getwd()
OUT  <- file.path(root, "文章配图")
RES  <- file.path(root, "运行结果")
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
  list(coef = as.numeric(b), sigma_hat = sqrt(sum((y - X %*% b)^2) / (length(y) - ncol(X))))
}
pearson <- function(x, y) {
  x <- x - mean(x); y <- y - mean(y)
  sum(x * y) / sqrt(sum(x^2) * sum(y^2))
}
r4 <- function(x) round(as.numeric(x), 4)

# ---------- 0) 教材对照数据 ----------
data(Howell1); data(WaffleDivorce)
wd <- WaffleDivorce %>%
  transmute(州 = as.character(Location), 离婚率 = Divorce,
            华夫饼店 = WaffleHouses, 人口 = Population, 南方 = South)
write_csv_fixed(wd, file.path(RES, "00_教材数据_华夫饼与离婚率.csv"))
# 皮尔逊相关：全体，以及按"南方/非南方"分层
r_all <- pearson(wd$华夫饼店, wd$离婚率)
r_south <- pearson(wd$华夫饼店[wd$南方 == 1], wd$离婚率[wd$南方 == 1])
r_north <- pearson(wd$华夫饼店[wd$南方 == 0], wd$离婚率[wd$南方 == 0])
wd_cor <- tibble(范围 = c("全体 50 个州", "南方各州", "非南方各州"),
                 华夫饼店与离婚率的相关 = r4(c(r_all, r_south, r_north)),
                 州数 = c(nrow(wd), sum(wd$南方 == 1), sum(wd$南方 == 0)))
write_csv_fixed(wd_cor, file.path(RES, "01_教材对照_分层相关.csv"))
print(wd_cor)

# ---------- 2) 我们的案例（一）：叉（混杂）----------
# 批次 → 接种与长势；接种的安排与批次相关（不完全随机）
treat_c <- c(1, 1, 1, 1, 0, 1, 0, 0, 0, 0)          # 10 个实验单元：处理（1）与对照（0）
batch_c <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)          # 批次 1 偏处理、批次 2 偏对照
resid_c <- c(0.6, -0.8, 1.1, -0.5, 0.3, 0.9, -1.2, 0.4, -0.6, 0.8)
growth_c <- 10 + 3.0 * treat_c + 2.5 * batch_c + resid_c
fork_df <- tibble(单元 = 1:10, 处理 = treat_c, 批次 = ifelse(batch_c == 1, "批次2", "批次1"), 长势 = growth_c)
write_csv_fixed(fork_df, file.path(RES, "02_我们的案例_叉.csv"))

f_naive <- ols_ls(cbind(1, treat_c), growth_c)
f_adj   <- ols_ls(cbind(1, treat_c, batch_c), growth_c)
fork_est <- tibble(
  写法 = c("不调整（长势 ~ 处理）", "把批次放进模型（长势 ~ 处理 + 批次）"),
  处理效应估计 = r4(c(f_naive$coef[2], f_adj$coef[2])),
  真值 = c(3.0, 3.0)
)
write_csv_fixed(fork_est, file.path(RES, "03_叉_调整前后.csv"))
print(fork_est)
# 处理与批次的相关（说明"不完全随机"）
r_tb <- pearson(treat_c, batch_c) %>% r4()
cat("处理与批次的相关 =", r_tb, "\n")

# ---------- 3) 我们的案例（二）：管（中介）----------
# 接种 → 菌群密度 → 长势；菌群密度是中介（也是处理后变量）
resid_d <- c(0.4, -0.6, 0.9, -0.3, 0.5, -0.8, 0.7, -0.2, 0.6, -0.9)
resid_g <- c(0.7, 0.5, -0.8, 0.9, -0.4, 0.3, -0.6, 0.8, -0.5, 0.2)
density_d <- 5 + 2.0 * treat_c + resid_d
growth_d  <- 8 + 1.5 * density_d + resid_g
pipe_df <- tibble(单元 = 1:10, 处理 = treat_c, 菌群密度 = round(density_d, 6), 长势 = round(growth_d, 6))
write_csv_fixed(pipe_df, file.path(RES, "04_我们的案例_管.csv"))
p_total <- ols_ls(cbind(1, treat_c), growth_d)                 # 总效应
p_direct <- ols_ls(cbind(1, treat_c, density_d), growth_d)      # 控制中介后的直接效应
pipe_est <- tibble(
  写法 = c("不控制中介（长势 ~ 处理）", "控制中介（长势 ~ 处理 + 菌群密度）"),
  处理效应估计 = r4(c(p_total$coef[2], p_direct$coef[2])),
  机制真值 = c(3.0, 0.0),
  说明 = c("估计的是总效应", "在适当识别假设下估计直接效应；本例机制里直接通路不存在")
)
write_csv_fixed(pipe_est, file.path(RES, "05_管_总效应与直接效应.csv"))
print(pipe_est)

# ---------- 4) 我们的案例（三）：对撞（选择）----------
# 处理与基线值都影响"是否被选中测序"；只看被选中的样本 → 两者出现负相关
base_b <- c(-1.8, -1.2, -0.7, -0.3, 0.1, 0.4, 0.9, 1.3, 1.7, 2.0)
treat_b <- c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
resid_s <- c(0.5, -0.4, 0.2, 0.8, -0.9, 0.3, -0.2, 0.6, -0.5, 0.1)
score_s <- 0.6 * treat_b + 1.2 * base_b + resid_s           # "被选中"的倾向
selected <- as.integer(score_s > median(score_s))           # 取一半样本"被选中"
coll_df <- tibble(单元 = 1:10, 处理 = treat_b, 基线值 = base_b,
                  选中 = selected, 倾向分 = round(score_s, 6))
write_csv_fixed(coll_df, file.path(RES, "06_我们的案例_对撞.csv"))
r_coll_all <- pearson(treat_b, base_b) %>% r4()
r_coll_sel <- pearson(treat_b[selected == 1], base_b[selected == 1]) %>% r4()
coll_tab <- tibble(
  范围 = c("全部 10 个单元", "只看被选中的 5 个单元"),
  处理与基线值的相关 = c(r_coll_all, r_coll_sel),
  单元数 = c(10L, sum(selected))
)
write_csv_fixed(coll_tab, file.path(RES, "07_对撞_条件化前后.csv"))
print(coll_tab)
# 入选者的处理/对照构成，以及"收紧筛选"（只取前 3 名）时的相关
sel_treat <- sum(treat_b[selected == 1] == 1)
sel_ctrl  <- sum(treat_b[selected == 1] == 0)
top3_idx  <- order(score_s, decreasing = TRUE)[1:3]
r_coll_top3 <- pearson(treat_b[top3_idx], base_b[top3_idx]) %>% r4()
sel_comp <- tibble(
  范围 = c("全部 10 个单元", "只看被选中的 5 个单元（处理 4、对照 1）", "只取分数最高的 3 个"),
  处理与基线值的相关 = c(r_coll_all, r_coll_sel, r_coll_top3),
  单元数 = c(10L, sum(selected), 3L)
)
write_csv_fixed(sel_comp, file.path(RES, "10_入选构成与前3名.csv"))
print(sel_comp)

# ---------- 5) 我们的案例（四）：后代——没有控制那个节点，却控制了它的后代 ----------
# 结构：处理 → 入选 ← 基线；入选 → 是否有报告。
# 只分析"有报告"的单元，等于通过后代间接按"入选"做了筛选。
rep_noise <- c(0.4, -0.6, 0.7, -0.2, 0.5, -0.8, 0.3, -0.1, 0.6, -0.5)
has_report <- as.integer(selected == 1 & rep_noise > -0.55)
desc_df <- tibble(单元 = 1:10, 处理 = treat_b, 基线值 = base_b,
                  入选 = selected, 报告扰动 = round(rep_noise, 6), 有报告 = has_report)
write_csv_fixed(desc_df, file.path(RES, "08_我们的案例_后代.csv"))
r_coll_rep <- pearson(treat_b[has_report == 1], base_b[has_report == 1]) %>% r4()
desc_est <- tibble(
  范围 = c("全部 10 个单元", "只看被选中的 5 个单元", "只看有报告的单元"),
  处理与基线值的相关 = c(r_coll_all, r_coll_sel, r_coll_rep),
  单元数 = c(10L, sum(selected), sum(has_report))
)
write_csv_fixed(desc_est, file.path(RES, "09_后代_间接筛选.csv"))
print(desc_est)

# ---------- 图 ----------
# 图 1：四种结构的小图（坐标直接写进边表，不靠 match 推，避免坐标成 NA）
mk_panel <- function(panel_name, nodes, edges) {
  n <- tibble(panel = panel_name, name = nodes$name, x = nodes$x, y = nodes$y)
  e <- as_tibble(edges)
  e$panel <- panel_name
  list(n = n, e = e)
}
P1 <- mk_panel("叉（混杂）",
  tibble(name = c("批次 Z", "处理 X", "长势 Y"), x = c(0.50, 0.22, 0.78), y = c(1.05, 0.0, 0.0)),
  tibble(from = c("批次 Z", "批次 Z", "处理 X"), to = c("处理 X", "长势 Y", "长势 Y"),
         x = c(0.50, 0.50, 0.22), y = c(1.05, 1.05, 0.0), xend = c(0.22, 0.78, 0.78), yend = c(0.0, 0.0, 0.0),
         style = c("dashed", "dashed", "solid")))
P2 <- mk_panel("管（中介）",
  tibble(name = c("处理 X", "菌群密度 M", "长势 Y"), x = c(0.16, 0.50, 0.84), y = c(0.5, 0.5, 0.5)),
  tibble(from = c("处理 X", "菌群密度 M"), to = c("菌群密度 M", "长势 Y"),
         x = c(0.16, 0.50), y = c(0.5, 0.5), xend = c(0.50, 0.84), yend = c(0.5, 0.5), style = c("solid", "solid")))
P3 <- mk_panel("对撞（选择）",
  tibble(name = c("处理 X", "基线值 B", "被选中 C"), x = c(0.22, 0.78, 0.50), y = c(0.0, 0.0, 1.05)),
  tibble(from = c("处理 X", "基线值 B"), to = c("被选中 C", "被选中 C"),
         x = c(0.22, 0.78), y = c(0.0, 0.0), xend = c(0.50, 0.50), yend = c(1.05, 1.05), style = c("solid", "solid")))
P4 <- mk_panel("对撞的后代",
  tibble(name = c("处理 X", "基线值 B", "被选中 C", "后代 R"),
         x = c(0.16, 0.84, 0.50, 0.50), y = c(0.0, 0.0, 0.72, 1.25)),
  tibble(from = c("处理 X", "基线值 B", "被选中 C"), to = c("被选中 C", "被选中 C", "后代 R"),
         x = c(0.16, 0.84, 0.50), y = c(0.0, 0.0, 0.72), xend = c(0.50, 0.50, 0.50), yend = c(0.72, 0.72, 1.25),
         style = c("solid", "solid", "solid")))
lvl <- c("叉（混杂）", "管（中介）", "对撞（选择）", "对撞的后代")
nodes <- bind_rows(P1$n, P2$n, P3$n, P4$n)
edges <- bind_rows(P1$e, P2$e, P3$e, P4$e)
nodes$panel <- factor(nodes$panel, levels = lvl)
edges$panel <- factor(edges$panel, levels = lvl)
p1 <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, linetype = style),
               arrow = arrow(length = unit(0.18, "cm")), colour = "#167d80", linewidth = 0.8) +
  geom_label(data = nodes, aes(x, y, label = name), size = 3.0, colour = "#163d48", fill = "#eef5f5") +
  facet_wrap(~ panel, nrow = 2) +
  scale_linetype_manual(values = c(solid = "solid", dashed = "dashed"), guide = "none") +
  coord_cartesian(xlim = c(0.02, 0.98), ylim = c(-0.25, 1.5)) +
  labs(title = "混杂的四种基本形状",
       subtitle = "虚线＝共同原因（本图绘图约定）；管＝中介；对撞＝共同结果；其后代＝携带对撞的信息") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#43585c"),
        plot.margin = margin(t = 30, r = 12, b = 8, l = 12),
        strip.text = element_text(size = 10, face = "bold", margin = margin(b = 4)))
save_fig("01-四种结构.png", p1, 7.4, 4.8)

# 图 2：叉——调整批次前后
p2 <- ggplot(fork_est, aes(写法, 处理效应估计)) +
  geom_col(fill = c("#e8b4a8", "#167d80"), width = 0.55) +
  geom_hline(yintercept = 3.0, linetype = "dashed", colour = "#c1462c") +
  geom_text(aes(label = sprintf("%.2f", 处理效应估计)), vjust = -0.5, size = 4) +
  labs(title = "叉（混杂）：不调整批次会把效应估偏",
       subtitle = "红线是真实处理效应 3.0（本例中批次 2 比批次 1 高 2.5，且批次 2 偏对照）",
       x = NULL, y = "处理效应估计") +
  ylim(0, 5.2)
save_fig("02-叉.png", p2, 7.0, 4.2)

# 图 3：管——总效应与直接效应
p3 <- ggplot(pipe_est, aes(写法, 处理效应估计)) +
  geom_col(fill = c("#167d80", "#e8b4a8"), width = 0.55) +
  geom_text(aes(label = sprintf("%.2f", 处理效应估计)), vjust = -0.5, size = 4) +
  labs(title = "管（中介）：控制中介会把总效应带走",
       subtitle = "本例生成机制：总效应真值 3、直接效应真值 0；控制中介后得到的 1.4807 是有限数据的估计偏离",
       x = NULL, y = "处理效应估计") +
  ylim(0, 4.4)
save_fig("03-管.png", p3, 7.0, 4.2)

# 图 4：对撞——全体与选中子集的散点
coll_long <- coll_df %>% mutate(组 = ifelse(选中 == 1, "被选中", "未被选中"))
p4 <- ggplot(coll_long, aes(基线值, 处理)) +
  geom_point(aes(colour = 组), size = 3.4) +
  geom_smooth(data = filter(coll_long, 组 == "被选中"), method = "lm", se = FALSE, colour = "#c1462c", linewidth = 0.9) +
  geom_smooth(method = "lm", se = FALSE, colour = "#9aa8a8", linetype = "dashed", linewidth = 0.9) +
  labs(title = "对撞（选择）：只看被选中的样本，就出现了假的负相关",
       subtitle = paste0("全体相关 = ", r_coll_all, "；只看被选中的 5 个 = ", r_coll_sel,
                         "（灰虚线：全体趋势；红实线：被选中子集）"),
       x = "基线值", y = "处理（1=接种，0=对照）", colour = NULL)
save_fig("04-对撞.png", p4, 7.4, 4.4)

# 图 5：后代——控制处理后变量
p5 <- ggplot(desc_est, aes(范围, 处理与基线值的相关)) +
  geom_col(fill = c("#e8b4a8", "#c1462c", "#167d80"), width = 0.55) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#9aa8a8") +
  geom_text(aes(label = sprintf("%.4f", 处理与基线值的相关)), vjust = -0.5, size = 3.8) +
  labs(title = "没控制那个节点，却控制了它的后代",
       subtitle = "只看「有报告」的单元＝通过后代间接按入选做了筛选；相关仍被拉动",
       x = NULL, y = "处理与基线值的相关") +
  ylim(min(desc_est$处理与基线值的相关) * 1.35, max(desc_est$处理与基线值的相关) * 1.2)
save_fig("05-后代.png", p5, 7.0, 4.2)

cat("\n完成：5 张图 + 10 份 CSV 已写入\n")
