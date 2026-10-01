# ============================================================
# 精读 06｜好控制与坏控制：后门准则、对撞偏倚、精度寄生虫、Table 2 谬误
# 对应：McElreath 2023 第 06 讲 Good and Bad Controls；教材第 6 章
# 用法： Rscript 代码/篇06_好控制与坏控制.R
#
# 约定：与 05 篇同一套工程约定——显式最小二乘（解正规方程）、显式相关与方差膨胀，
#       全部确定性；三语言输出逐字节一致。
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
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
  res <- y - X %*% b
  list(coef = as.numeric(b), sigma_hat = sqrt(sum(res^2) / (length(y) - ncol(X))),
       r2 = 1 - sum(res^2) / sum((y - mean(y))^2))
}
r4 <- function(x) round(as.numeric(x), 4)

# ---------- 0) 工具：把扰动正交化，让最小二乘能精确复现真值 ----------
orth_resid <- function(v, X) {
  v - X %*% solve(t(X) %*% X, t(X) %*% v)
}

# ---------- 1) 我们的案例：12 个实验单元 ----------
n <- 12
treat_x <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
X0 <- cbind(1, treat_x)
# 混杂 z 的扰动：与截距、处理都正交 → z = 0.7·处理 + 扰动，之后 OLS 能精确还原
resid_z <- orth_resid(c(0.4, -0.3, 0.2, -0.5, 0.6, -0.2, 0.3, -0.4, 0.5, -0.6, 0.1, 0.4), X0)
z_conf  <- 0.7 * treat_x + resid_z
# 结果：真实处理效应 3.0、混杂效应 2.0；扰动与 [1, 处理, 混杂] 正交
resid_y <- orth_resid(c(0.5, -0.4, 0.3, -0.6, 0.2, -0.3, 0.6, -0.5, 0.4, -0.2, 0.1, 0.5),
                      cbind(1, treat_x, z_conf))
y_out   <- 3.0 * treat_x + 2.0 * z_conf + resid_y
# 对撞 w：被处理与结果共同影响（治疗后观察到的指标）；扰动与 [1, 处理, 结果] 正交
resid_w <- orth_resid(c(0.3, 0.5, -0.4, 0.2, -0.5, 0.4, -0.2, 0.6, -0.3, 0.1, -0.6, 0.4),
                      cbind(1, treat_x, y_out))
w_coll  <- 1.0 * y_out + 0.5 * treat_x + resid_w

dat <- tibble(单元 = 1:n, 处理 = treat_x, 混杂 = round(z_conf, 4),
              结果 = round(y_out, 4), 对撞 = round(w_coll, 4))
write_csv_fixed(dat, file.path(RES, "01_我们的案例.csv"))

# 三行关键回归：不控制 / 只控制混杂（好的控制）/ 控制混杂与对撞（坏的控制）
m_none <- ols_ls(cbind(1, treat_x), y_out)
m_good <- ols_ls(cbind(1, treat_x, z_conf), y_out)
m_bad  <- ols_ls(cbind(1, treat_x, z_conf, w_coll), y_out)
ctrl_tab <- tibble(
  写法 = c("不控制任何变量", "只控制混杂 z（好的控制）", "同时控制混杂 z 与对撞 w（坏的控制）"),
  处理效应估计 = r4(c(m_none$coef[2], m_good$coef[2], m_bad$coef[2])),
  真值 = 3.0
)
write_csv_fixed(ctrl_tab, file.path(RES, "02_好控制与坏控制.csv"))
print(ctrl_tab)

# ---------- 2) 对撞偏倚：处理对结果真的没有效应时 ----------
y_zero  <- orth_resid(c(0.5, -0.4, 0.3, -0.6, 0.2, -0.3, 0.6, -0.5, 0.4, -0.2, 0.1, 0.5), X0)
resid_w2 <- orth_resid(c(0.2, -0.3, 0.4, -0.2, 0.5, -0.4, 0.3, -0.5, 0.1, 0.4, -0.1, 0.3),
                       cbind(1, treat_x, y_zero))
w_coll2  <- 1.2 * treat_x + 1.5 * y_zero + resid_w2
z0 <- ols_ls(cbind(1, treat_x), y_zero)
z1 <- ols_ls(cbind(1, treat_x, w_coll2), y_zero)
coll_tab <- tibble(
  写法 = c("不控制对撞（结果 ~ 处理）", "控制对撞（结果 ~ 处理 + 对撞）"),
  处理效应估计 = r4(c(z0$coef[2], z1$coef[2])),
  真值 = 0.0
)
write_csv_fixed(coll_tab, file.path(RES, "03_对撞偏倚.csv"))
print(coll_tab)

# ---------- 3) 精度寄生虫：与处理相关但不帮助预测结果的变量 ----------
resid_p <- orth_resid(c(0.1, -0.2, 0.3, -0.1, 0.2, -0.3, 0.1, 0.2, -0.2, 0.3, -0.1, 0.2),
                      cbind(1, treat_x, y_out))
parasite <- 0.9 * treat_x + resid_p        # 与处理强相关，但对结果没有额外解释力
r2_xz <- ols_ls(cbind(1, parasite), treat_x)$r2
se_a <- ols_ls(cbind(1, treat_x), y_out)$sigma_hat / sqrt(sum((treat_x - mean(treat_x))^2))
se_b <- ols_ls(cbind(1, treat_x, parasite), y_out)$sigma_hat /
  sqrt(sum((treat_x - mean(treat_x))^2) * (1 - r2_xz))
par_tab <- tibble(
  量 = c("处理与寄生虫的决定系数 R²", "不控制时的系数标准误", "控制寄生虫后的系数标准误",
         "标准误放大倍数", "理论值 1/sqrt(1-R²)"),
  值 = r4(c(r2_xz, se_a, se_b, se_b / se_a, 1 / sqrt(1 - r2_xz)))
)
write_csv_fixed(par_tab, file.path(RES, "04_精度寄生虫.csv"))
print(par_tab)

# ---------- 4) 偏倚放大：存在不可测混杂时，控制"对撞"会比不控制更糟 ----------
# U 是未被观测到的混杂（处理与结果都被它影响），z 是已观测到的混杂
u_hidden <- 0.6 * treat_x + orth_resid(c(0.2, -0.4, 0.5, -0.1, 0.3, -0.5, 0.4, -0.2, 0.1, 0.5, -0.3, 0.2), X0)
resid_amp <- orth_resid(c(0.4, -0.3, 0.2, -0.5, 0.6, -0.1, 0.3, -0.4, 0.5, -0.6, 0.1, 0.3),
                        cbind(1, treat_x, z_conf, u_hidden))
y_amp <- 3.0 * treat_x + 2.0 * z_conf + 1.5 * u_hidden + resid_amp
w_amp <- 1.0 * y_amp + 0.5 * treat_x +
  orth_resid(c(0.3, 0.5, -0.4, 0.2, -0.5, 0.4, -0.2, 0.6, -0.3, 0.1, -0.6, 0.4),
             cbind(1, treat_x, y_amp))
a0 <- ols_ls(cbind(1, treat_x), y_amp)$coef[2]                 # 什么都不控制
a1 <- ols_ls(cbind(1, treat_x, z_conf), y_amp)$coef[2]         # 只控制已观测的混杂
a2 <- ols_ls(cbind(1, treat_x, z_conf, w_amp), y_amp)$coef[2]  # 再控制"后果 w"
amp_hidden <- tibble(
  写法 = c("不控制任何变量", "只控制混杂 z", "控制混杂 z 与后果 w"),
  处理效应估计 = r4(c(a0, a1, a2)),
  真值 = 3.0,
  偏离真值的量 = r4(abs(c(a0, a1, a2) - 3.0))
)
write_csv_fixed(amp_hidden, file.path(RES, "05_偏倚放大.csv"))
print(amp_hidden)

# ---------- 图 ----------
mk_panel <- function(panel_name, nodes, edges) {
  list(n = tibble(panel = panel_name, name = nodes$name, x = nodes$x, y = nodes$y),
       e = mutate(as_tibble(edges), panel = panel_name))
}
A <- mk_panel("好的控制：z 是混杂",
  tibble(name = c("混杂 z", "处理 x", "结果 y"), x = c(0.50, 0.18, 0.82), y = c(1.05, 0.0, 0.0)),
  tibble(x = c(0.50, 0.50, 0.18), y = c(1.05, 1.05, 0.0), xend = c(0.18, 0.82, 0.82), yend = c(0.0, 0.0, 0.0),
         style = c("dashed", "dashed", "solid")))
B <- mk_panel("坏的控制：w 是对撞",
  tibble(name = c("处理 x", "结果 y", "对撞 w"), x = c(0.18, 0.82, 0.50), y = c(0.0, 0.0, 1.05)),
  tibble(x = c(0.18, 0.82, 0.50), y = c(0.0, 0.0, 1.05), xend = c(0.50, 0.50, 0.50), yend = c(1.05, 1.05, 1.05),
         style = c("solid", "solid", "solid")))
C <- mk_panel("精度寄生虫：只影响结果",
  tibble(name = c("处理 x", "结果 y", "寄生虫 p"), x = c(0.18, 0.82, 0.50), y = c(0.5, 0.5, 1.05)),
  tibble(x = c(0.18, 0.82, 0.50), y = c(0.5, 0.5, 1.05), xend = c(0.82, 0.82, 0.82), yend = c(0.5, 0.5, 0.5),
         style = c("solid", "solid", "solid")))
lvl <- c("好的控制：z 是混杂", "坏的控制：w 是对撞", "精度寄生虫：只影响结果")
nodes <- bind_rows(A$n, B$n, C$n); edges <- bind_rows(A$e, B$e, C$e)
nodes$panel <- factor(nodes$panel, levels = lvl); edges$panel <- factor(edges$panel, levels = lvl)
p1 <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, linetype = style),
               arrow = arrow(length = unit(0.18, "cm")), colour = "#167d80", linewidth = 0.8) +
  geom_label(data = nodes, aes(x, y, label = name), size = 3.0, colour = "#163d48", fill = "#eef5f5") +
  facet_wrap(~ panel, nrow = 1) +
  scale_linetype_manual(values = c(solid = "solid", dashed = "dashed"), guide = "none") +
  coord_cartesian(xlim = c(0.0, 1.0), ylim = c(-0.3, 1.35)) +
  labs(title = "三种「要不要控制它」的情形",
       subtitle = "虚线＝控制它是对的；实线指向它＝控制它会出错（对撞／寄生虫）") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#43585c"),
        plot.margin = margin(t = 30, r = 12, b = 8, l = 12),
        strip.text = element_text(size = 9.5, face = "bold", margin = margin(b = 4)))
save_fig("01-三种情形.png", p1, 7.6, 3.6)

p2 <- ggplot(ctrl_tab, aes(写法, 处理效应估计)) +
  geom_col(fill = c("#e8b4a8", "#167d80", "#c1462c"), width = 0.55) +
  geom_hline(yintercept = 3.0, linetype = "dashed", colour = "#163d48") +
  geom_text(aes(label = sprintf("%.2f", 处理效应估计)), vjust = -0.5, size = 4) +
  labs(title = "好控制、不控制、坏控制：同一份数据三种答案",
       subtitle = "虚线是真实处理效应 3.0",
       x = NULL, y = "处理效应估计") +
  ylim(0, 5.6)
save_fig("02-三种答案.png", p2, 7.4, 4.2)

p3 <- ggplot(coll_tab, aes(写法, 处理效应估计)) +
  geom_col(fill = c("#167d80", "#c1462c"), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#163d48") +
  geom_text(aes(label = sprintf("%.2f", 处理效应估计)), vjust = -0.5, size = 4, colour = "#163d48") +
  labs(title = "处理真的没有效应时，控制对撞会凭空造出一个",
       subtitle = "真值 0：不控制时估计接近 0，把对撞放进来就出现明显的负值",
       x = NULL, y = "处理效应估计") +
  ylim(-1.2, 0.6)
save_fig("03-对撞偏倚.png", p3, 7.4, 4.2)

p4 <- ggplot(par_tab[2:3, ], aes(量, 值)) +
  geom_col(fill = c("#167d80", "#c1462c"), width = 0.5) +
  geom_text(aes(label = sprintf("%.4f", 值)), vjust = -0.5, size = 4) +
  labs(title = "精度寄生虫：控制它不会减小偏倚，只会放大方差",
       subtitle = paste0("系数标准误从 ", par_tab$值[2], " 变成 ", par_tab$值[3],
                         "（放大 ", par_tab$值[4], " 倍，理论值 1/sqrt(1-R²) = ", par_tab$值[5], "）"),
       x = NULL, y = "系数标准误") +
  ylim(0, max(par_tab$值[2:3]) * 1.35)
save_fig("04-精度寄生虫.png", p4, 7.4, 4.2)

p5 <- ggplot(amp_hidden, aes(写法, 处理效应估计)) +
  geom_col(fill = c("#e8b4a8", "#167d80", "#c1462c"), width = 0.55) +
  geom_hline(yintercept = 3.0, linetype = "dashed", colour = "#163d48") +
  geom_text(aes(label = sprintf("%.2f", 处理效应估计)), vjust = -0.5, size = 4) +
  labs(title = "有不可测混杂时：控制「后果」会比不控制更糟",
       subtitle = "虚线是真值 3.0；三根柱依次是：什么都不控制、只控制已观测混杂 z、再控制后果 w",
       x = NULL, y = "处理效应估计") +
  ylim(0, max(amp_hidden$处理效应估计) * 1.25)
save_fig("05-偏倚放大.png", p5, 7.4, 4.2)

cat("\n完成：5 张图 + 5 份 CSV 已写入\n")
