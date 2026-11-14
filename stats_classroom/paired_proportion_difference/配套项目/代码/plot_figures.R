# 统计课堂27：配对率差教学图（R/ggplot2）
# 本文件可独立运行，不依赖已有CSV；也可由paired_rate.R自动调用。
# 图形对象默认显示；改save_figures <- TRUE后才保存到运行结果/R/图像目录。

# 问题/条件：定位项目并加载实际绘图库 ----
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("需要安装ggplot2后再运行本文件。")
library(ggplot2)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
if (nzchar(script_path) && !file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
script_dir <- if (nzchar(script_path) && file.exists(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
candidates <- unique(c(getwd(), script_dir, dirname(script_dir)))
valid <- vapply(candidates, function(path) dir.exists(file.path(path, "数据")), logical(1))
if (!any(valid)) stop("找不到项目根目录；请从项目根或代码目录运行。")
root <- normalizePath(candidates[which(valid)[1]], mustWork = TRUE)

# 手动计算：自编卷积生成绘图数据 ----
# 每次循环加入1人；exact_probability是精确PMF，normal_cell_probability是正态中点单格概率。
distribution_for_plot <- function(n, q_plus, q_minus) {
  q_zero <- 1 - q_plus - q_minus
  probability <- 1
  for (person in seq_len(n)) {
    next_probability <- numeric(length(probability) + 2L)
    index <- seq_along(probability)
    next_probability[index] <- next_probability[index] + q_minus * probability
    next_probability[index + 1L] <- next_probability[index + 1L] + q_zero * probability
    next_probability[index + 2L] <- next_probability[index + 2L] + q_plus * probability
    probability <- next_probability
  }
  k <- -n:n
  D <- k / n
  theta <- q_plus - q_minus
  SE_D <- sqrt((q_plus + q_minus - theta^2) / n)
  edges <- (D[-length(D)] + D[-1]) / 2
  normal_cell_probability <- diff(c(0, pnorm(edges, theta, SE_D), 1))
  data.frame(k = k, D = D, exact_probability = probability,
             normal_cell_probability = normal_cell_probability)
}

# 条件：准备可直接查看的绘图数据表 ----
normal20 <- distribution_for_plot(20, 0.05, 0.25)
normal1000 <- distribution_for_plot(1000, 0.05, 0.25)
rare1000 <- distribution_for_plot(1000, 0.001, 0.002)
normal20
normal1000
rare1000

# 手动计算：图形范围内的显示数据 ----
# 柱高是精确概率；点高是同一D的正态单格概率，不是连续密度。
plot_data_20 <- subset(normal20, D >= -0.65 & D <= 0.25)
plot_data_1000 <- subset(normal1000, D >= -0.27 & D <= -0.13)
plot_data_rare <- subset(rare1000, D >= -0.009 & D <= 0.007)
plot_data_20$scenario <- "n=20"
plot_data_1000$scenario <- "n=1000 (local x-axis)"
plot_data_20$bar_width <- 0.04
plot_data_1000$bar_width <- 0.0008
plot_data_normal <- rbind(plot_data_20, plot_data_1000)
plot_data_normal$scenario <- factor(plot_data_normal$scenario, levels = c("n=20", "n=1000 (local x-axis)"))
plot_data_20
plot_data_1000
plot_data_rare

# 调用函数：ggplot2图1，精确分布和正态单格近似 ----
# geom_col画精确概率，geom_point叠加正态单格概率，facet_wrap分场景。
p1 <- ggplot(plot_data_normal, aes(x = D)) +
  geom_col(aes(y = exact_probability, width = bar_width, fill = "Exact PMF"), alpha = 0.8) +
  geom_point(aes(y = normal_cell_probability, colour = "Normal cell probability"), size = 1.4) +
  scale_fill_manual(values = c("Exact PMF" = "#4C78A8"), name = NULL) +
  scale_colour_manual(values = c("Normal cell probability" = "#E45756"), name = NULL) +
  facet_wrap(~scenario, scales = "free", ncol = 1) +
  labs(title = "Exact and normal-cell probabilities for paired rate difference",
       x = "Rate difference D=k/n (post - pre)", y = "Probability at one D value") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
print(p1)  # Source运行时显示图形对象；逐节也可输入p1查看。

# 调用函数：ggplot2图2，变化稀少时的局部形状 ----
# 这里q+=.001、q−=.002，用来观察稀有变化时的近似边界。
p2 <- ggplot(plot_data_rare, aes(x = D)) +
  geom_col(aes(y = exact_probability, fill = "Exact PMF"), alpha = 0.8, width = 0.0008) +
  geom_point(aes(y = normal_cell_probability, colour = "Normal cell probability"), size = 1.4) +
  scale_fill_manual(values = c("Exact PMF" = "#4C78A8"), name = NULL) +
  scale_colour_manual(values = c("Normal cell probability" = "#E45756"), name = NULL) +
  labs(title = "Paired rate difference when changes are rare (n=1000)",
       x = "Rate difference D=k/n (post - pre)", y = "Probability at one D value") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
print(p2)  # 显式显示图形对象。

# 可选保存：学生可改TRUE；保存不会影响前面的计算和显示 ----
save_figures <- FALSE
if (save_figures) {
  fig_dir <- file.path(root, "运行结果", "R", "图像")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(fig_dir, "图1_配对率差精确与近似.png"), p1, width = 6.6, height = 8.5, dpi = 190)
  ggsave(file.path(fig_dir, "图2_变化稀少.png"), p2, width = 6.6, height = 4.8, dpi = 190)
}

# 思考/可选练习：改q+、q−或横轴筛选范围，先预测图形变化再运行。
