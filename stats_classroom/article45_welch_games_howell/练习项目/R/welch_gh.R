# 0. 准备 ----
# 打开统计课堂45.Rproj，在同一新会话按节运行；也可Source整文件。
# 教学模拟数据及字段说明见README；UTF-8；只依赖R自带stats与graphics。
alpha <- 0.05
confidence <- 1 - alpha
unit_scale <- 1  # 第二遍改为0.133322，把mmHg换成kPa后完整重算。
save_outputs <- identical(Sys.getenv("SAVE_OUTPUTS"), "true")
unit_label <- if (unit_scale == 1) "mmHg" else if (unit_scale == 0.133322) "kPa" else "换算单位"
# 项目根目录有data；若当前在R子目录，则向上一层找。
project_dir <- if (dir.exists("data")) "." else ".."
output_dir <- file.path(project_dir, "运行结果", "R")
figure_dir <- file.path(project_dir, "图形", "R")

# 1. 数据：每行一人，先看字段与人数 ----
h <- read.csv(file.path(project_dir, "data", "hetero_normal.csv"))
h$value <- h$value * unit_scale  # 单位只在这里转换一次。
h$group <- factor(h$group, levels = c("G1", "G2", "G3"))
print(head(h))
print(table(h$group))
# split按方案分开，sapply对每组执行同一种计算；sd使用n-1分母。
groups <- split(h$value, h$group)
names_g <- names(groups)
ni <- sapply(groups, length)
mi <- sapply(groups, mean)
si <- sapply(groups, sd)
k <- length(groups)
if (any(ni < 2) || any(!is.finite(si)) || any(si <= 0)) {
  stop("每组至少2个完整观测且样本方差必须为正。")
}

# 2. 均数精度与Welch整体比较 ----
v <- si^2 / ni
w <- 1 / v
W <- sum(w)
weighted_mean <- sum(w * mi) / W
summary_tab <- data.frame(group = names_g, n = ni, mean = mi, sd = si,
                          variance = si^2, v = v, w = w, row.names = NULL)
print(summary_tab)
# A是精度加权的均数偏离，C将权重比例与各组方差的自由度汇总。
A <- sum(w * (mi - weighted_mean)^2) / (k - 1)
C <- sum((1 - w / W)^2 / (ni - 1))
f_welch <- A / (1 + 2 * (k - 2) * C / (k^2 - 1))
df1 <- k - 1
df2 <- (k^2 - 1) / (3 * C)
# pf的lower.tail=FALSE取近似F分布的右尾。
p_welch <- pf(f_welch, df1, df2, lower.tail = FALSE)
# qf输入左侧累计概率1-alpha，返回右尾alpha对应的临界F；小数自由度直接代入。
f_critical <- qf(confidence, df1, df2)
welch_tab <- data.frame(weighted_mean, W, A, C, F = f_welch, df1, df2,
                        F_critical = f_critical, P = p_welch)
print(welch_tab)
# 公式指定响应与方案；FALSE明确使用Welch，输出F、自由度及P。
welch_function <- oneway.test(value ~ group, data = h, var.equal = FALSE)
print(welch_function)
verdict <- if (p_welch < alpha) "拒绝" else "未拒绝"
cat("在α=", alpha, "下", verdict, "三组总体均数相等；整体检验不指出具体组对。\n", sep = "")

# 3. GH：后组减前组，每对自己的标准误与近似自由度 ----
rows <- list()
row_id <- 0
for (i in seq_len(k - 1)) {
  for (j in seq.int(i + 1, k)) {
    row_id <- row_id + 1
    difference <- mi[j] - mi[i]
    se <- sqrt(v[i] + v[j])
    # i、j只表示当前两个组；自由度不能使用整体Welch的df2。
    nu <- (v[i] + v[j])^2 / (v[i]^2 / (ni[i] - 1) + v[j]^2 / (ni[j] - 1))
    # 这就是当前组对的Welch t；GH把其绝对值转成q尺度。
    t_statistic <- difference / se
    q <- sqrt(2) * abs(t_statistic)
    # k=2只有一对，q的右尾与Welch t双侧P对应；本题GH仍用k=3。
    p_single <- 2 * pt(-abs(t_statistic), df = nu)
    p_q_two <- ptukey(q, nmeans = 2, df = nu, lower.tail = FALSE)
    # k一直是全部组数；P为学生化极差右尾，已包含多重比较调整。
    p_adj <- ptukey(q, nmeans = k, df = nu, lower.tail = FALSE)
    # qtukey输入置信水平得到q临界值，除根号2还原差值SE尺度。
    q_critical <- qtukey(confidence, nmeans = k, df = nu)
    half_width <- q_critical * se / sqrt(2)
    rows[[row_id]] <- data.frame(group1 = names_g[i], group2 = names_g[j],
                                difference, SE = se, df = nu, q, q_critical,
                                CI_low = difference - half_width,
                                CI_high = difference + half_width, P_adj = p_adj,
                                t_statistic = t_statistic, P_single = p_single, P_q_k2 = p_q_two,
                                row.names = NULL)
  }
}
gh_tab <- do.call(rbind, rows)
cat("P_single为单对Welch双侧P，P_q_k2用于两组q尺度对照，本题全部三对的调整结果为P_adj。\n")
print(gh_tab)
# 判断从本轮结果生成，改变数据或单位后不会打印原例固定结论。
for (r in seq_len(nrow(gh_tab))) {
  verdict <- if (gh_tab$P_adj[r] < alpha) "支持总体均数有差异" else "证据不足以确认总体均数有差异"
  cat(sprintf("%s−%s: 差值=%.4f %s，近似%.0f%%同时CI=[%.4f, %.4f] %s，调整P=%.6g；%s。\n",
              gh_tab$group2[r], gh_tab$group1[r], gh_tab$difference[r], unit_label,
              100 * confidence, gh_tab$CI_low[r], gh_tab$CI_high[r], unit_label,
              gh_tab$P_adj[r], verdict))
}

# 4. 图形：个体和均数，以及以0为参照的三对区间 ----
# 选择本机中文字体；其他系统需安装对应字体，不能靠UTF-8代替字体。
plot_font <- switch(Sys.info()[["sysname"]],
                    Darwin = "Heiti TC", Windows = "Microsoft YaHei",
                    "Noto Sans CJK SC")
if (.Platform$OS.type == "windows") {
  windowsFonts(Chinese = windowsFont(plot_font))
  plot_font <- "Chinese"
}
# 将重复排版放入自编绘图函数，便于在屏幕和PNG使用相同内容。
draw_raw <- function() {
  par(family = plot_font, mar = c(5.1, 5.0, 4.1, 2.1))
  plot(NA, xlim = c(0.5, k + 0.5), ylim = range(h$value), xaxt = "n",
       xlab = "降压方案", ylab = paste0("收缩压下降值（", unit_label, "）"),
       main = "个体下降值与各组均数")
  for (i in seq_len(k)) {
    points(i + seq(-0.15, 0.15, length.out = ni[i]), groups[[i]], pch = 16, col = i + 1)
  }
  points(seq_len(k), mi, pch = 18, cex = 1.5)
  axis(1, at = seq_len(k), labels = names_g)
  legend("topleft", "组均数", pch = 18, bty = "n")
}
draw_ci <- function() {
  par(family = plot_font, mar = c(5.1, 6.5, 4.1, 2.1))
  positions <- rev(seq_len(nrow(gh_tab)))
  plot(gh_tab$difference, positions, xlim = range(c(0, gh_tab$CI_low, gh_tab$CI_high)),
       ylim = c(0.5, nrow(gh_tab) + 0.5), yaxt = "n", pch = 16, ylab = "",
       xlab = paste0("平均下降值之差（", unit_label, "）"),
       main = sprintf("均数差与近似%.0f%%同时区间", 100 * confidence))
  mtext("方案比较", side = 2, line = 5.0)
  segments(gh_tab$CI_low, positions, gh_tab$CI_high, positions)
  axis(2, at = positions, labels = paste0(gh_tab$group2, "−", gh_tab$group1), las = 1)
  abline(v = 0, lty = 2, col = "grey")
}
# IDE中直接看图；命令行用保存开关输出PNG，避免默认PDF设备的字体限制。
if (interactive()) {
  draw_raw()
  draw_ci()
}

# 5. 可选保存：IDE首次练习默认关闭；开关在第0节 ----
if (save_outputs) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_tab, file.path(output_dir, "group_summary.csv"), row.names = FALSE)
  write.csv(welch_tab, file.path(output_dir, "welch_intermediate.csv"), row.names = FALSE)
  write.csv(gh_tab, file.path(output_dir, "games_howell.csv"), row.names = FALSE)
  png(file.path(figure_dir, "raw_points_means.png"), width = 1100, height = 800, res = 160)
  draw_raw()
  dev.off()
  png(file.path(figure_dir, "games_howell_ci.png"), width = 1100, height = 800, res = 160)
  draw_ci()
  dev.off()
}
