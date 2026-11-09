# 配对样本率差：三张教学图（基础 R，无额外依赖）
# 本脚本可由 Rscript 执行，也可在 RStudio 中 Source；输出均写入项目目录。

#---- 项目路径与字体 ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root_candidates <- unique(c(script_dir, getwd(), dirname(script_dir)))
root <- root_candidates[vapply(root_candidates, function(x) {
  dir.exists(file.path(x, "数据")) && dir.exists(file.path(x, "运行结果"))
}, logical(1))][1]
if (is.na(root) || !nzchar(root)) stop("找不到项目根目录（需要同时有‘数据’和‘运行结果’）。")
fig_dir <- file.path(dirname(root), "文章配图")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
result_dir <- file.path(root, "运行结果", "R")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# Cairo 使用系统 CJK 字体，避免中文标题和坐标轴乱码；可按本机字体替换。
# 不把 fc-match 或其他外部命令作为运行前提。
font_family <- switch(Sys.info()[["sysname"]],
                      Darwin = "Arial Unicode MS",
                      Windows = "Microsoft YaHei",
                      Linux = "Noto Sans CJK SC",
                      "sans")

#---- 通用分布图函数 ----
read_distribution <- function(name) {
  path <- file.path(result_dir, paste0("distribution_", name, ".csv"))
  if (!file.exists(path)) stop("找不到运行结果：", path)
  read.csv(path, check.names = FALSE)
}

draw_distribution <- function(dat, title, limits) {
  shown <- dat$D >= limits[1] & dat$D <= limits[2]
  d <- dat[shown, , drop = FALSE]
  if (!nrow(d)) stop("指定横轴范围内没有分布数据。")
  ymax <- 1.25 * max(d$exact_probability, d$normal_cell_probability)
  plot(d$D, d$exact_probability, type = "n", xlim = limits, ylim = c(0, ymax),
       xlab = "率差 D=k/n（后测减前测）", ylab = "单个取值的概率", main = title)
  abline(h = pretty(c(0, ymax)), col = "grey90", lty = 1)
  width <- min(diff(sort(unique(d$D)))) * 0.82
  rect(d$D - width / 2, 0, d$D + width / 2, d$exact_probability,
       col = adjustcolor("#4C78A8", alpha.f = 0.8), border = NA)
  points(d$D, d$normal_cell_probability, pch = 16, cex = 0.65, col = "#E45756")
  legend("topright", legend = c("精确概率", "正态单格概率"),
         fill = c(adjustcolor("#4C78A8", 0.8), NA), border = c(NA, NA),
         pch = c(NA, 16), col = c("#4C78A8", "#E45756"), bty = "n", pt.cex = 0.8)
}

#---- 图1：精确分布与正态单格近似 ----
png(file.path(fig_dir, "图1_配对率差精确与近似.png"), width = 6.6, height = 8.5,
    units = "in", res = 190, type = "cairo", family = font_family)
par(mfrow = c(2, 1), mar = c(4.5, 5, 3.2, 1.2), oma = c(0, 0, 0, 0))
draw_distribution(read_distribution("normal_n20"), "给定变化概率，n=20", c(-0.65, 0.25))
draw_distribution(read_distribution("normal_n1000"), "同一变化概率，n=1000（横轴局部放大）", c(-0.27, -0.13))
dev.off()
grDevices::cairo_pdf(file.path(fig_dir, "图1_配对率差精确与近似.pdf"), width = 6.6, height = 8.5,
                     family = font_family)
par(mfrow = c(2, 1), mar = c(4.5, 5, 3.2, 1.2))
draw_distribution(read_distribution("normal_n20"), "给定变化概率，n=20", c(-0.65, 0.25))
draw_distribution(read_distribution("normal_n1000"), "同一变化概率，n=1000（横轴局部放大）", c(-0.27, -0.13))
dev.off()

#---- 图2：变化稀少时的局部形状 ----
png(file.path(fig_dir, "图2_变化稀少.png"), width = 6.6, height = 4.8,
    units = "in", res = 190, type = "cairo", family = font_family)
par(mar = c(4.5, 5, 3.2, 1.2))
draw_distribution(read_distribution("rare_n1000"), "n=1000，平均只有3人发生变化", c(-0.009, 0.007))
dev.off()
grDevices::cairo_pdf(file.path(fig_dir, "图2_变化稀少.pdf"), width = 6.6, height = 4.8,
                     family = font_family)
par(mar = c(4.5, 5, 3.2, 1.2))
draw_distribution(read_distribution("rare_n1000"), "n=1000，平均只有3人发生变化", c(-0.009, 0.007))
dev.off()

#---- 图3：Wang 下限的两步网格搜索 ----
# 141 组真实尾部计数逐项求和；不读取 Python 曲线，也不按 D 简化排序。
counts <- read.csv(file.path(root, "数据", "wang_lower_tail_n20.csv"))[, c("n_plus", "n_zero", "n_minus")]
counts <- as.matrix(counts)
n <- sum(counts[1, ])
stopifnot(all(rowSums(counts) == n), nrow(counts) == 141)
theta <- -0.45
upper <- 1 - abs(theta)
log_coeff <- lgamma(n + 1) - rowSums(lgamma(counts + 1))
xlogy <- function(k, p) ifelse(k == 0, 0, k * log(p))
tail_probability <- function(p_t) {
  vapply(p_t, function(p) {
    probs <- c((1 - p + theta) / 2, p, (1 - p - theta) / 2)
    if (min(probs) < -1e-12) stop("p_T 超出合法范围。")
    probs <- pmax(probs, 0)
    log_terms <- log_coeff + rowSums(vapply(seq_len(3), function(j) xlogy(counts[, j], probs[j]), numeric(nrow(counts))))
    sum(exp(log_terms))
  }, numeric(1))
}

coarse <- seq(0, 0.55, by = 0.055)
fine <- seq(0.44, 0.55, by = 0.0055)
reference <- seq(0, upper, length.out = 4001)
y_coarse <- tail_probability(coarse)
y_fine <- tail_probability(fine)
y_reference <- tail_probability(reference)
coarse_best <- which.max(y_coarse)
fine_best <- which.max(y_fine)
stopifnot(abs(coarse[coarse_best] - 0.495) < 1e-12,
          abs(y_coarse[coarse_best] - 0.02462737352027) < 1e-12,
          abs(fine[fine_best] - 0.4895) < 1e-12,
          abs(y_fine[fine_best] - 0.02463073857116) < 1e-12)
write_grid <- function(name, x, y) write.csv(data.frame(theta = theta, p_t = x, tail_probability = y),
                                             file.path(result_dir, paste0("wang_grid_", name, ".csv")), row.names = FALSE)
write_grid("coarse", coarse, y_coarse)
write_grid("fine", fine, y_fine)
write_grid("reference", reference, y_reference)
summary <- list(theta = theta, tail_count = nrow(counts), p_t_range = c(0, upper),
                coarse_points = length(coarse), fine_points = length(fine),
                coarse_step = coarse[2] - coarse[1], fine_step = fine[2] - fine[1],
                coarse_peak_p_t = coarse[coarse_best], coarse_peak_probability = y_coarse[coarse_best],
                fine_range = c(0.44, 0.55), fine_peak_p_t = fine[fine_best],
                fine_peak_probability = y_fine[fine_best])
json_scalar <- function(x) if (length(x) == 1) format(x, digits = 16, scientific = FALSE) else paste0("[", paste(vapply(x, json_scalar, character(1)), collapse = ","), "]")
json_text <- paste0("{\n  \"theta\":", json_scalar(summary$theta), ",\n  \"tail_count\":", summary$tail_count,
                    ",\n  \"p_t_range\":", json_scalar(summary$p_t_range), ",\n  \"coarse_points\":", summary$coarse_points,
                    ",\n  \"fine_points\":", summary$fine_points, ",\n  \"coarse_step\":", json_scalar(summary$coarse_step),
                    ",\n  \"fine_step\":", json_scalar(summary$fine_step), ",\n  \"coarse_peak_p_t\":", json_scalar(summary$coarse_peak_p_t),
                    ",\n  \"coarse_peak_probability\":", json_scalar(summary$coarse_peak_probability),
                    ",\n  \"fine_range\":", json_scalar(summary$fine_range), ",\n  \"fine_peak_p_t\":", json_scalar(summary$fine_peak_p_t),
                    ",\n  \"fine_peak_probability\":", json_scalar(summary$fine_peak_probability), "\n}\n")
writeLines(json_text, file.path(result_dir, "wang_grid_summary.json"), useBytes = TRUE)

draw_grid <- function() {
  par(mfrow = c(2, 1), mar = c(4.5, 5, 3.3, 1.2))
  plot(reference, y_reference, type = "n", col = "#6f7c88", lwd = 1.7,
       xlim = c(0, upper), ylim = c(0, max(y_reference) * 1.28),
       main = "上：固定候选 θ = -0.45，先在整个范围粗搜",
       xlab = expression(paste("两次结果相同的总体概率 ", p[T])), ylab = "单侧尾部概率")
  usr <- par("usr")
  rect(0.44, usr[3], 0.55, usr[4], col = adjustcolor("#df6b25", 0.12), border = NA)
  grid(col = "grey88")
  lines(reference, y_reference, col = "#6f7c88", lwd = 1.7)
  points(coarse, y_coarse, pch = 16, col = "#1266ab")
  points(coarse[coarse_best], y_coarse[coarse_best], pch = 1, col = "#df6b25", cex = 1.8, lwd = 2)
  abline(v = c(0.44, 0.55), col = adjustcolor("#df6b25", 0.3), lty = 3)
  legend("bottomleft", c("密集计算的参考曲线", "粗网格：11个试算点"), col = c("#6f7c88", "#1266ab"), lty = c(1, NA), pch = c(NA, 16), bty = "n")
  text(0.38, max(y_reference) * 1.05, bquote(atop("粗网格最高点", p[T] == .(sprintf("%.3f", coarse[coarse_best])))), col = "#a42d2d")
  local_x <- seq(0.44, 0.55, length.out = 1601)
  local_y <- tail_probability(local_x)
  spread <- max(y_fine) - min(y_fine)
  plot(local_x, local_y, type = "l", col = "#6f7c88", lwd = 1.7,
       xlim = c(0.44, 0.55), ylim = c(max(0, min(y_fine) - 0.18 * spread), max(y_fine) + 0.48 * spread),
       main = "下：放大 [0.440, 0.550]，再细搜",
       xlab = expression(paste("两次结果相同的总体概率 ", p[T])), ylab = "单侧尾部概率")
  grid(col = "grey88")
  points(fine, y_fine, pch = 16, col = "#df6b25")
  points(coarse[coarse >= 0.44 & coarse <= 0.55], y_coarse[coarse >= 0.44 & coarse <= 0.55], pch = 16, col = "#1266ab", cex = 1.2)
  points(fine[fine_best], y_fine[fine_best], pch = 1, col = "#a42d2d", cex = 1.8, lwd = 2)
  legend("bottomleft", c("密集计算的参考曲线", "细网格：21个试算点", "该区域内粗网格点"), col = c("#6f7c88", "#df6b25", "#1266ab"), lty = c(1, NA, NA), pch = c(NA, 16, 16), bty = "n")
  text(0.497, max(y_fine) + 0.18 * spread, bquote(atop(paste("细网格最高点 ", p[T], "=", .(sprintf("%.4f", fine[fine_best]))), paste("尾部概率=", .(sprintf("%.6f", y_fine[fine_best]))))), col = "#a42d2d")
}
png(file.path(fig_dir, "图3_未知概率的两步网格搜索.png"), width = 7.1, height = 8.8,
    units = "in", res = 210, type = "cairo", family = font_family)
draw_grid()
dev.off()
grDevices::cairo_pdf(file.path(fig_dir, "图3_未知概率的两步网格搜索.pdf"), width = 7.1, height = 8.8,
                     family = font_family)
draw_grid()
dev.off()

message(sprintf("已生成三张 R 图；粗峰 %.4f / %.14f，细峰 %.4f / %.14f。", coarse[coarse_best], y_coarse[coarse_best], fine[fine_best], y_fine[fine_best]))
