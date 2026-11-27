# 用R重画统计课堂45的两幅公众号数据图；仅依赖base R。
# 在RStudio项目根目录运行：source("R/plot_article.R"); draw_article()
# data_file为共享教学模拟CSV，out_dir为图像保存目录，coordinates_file为图2坐标。
draw_article <- function(data_file = "data/hetero_normal.csv", out_dir = "文章配图",
                         coordinates_file = file.path(out_dir, "区间坐标.csv")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(coordinates_file), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(data_file)) stop("找不到共享CSV：", data_file)
  h <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
  names_g <- c("G1", "G2", "G3")
  groups <- lapply(names_g, function(z) h$value[h$group == z])
  ni <- vapply(groups, length, integer(1))
  mi <- vapply(groups, mean, numeric(1))
  v <- vapply(seq_along(groups), function(i) var(groups[[i]]) / ni[i], numeric(1))
  if (any(ni < 2) || any(!is.finite(v)) || any(v <= 0)) stop("各组需至少2个观测且方差为正。")
  colors <- c("#2377B5", "#39A7A0", "#E39437")
  font <- switch(Sys.info()[["sysname"]], Darwin = "Heiti TC",
                 Windows = "Microsoft YaHei", "Noto Sans CJK SC")
  png_open <- function(path, width, height) {
    png(path, width = width, height = height, res = 300, type = "cairo")
  }

  # 图1：个体点与组均数横线。
  png_open(file.path(out_dir, "图01_原始数据.png"), 1098, 1110)
  op <- par(mar = c(4.5, 5.0, 3.2, 1.0), family = font, mgp = c(2.7, 0.8, 0))
  ylim <- range(h$value)
  ylim <- c(floor(ylim[1] - 1), ceiling(ylim[2] + 1))
  plot(NA, xlim = c(0.5, 3.5), ylim = ylim, xaxt = "n", xlab = "降压方案", ylab = "收缩压下降值（mmHg）", main = "个体下降值与各组均数", cex.main = 1.15, cex.lab = 1.05)
  abline(h = pretty(ylim), col = "grey90", lty = 1)
  for (i in seq_along(groups)) {
    xpos <- i + seq(-0.16, 0.16, length.out = ni[i])
    points(xpos, groups[[i]], pch = 16, col = colors[i], cex = 1.0)
    segments(i - 0.27, mi[i], i + 0.27, mi[i], col = "#142F4C", lwd = 3)
  }
  axis(1, at = 1:3, labels = names_g)
  box(bty = "l")
  par(op)
  dev.off()

  # 图2：Games–Howell后组减前组的近似95%同时区间。
  # 每行是一对前、后组的下标：G2−G1、G3−G1、G3−G2。
  pairs <- matrix(c(1, 2, 1, 3, 2, 3), ncol = 2, byrow = TRUE)
  rows <- vector("list", nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1]
    j <- pairs[r, 2]
    se <- sqrt(v[i] + v[j])
    df <- (v[i] + v[j])^2 / (v[i]^2 / (ni[i] - 1) + v[j]^2 / (ni[j] - 1))
    difference <- mi[j] - mi[i]
    half <- qtukey(0.95, nmeans = 3, df = df) * se / sqrt(2)
    rows[[r]] <- data.frame(comparison = paste0(names_g[j], "−", names_g[i]),
                            difference = difference, SE = se, df = df,
                            lower = difference - half, upper = difference + half)
  }
  # 把三对的命名结果表按行合并，保留数值列的类型。
  ci <- do.call(rbind, rows)
  write.csv(ci, coordinates_file, row.names = FALSE, fileEncoding = "UTF-8")
  png_open(file.path(out_dir, "图02_同时区间.png"), 1098, 930)
  op <- par(mar = c(4.5, 5.0, 3.2, 1.0), family = font, mgp = c(2.7, 0.8, 0))
  xrange <- range(c(0, ci$lower, ci$upper))
  pad <- max(diff(xrange) * 0.08, 0.5)
  plot(ci$difference, rev(seq_len(nrow(ci))), type = "n", xlim = c(xrange[1] - pad, xrange[2] + pad), ylim = c(0.5, 3.5), yaxt = "n", xlab = "平均下降值之差（mmHg）", ylab = "", main = "均数差与95%同时区间", cex.main = 1.15, cex.lab = 1.05)
  abline(v = 0, col = "#687B8D", lty = 2, lwd = 1.2)
  for (r in seq_len(nrow(ci))) {
    y <- nrow(ci) - r + 1
    segments(ci$lower[r], y, ci$upper[r], y, col = colors[r], lwd = 3)
    segments(ci$lower[r], y - 0.10, ci$lower[r], y + 0.10, col = colors[r], lwd = 2)
    segments(ci$upper[r], y - 0.10, ci$upper[r], y + 0.10, col = colors[r], lwd = 2)
    points(ci$difference[r], y, pch = 16, col = colors[r], cex = 1.15)
  }
  axis(2, at = rev(seq_len(nrow(ci))), labels = ci$comparison, las = 1)
  box(bty = "l")
  par(op)
  dev.off()
  cat("已生成R文章图与区间坐标：", out_dir, "\n", sep = "")
}

# Rscript整文件入口：按本脚本位置找到项目，解压后可从任意工作目录运行。
if (sys.nframe() == 0L) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_file <- sub("^--file=", "", script_arg[1])
  project <- dirname(dirname(normalizePath(script_file)))
  draw_article(data_file = file.path(project, "data", "hetero_normal.csv"),
               out_dir = file.path(project, "文章配图"))
}
