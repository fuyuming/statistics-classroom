# 统计课堂23：逐组检查分布，并用同一张图展示样本与正态QQ图。
# 函数只负责计算和绘图；不会自动读文件、运行函数或改变工作目录。

check_distribution <- function(dat, output_dir = "运行结果/R",
                               figure_dir = "文章配图") {
  required <- c("id", "group", "days")
  if (!is.data.frame(dat) || !all(required %in% names(dat))) {
    stop("dat必须是含id、group、days三列的数据框。")
  }
  if (anyNA(dat[, required])) {
    stop("id、group、days不能含缺失值。")
  }
  if (!is.numeric(dat$days) || any(!is.finite(dat$days))) {
    stop("days必须是有限数值。")
  }
  labels <- c("A", "B", "C")
  dat$group <- as.character(dat$group)
  if (!all(dat$group %in% labels)) {
    stop("group只能包含A、B、C。")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # R函数接受每组3到5000个观测；这里也给出清楚的错误提示。
  n_by_group <- vapply(labels, function(g) sum(dat$group == g), integer(1))
  if (any(n_by_group < 3L) || any(n_by_group > 5000L)) {
    stop("每组样本量必须在3到5000之间，才能进行Shapiro-Wilk检验。")
  }
  normality <- do.call(rbind, lapply(labels, function(g) {
    x <- dat$days[dat$group == g]
    sw <- shapiro.test(x)
    data.frame(group = g, n = length(x), W = unname(sw$statistic),
               P = sw$p.value, check.names = FALSE)
  }))
  row.names(normality) <- NULL
  write.csv(normality, file.path(output_dir, "normality.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # qqnorm(plot.it=FALSE)返回的坐标原样保存，便于课堂复核图上的点。
  qq_parts <- lapply(labels, function(g) {
    x <- dat$days[dat$group == g]
    qq <- qqnorm(x, plot.it = FALSE)
    qx <- as.numeric(qq$x)
    qy <- as.numeric(qq$y)
    # qqline默认以理论正态分位数qnorm(.25/.75)和样本四分位数定线。
    q_theoretical <- qnorm(c(0.25, 0.75))
    q_sample <- as.numeric(quantile(x, c(0.25, 0.75), names = FALSE, type = 7))
    q1x <- q_theoretical[1]
    q3x <- q_theoretical[2]
    q1y <- q_sample[1]
    q3y <- q_sample[2]
    slope <- (q3y - q1y) / (q3x - q1x)
    intercept <- q1y - slope * q1x
    stopifnot(isTRUE(all.equal(c(q1x, q3x), qnorm(c(0.25, 0.75)))))
    list(group = g,
         coordinates = data.frame(group = g, theoretical = qx,
                                   sample = qy),
         line = data.frame(group = g, intercept = intercept, slope = slope,
                           q1_theoretical = q1x, q3_theoretical = q3x))
  })
  qq_coordinates <- do.call(rbind, lapply(qq_parts, `[[`, "coordinates"))
  qq_lines <- do.call(rbind, lapply(qq_parts, `[[`, "line"))
  row.names(qq_coordinates) <- NULL
  row.names(qq_lines) <- NULL
  write.csv(qq_coordinates, file.path(output_dir, "qq_coordinates.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(qq_lines, file.path(output_dir, "qq_lines.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # 不用随机抖动：同一天数按出现次序向上排开，因此每次运行坐标完全一致。
  colours <- c(A = "#2F6DB0", B = "#E68A2E", C = "#3A9D5D")
  qq_x_all <- qq_coordinates$theoretical
  qq_x_lim <- max(2.2, max(abs(qq_x_all), na.rm = TRUE) + 0.25)
  png(file.path(figure_dir, "04-各组分布与QQ图.png"), width = 1000,
      height = 2400, res = 120, type = "cairo-png", family = "Noto Sans CJK SC")
  oldpar <- par(no.readonly = TRUE)
  on.exit({par(oldpar); dev.off()}, add = TRUE)
  par(mfrow = c(4, 1), mar = c(4.8, 5.7, 3.1, 1.3),
      oma = c(0, 0, 0, 0), family = "Noto Sans CJK SC", cex = 1.55,
      mgp = c(2.6, 0.8, 0), tcl = -0.3)

  # 顶部合并图共用0~22天横轴；同一天数只在所属组的行附近确定性错开。
  group_row <- setNames(c(3, 2, 1), labels)
  plot(NA, xlim = c(0, 22), ylim = c(0.4, 3.6), xaxt = "n", yaxt = "n",
       xlab = "存活天数", ylab = "组别", main = "各组样本分布")
  axis(1, at = seq(0, 22, by = 5), cex.axis = 1.05)
  axis(2, at = unname(group_row), labels = labels, las = 1, cex.axis = 1.1)
  abline(h = unname(group_row), col = "grey90", lty = 1)
  for (g in labels) {
    x <- dat$days[dat$group == g]
    occurrence <- ave(seq_along(x), x, FUN = seq_along)
    # 最大偏移约0.15，既能看出重复值，又不会混入相邻组。
    offsets <- (occurrence - (ave(occurrence, x, FUN = max) + 1) / 2) * 0.075
    y <- group_row[g] + offsets
    points(x, y, pch = 21, bg = colours[g], col = "white", cex = 1.35)
    med <- median(x)
    segments(med, group_row[g] - 0.17, med, group_row[g] + 0.17,
             col = colours[g], lwd = 3)
    text(med, group_row[g] + 0.24, sprintf("中位数=%g", med),
         col = colours[g], cex = 0.85, pos = 2)
  }

  for (g in labels) {
    part <- qq_parts[[match(g, labels)]]
    coords <- part$coordinates
    line <- part$line
    plot(coords$theoretical, coords$sample, xlim = c(-qq_x_lim, qq_x_lim),
         ylim = c(0, 22), xlab = "正态理论分位数", ylab = "存活天数",
         main = sprintf("%s组  正态QQ图", g), col = colours[g], pch = 21,
         bg = colours[g], cex = 1.25)
    abline(a = line$intercept, b = line$slope, col = colours[g],
           lty = 2, lwd = 1.5)
  }
  invisible(list(normality = normality, qq_coordinates = qq_coordinates,
                 qq_lines = qq_lines))
}
