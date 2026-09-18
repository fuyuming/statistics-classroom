# 2026-09-18｜第二次课·统计描述：标准正态与 u（Z）：跨量纲比较
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-12-数据标准化.R
# 运行：在RStudio打开本文件后Source，或命令行 Rscript 本文件。
# 依赖：只用到 base R 与 showtext（绘图显示中文）。缺包时执行：
#   install.packages("showtext", repos="https://cloud.r-project.org")
# 数据：data/newborn_2025.csv（10 名顺产新生儿的身高 cm 与体重 kg，本文件夹已有）
# 要点：Z = (x - μ) / σ 只改变位置与刻度，不改变分布形状；用来把不同量纲放到同一把尺子上。
# options(digits=3) 只影响显示的有效数字，不改变原数据；统计解释要保留单位和分母。

local({
  # 优先定位被 source 的文件，其次 Rscript 参数和 RStudio 当前编辑文件。
  paths <- Filter(function(x) is.character(x) && length(x) == 1 && nzchar(x),
                  lapply(sys.frames(), function(x) x$ofile))
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  path <- if (length(paths)) tail(paths, 1)[[1]] else if (length(args)) sub("^--file=", "", args[1]) else ""
  if (!nzchar(path) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(path)) {
    root <- dirname(normalizePath(path, mustWork = TRUE))
    if (basename(root) == "R" && file.exists(file.path(dirname(root), "医学统计学.Rproj"))) root <- dirname(root)
    setwd(root)
  }
  packages <- c("showtext")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("首次运行，正在安装缺失的包：", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  failed <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(failed)) stop("以下包尚未安装成功，请检查网络及安装日志：", paste(failed, collapse = ", "))
})

library(showtext)
showtext_auto()

options(digits = 3)
options(scipen = 200)
# === 结果保存开关：默认打开（结果直接存到本文件夹）===
# 只想到控制台看、不生成文件时，先运行 options(medstats.export = FALSE) 再跑本脚本。
if (is.null(getOption("medstats.export"))) options(medstats.export = TRUE)

# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)


options(medstats.export = isTRUE(getOption("medstats.export")))
# 画图与导出：RStudio 里直接显示；默认同时把图存成 output/R/<name>.pdf
show_fig <- function(draw, name, width = 9, height = 4.5) {
  if (interactive()) draw()
  if (isTRUE(getOption("medstats.export"))) {
    dir.create("output/R", showWarnings = FALSE, recursive = TRUE)
    pdf(file.path("output/R", paste0(name, ".pdf")), width = width, height = height)
    draw(); invisible(dev.off())
    cat(sprintf("已保存图形：output/R/%s.pdf\n", name))
  }
}

cat("\n===== 1. 先看问题：单位不同，能不能直接比大小？ =====\n")
newborn <- read.csv("./data/newborn_2025.csv")
cat("数据：", nrow(newborn), "名新生儿；变量：",
    paste(setdiff(names(newborn), "number"), collapse = "、"), "\n")
print(summary(newborn[, c("height_cm", "weight_kg")]))
cat("身高约 42.98~51.38 cm，体重约 3.13~3.85 kg。",
    "两者的数值大小、单位和波动都不同，直接比较“哪个偏离得更多”没有意义。\n")

cat("\n===== 2. Z 标准化：把两项指标放到同一把尺子上 =====\n")
Z <- scale(newborn[, c("height_cm", "weight_kg")])   # 按列：减均数、除标准差
cat(sprintf("标准化后：两列的均数为 %s，标准差为 %s\n",
            paste(sprintf("%.3f", colMeans(Z)), collapse = ", "),
            paste(sprintf("%.3f", apply(Z, 2, sd)), collapse = ", ")))
cat("也就是说：Z 告诉我们“这个观测值在均数之上/之下多少个标准差”。\n")
print(data.frame(id = newborn$number,
                 z_height = round(Z[, "height_cm"], 2),
                 z_weight = round(Z[, "weight_kg"], 2)))

cat("\n比较示例：\n")
# 同一名新生儿在两项指标上的相对位置，谁的 |Z| 更大，谁的偏离就更远。
i <- which.max(abs(Z[, "weight_kg"]))
dh <- Z[i, "height_cm"]; dw <- Z[i, "weight_kg"]
cat(sprintf("第 %d 名新生儿：身高 Z = %.2f，体重 Z = %.2f；", newborn$number[i], dh, dw),
    if (abs(dw) > abs(dh)) "相对本组，这名婴儿在体重上偏离得更远。\n"
    else "相对本组，这名婴儿在身高上偏离得更远。\n")

cat("\n===== 3. 标准化不改变分布形状 =====\n")
# 形状指标口径与 Python/MATLAB 共用模块一致：g1、g2 用中心矩（不校正），
# G1、G2 为小样本校正版；两者都可用，但报告时要写明用的是哪一套。
shape_moments <- function(x) {
  x <- as.numeric(x); n <- length(x); dx <- x - mean(x); m2 <- mean(dx^2)
  g1 <- mean(dx^3) / m2^1.5
  g2 <- mean(dx^4) / m2^2 - 3
  G1 <- g1 * sqrt(n * (n - 1)) / (n - 2)
  G2 <- ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))
  c(g1 = g1, G1 = G1, g2 = g2, G2 = G2)
}
for (v in c("height_cm", "weight_kg")) {
  a <- shape_moments(newborn[[v]]); b <- shape_moments(Z[, v])
  cat(sprintf("%s：原始 g1=%.3f、G2=%.3f ｜ 标准化后 g1=%.3f、G2=%.3f\n",
              v, a["g1"], a["G2"], b["g1"], b["G2"]))
  stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-12)))
}
cat("结论：偏态还是偏态，厚尾还是厚尾。Z 只是换了刻度，不会把右偏资料变成正态资料。\n")

cat("\n===== 4. 附：这几个数后面还会用到 =====\n")
# 变异系数 CV = SD / 均数 × 100%（比较相对波动）、几何均数与 GSD 等指标属下一讲内容，
# 这里只把数字列出来备查；本讲的重点是 Z 只换刻度、不改变分布形状。
for (v in c("height_cm", "weight_kg")) {
  cat(sprintf("%s：SD=%.3f，均数=%.2f（比值 %.2f%%＝CV，下一讲再讲）\n", v,
              sd(newborn[[v]]), mean(newborn[[v]]),
              100 * sd(newborn[[v]]) / mean(newborn[[v]])))
}

cat("\n===== 5. 图形对照 =====\n")
fig_compare <- function() {
  op <- par(mfrow = c(1, 2))
  boxplot(newborn[, c("height_cm", "weight_kg")], col = "lightblue",
          main = "原始尺度：单位不同，纵轴不能共用", ylab = "原始数值")
  boxplot(Z, col = "lightgreen",
          main = "Z 标准化后：纵轴同为“标准差个数”", ylab = "Z 值")
  abline(h = 0, col = "red", lty = 2)
  par(op)
}
show_fig(fig_compare, "2-12-数据标准化")

cat("\n提示：\n")
cat("1) Z 是有方向的：正号表示在均数之上，负号表示在均数之下；数值大小表示“多少个标准差”。\n")
cat("2) 只有当零点有意义、量纲可线性缩放时才谈 Z；等级资料先转成有序因子，不要直接标准化。\n")
cat("3) 标准化不改变形状，所以右偏资料该用中位数/几何均数时，标准化以后依然要这样描述。\n")

cat("\n脚本 2-12-数据标准化.R 运行完毕。\n")
