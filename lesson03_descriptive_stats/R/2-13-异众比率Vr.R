# 2026-09-18｜第二次课·统计描述：分类数据的离散程度：异众比率 V_r
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-13-异众比率Vr.R
# 运行：在RStudio打开本文件后Source，或命令行 Rscript 2-13-异众比率Vr.R
# 说明：本脚本只用 base R（不依赖 tidyverse）。
# 公式：V_r ＝ 1 − 众类构成比 ＝ (Σf_i − f_m) / Σf_i，出自李春喜、姜丽娜、邵云、
#       张黛静、马建辉《生物统计学（第六版）》（科学出版社，2023）第二章 第三节 式 (2-22)。
#       颜艳、王彤《医学统计学》第5版未列这一指标（其第五章只有频数、构成比、率、相对比与标准化）。
# 数据：① 例 2-1 的 80 名患者诊室分布（第 10 页）② 同一批数据的性别分布

## ── 图形中文显示（可选）──────────────────────────────────────
# 与其他脚本一致：装了 showtext 就让图形里的中文正常显示；没装也不影响计算与结论。
# 用 try() 包住：字体文件缺失时只跳过显示设置，不会中断脚本。
if (requireNamespace("showtext", quietly = TRUE)) {
  try(showtext::showtext_auto(), silent = TRUE)
}

## ── 一、诊室分布 ────────────────────────────────────────────────
clinic <- c(A = 27, B = 17, C = 21, D = 15)
n <- sum(clinic)
percent <- clinic / n * 100
mode_class <- names(clinic)[which.max(clinic)]
vr <- 1 - max(clinic) / n

print(data.frame(人数 = clinic, 构成比 = round(percent, 2)))
cat(sprintf("n = %d｜众类 %s（构成比 %.2f%%）→ V_r = 1 − %.4f = %.4f = %.2f%%\n",
            n, mode_class, max(clinic) / n * 100, max(clinic) / n, vr, vr * 100))

## ── 二、四种情形的对照 ─────────────────────────────────────────
vrfun <- function(x) 1 - max(x) / sum(x)
contrast <- list(
  "四类完全均匀 20/20/20/20"   = c(20, 20, 20, 20),
  "全部集中在 A 80/0/0/0"     = c(80, 0, 0, 0),
  "众类同为 27 的近均匀 27/17/21/15" = c(27, 17, 21, 15),
  "众类同为 27 的极端不均 27/26/26/1" = c(27, 26, 26, 1),
  "例2-1 性别 女44/男36"       = c(44, 36)
)
res <- do.call(rbind, lapply(names(contrast), function(k) {
  x <- contrast[[k]]
  mode_label <- names(x)[which.max(x)]              # 未命名向量返回 NULL，下面统一成位置
  if (is.null(mode_label) || !nzchar(mode_label)) mode_label <- paste0("第", which.max(x), "类")
  data.frame(情形 = k, 众类 = mode_label,
             众类构成比 = round(max(x) / sum(x) * 100, 2),
             V_r = round(vrfun(x) * 100, 2), check.names = FALSE)
}))
print(res, row.names = FALSE)

# 四类完全均匀时 V_r 取到上限 1 − 1/K（K＝4 → 75%）；全在一类时 V_r ＝ 0。
stopifnot(abs(vr - 0.6625) < 1e-12,
          abs(vrfun(c(20, 20, 20, 20)) - 0.75) < 1e-12,
          abs(vrfun(c(80, 0, 0, 0))) < 1e-12,
          # 局限：V_r 只用到众类一个频数——这两个分布差得很远，V_r 却完全相同
          abs(vrfun(c(27, 17, 21, 15)) - vrfun(c(27, 26, 26, 1))) < 1e-12)

## ── 三、图形核对 ────────────────────────────────────────────────
draw_vr <- function() {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  barplot(clinic, ylab = "人数", main = sprintf("诊室分布 n=%d，V_r=%.2f%%", n, vr * 100))
  barplot(c(20, 20, 20, 20), ylab = "人数", main = "四类均匀：V_r=75%（上限）")
  par(op)
}
# 交互式（RStudio）里直接显示；Rscript 批处理时只保存 PNG，不生成 Rplots.pdf
if (interactive()) draw_vr()

# === 结果保存（默认开）：图与数字存到本文件夹 ===
# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)

# 只想到控制台看时，先运行 options(medstats.export = FALSE)。
if (isTRUE(getOption("medstats.export", TRUE))) {
  png("./output/R/2-13-异众比率-图.png", width = 1800, height = 1100, res = 150); draw_vr(); dev.off()
  write.csv(data.frame(诊室 = names(clinic), 人数 = as.integer(clinic),
                       构成比 = round(percent, 2)), "./output/R/2-13-异众比率-诊室分布.csv", row.names = FALSE)
  write.csv(res, "./output/R/2-13-异众比率-四种情形对照.csv", row.names = FALSE)
  cat("已保存：2-13-异众比率-图.png、2-13-异众比率-诊室分布.csv、2-13-异众比率-四种情形对照.csv\n")
}

cat("\n要点：分类数据也有离散程度——不在众类里的人越多，众数的代表性越差。\n")
cat("报告时同时给出各类频数与构成比（或条形图）；V_r 看不见众类之外的分布形态。\n")
