# 2026-09-18｜第二次课·统计描述：相对数三兄弟、应用相对数的注意事项、率的标准化
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-3-相对数与率的标准化.R
#   （相对数三兄弟 → 应用相对数的注意事项 → 率的标准化）
# 运行：在RStudio打开本文件后Source，或命令行 Rscript 2-3-相对数与率的标准化.R
# 说明：本脚本只用 base R（不依赖 tidyverse），无网环境也能直接跑；数字打印在控制台。
# 数据：① 讲义表 3.3 某地某年某肿瘤患病情况（226 例病例＋各年龄组人口数）
#       ② 教材例 5-3 出生性别比（506 名男婴、470 名女婴）
#       ③ 第一次课结石案例的分层数据＋共同标准构成 357/343（直接标准化）

## ── 图形中文显示（可选）──────────────────────────────────────
# 与其他脚本一致：装了 showtext 就让图形里的中文正常显示；没装也不影响计算与结论。
# 用 try() 包住：字体文件缺失时只跳过显示设置，不会中断脚本。
if (requireNamespace("showtext", quietly = TRUE)) {
  try(showtext::showtext_auto(), silent = TRUE)
}

## ── 一、构成比（结构相对数）与患病率（强度相对数）────────────────
ages  <- c("0~", "30~", "40~", "50~", "60~")
cases <- c(8, 21, 53, 84, 60)                          # 各年龄组患者数
pop   <- c(1012321, 506534, 574637, 592340, 201765)    # 各年龄组人口数

const_pct <- cases / sum(cases) * 100                  # 分母＝病例总数 226
prev_1e5  <- cases / pop * 1e5                         # 分母＝该年龄组人口数

tab <- data.frame(年龄组 = ages, 患者数 = cases, 患者构成比 = round(const_pct, 2),
                  人口数 = pop, 患病率每10万 = round(prev_1e5, 2), check.names = FALSE)
print(tab)
cat(sprintf("合计：构成比 %.2f%%｜患病率 %.2f/10万\n",
            sum(const_pct), sum(cases) / sum(pop) * 1e5))
# 只看构成比（第4列）会得出"50~ 组最容易患病"；换成患病率（第5列）才知道 60~ 组最高。
stopifnot(abs(const_pct[4] - 37.17) < 0.01, abs(prev_1e5[5] - 29.74) < 0.01)

## ── 二、相对比（教材例 5-3）────────────────────────────────────
sex_ratio <- 506 / 470 * 100                           # 男婴数 / 女婴数 × 100
cat(sprintf("出生性别比 = 506/470×100 = %.2f（每 100 名女婴约对应 108 名男婴）\n", sex_ratio))
stopifnot(abs(sex_ratio - 107.66) < 0.01)

## ── 三、率的标准化（接第一次课结石案例）────────────────────────
# 分层数据：每层"成功数 / 层内人数"
layer <- list(
  小结石 = list(A = c(81, 87),  B = c(234, 270)),
  大结石 = list(A = c(192, 263), B = c(55, 80))
)
std_pop <- c(小结石 = 357, 大结石 = 343)               # 共同标准构成（总标准人数 700）

rate <- lapply(layer, function(d) sapply(d, function(v) v[1] / v[2]))
crude <- sapply(c("A", "B"), function(g)
  sum(sapply(layer, function(d) d[[g]][1])) / sum(sapply(layer, function(d) d[[g]][2])))
stand <- sapply(c("A", "B"), function(g)
  sum(sapply(names(layer), function(s) rate[[s]][[g]] * std_pop[[s]])) / sum(std_pop))

cat("分层成功率：\n"); print(round(do.call(rbind, rate), 4))
cat(sprintf("粗率     A/B = %.4f / %.4f（%.2f%% / %.2f%%）\n",
            crude[["A"]], crude[["B"]], crude[["A"]] * 100, crude[["B"]] * 100))
cat(sprintf("标准化率 A/B = %.4f / %.4f（%.2f%% / %.2f%%），标准构成 %d/%d\n",
            stand[["A"]], stand[["B"]], stand[["A"]] * 100, stand[["B"]] * 100,
            std_pop[["小结石"]], std_pop[["大结石"]]))
# 每层 A 都更高，粗率却是 B 更高——差别来自两组的结石大小构成不同；统一权重后 A 更高。
stopifnot(abs(crude[["A"]] - 0.78) < 1e-9,
          abs(crude[["B"]] - 0.8257142857142857) < 1e-9,
          abs(stand[["A"]] - 0.8325) < 1e-4,
          abs(stand[["B"]] - 0.7789) < 1e-4)

## ── 四、图形核对（base graphics）──────────────────────────────
draw_rates <- function() {
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  barplot(const_pct, names.arg = ages, ylab = "患者构成比 (%)",
          main = "只看构成：50~ 组病例最多")
  barplot(prev_1e5, names.arg = ages, ylab = "患病率 (1/10万)",
          main = "换成率：60~ 组患病率最高")
  barplot(rbind(crude, stand), beside = TRUE, names.arg = c("A 疗法", "B 疗法"),
          ylab = "成功率", main = "粗率 vs 标准化率",
          legend.text = c("粗率", "标准化率"), args.legend = list(x = "topright", bty = "n"))
  par(op)
}
# 交互式（RStudio）里直接显示；Rscript 批处理时只保存 PNG，不生成 Rplots.pdf
if (interactive()) draw_rates()

# === 结果保存（默认开）：图与数字存到本文件夹 ===
# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)

# 只想到控制台看时，先运行 options(medstats.export = FALSE)。
if (isTRUE(getOption("medstats.export", TRUE))) {
  png("./output/R/2-3-相对数与率的标准化-图.png", width = 2200, height = 1000, res = 150); draw_rates(); dev.off()
  write.csv(data.frame(年龄组 = ages, 患者数 = cases, 患者构成比 = round(const_pct, 2),
                       人口数 = pop, 患病率每10万 = round(prev_1e5, 2)),
            "./output/R/2-3-构成比与患病率.csv", row.names = FALSE)
  write.csv(data.frame(疗法 = c("A", "B"), 粗率 = round(crude, 4), 标准化率 = round(stand, 4)),
            "./output/R/2-3-粗率与标准化率.csv", row.names = FALSE)
  cat("已保存：2-3-相对数与率的标准化-图.png、2-3-构成比与患病率.csv、2-3-粗率与标准化率.csv\n")
}

cat("\n要点：分子与分母决定它在回答什么问题；比例的含义取决于分母，构成不同时用标准化统一权重。\n")
cat("标准化只统一了已分层的因素（本例为结石大小），未测因素与研究设计限制仍需在解释中保留。\n")
