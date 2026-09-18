# 2026-09-18｜第二次课·统计描述：正态分布（μ 与 σ）
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-11-正态分布.R
# 运行：在RStudio打开本文件后Source，或命令行 Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("showtext"), repos="https://cloud.r-project.org")
# 结果默认保存到 output/R/；只想到控制台看时先运行 options(medstats.export = FALSE)。
# 说明：本次课讲正态曲线的形状（μ 与 σ）、从直方图到光滑曲线，以及 d/p/q/r 四种读法；
#       覆盖比例（68-95-99.7）、偏态资料（MIC）与正态性检验留到下一讲。
local({
  # 优先定位被source的文件，其次Rscript参数和RStudio当前编辑文件。
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
})

# === 结果保存开关：默认打开（结果直接存到 output/R/）===
# 只想到控制台看、不生成文件时，先运行 options(medstats.export = FALSE) 再跑本脚本。
if (is.null(getOption("medstats.export"))) options(medstats.export = TRUE)

# 防止中文出现乱码
if (requireNamespace("showtext", quietly = TRUE)) try(showtext::showtext_auto(), silent = TRUE)

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

cat("\n===== 1. 正态曲线的两个参数：μ 决定位置，σ 决定胖瘦=====\n")
# 左：σ 相同、改变 μ → 曲线整体平移；右：μ 相同、改变 σ → σ 越大越矮胖。
# 这正是密度函数里 μ 只出现在 (x−μ)、σ 既在分母 √(2π) 前又在指数上的直接结果。
fig_params <- function() {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  x <- seq(-5, 12.5, length.out = 1000)
  plot(x, dnorm(x, 0, 1), type = "l", lwd = 2, col = "steelblue", ylim = c(0, 0.45),
       main = "σ 相同，改变 μ：曲线整体平移", xlab = "取值", ylab = "概率密度")
  lines(x, dnorm(x, 3, 1), lwd = 2, col = "darkorange")
  lines(x, dnorm(x, 6, 1), lwd = 2, col = "forestgreen")
  legend("topright", legend = c("μ=0, σ=1", "μ=3, σ=1", "μ=6, σ=1"),
         col = c("steelblue", "darkorange", "forestgreen"), lwd = 2, bty = "n")
  plot(x, dnorm(x, 0, 0.6), type = "l", lwd = 2, col = "steelblue", ylim = c(0, 0.7),
       main = "μ 相同，改变 σ：σ 越大越矮胖", xlab = "取值", ylab = "概率密度")
  lines(x, dnorm(x, 0, 1.0), lwd = 2, col = "darkorange")
  lines(x, dnorm(x, 0, 2.0), lwd = 2, col = "forestgreen")
  legend("topright", legend = c("μ=0, σ=0.6", "μ=0, σ=1.0", "μ=0, σ=2.0"),
         col = c("steelblue", "darkorange", "forestgreen"), lwd = 2, bty = "n")
  par(op)
}
show_fig(fig_params, "2-11-1-正态曲线的参数", width = 10, height = 4.5)
cat(sprintf("峰高＝1/(σ√(2π))：σ=0.6 时 %.3f，σ=1 时 %.3f，σ=2 时 %.3f；曲线在 x＝μ±σ 各有拐点。\n",
            1 / (0.6 * sqrt(2 * pi)), 1 / (1 * sqrt(2 * pi)), 1 / (2 * sqrt(2 * pi))))
# 密度函数（教材式 2-17）：f(x) = 1/(σ√(2π)) · exp[−(x−μ)²/(2σ²)]

cat("\n===== 2. 从直方图到光滑曲线：正态分布的由来=====\n")
# 同一批数据（例2-3，138 名成年女子红细胞数）：组段越细，直方图顶端越接近光滑钟形曲线。
# 纵坐标换成频率密度（频率 ÷ 组距）后，每个直方的面积＝该组频率，总面积＝1。
rbc0 <- read.csv("./data/rbc.csv")$rbc
cat(sprintf("红细胞：n=%d，min=%.2f，max=%.2f，均数=%.3f，SD=%.3f\n",
            length(rbc0), min(rbc0), max(rbc0), mean(rbc0), sd(rbc0)))
fig_hist <- function() {
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  hist(rbc0, breaks = seq(3.0, 5.7, by = 0.3), probability = TRUE, col = "lightblue",
       border = "white", ylim = c(0, 1.0), main = "a. 粗分组（组距 0.3）",
       xlab = "红细胞计数（×10¹²/L）", ylab = "频率密度")
  hist(rbc0, breaks = seq(3.0, 5.7, by = 0.1), probability = TRUE, col = "lightblue",
       border = "white", ylim = c(0, 1.0), main = "b. 细分组（组距 0.1）",
       xlab = "红细胞计数（×10¹²/L）", ylab = "频率密度")
  hist(rbc0, breaks = seq(3.0, 5.7, by = 0.1), probability = TRUE, col = "lightblue",
       border = "white", ylim = c(0, 1.0), main = "c. 组段不断分细 → 光滑钟形曲线",
       xlab = "红细胞计数（×10¹²/L）", ylab = "频率密度")
  grid_x <- seq(3.0, 5.7, length.out = 400)
  lines(grid_x, dnorm(grid_x, mean(rbc0), sd(rbc0)), col = "red", lwd = 2)
  par(op)
}
show_fig(fig_hist, "2-11-2-从直方图到光滑曲线", width = 12, height = 4)
cat("读图：组距从 0.3 减到 0.1，直方顶端逐渐连成两头低、中间高、左右对称的钟形——这就是正态分布的形状。\n")
cat("颜艳、王彤《医学统计学》第5版 第一章也指出：随机误差一般呈正态分布。\n")

cat("\n===== 3. 用课堂数据看正态：红细胞（近正态）=====\n")
rbc <- read.csv("./data/rbc.csv")$rbc
n <- length(rbc); m <- mean(rbc); s <- sd(rbc)
cat(sprintf("红细胞：n=%d，均数=%.3f，标准差=%.3f，中位数=%.3f\n", n, m, s, median(rbc)))
fig_rbc <- function() {
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  hist(rbc, breaks = seq(3.0, 5.7, by = 0.2), probability = TRUE, col = "lightblue",
       border = "white", main = "138 名成年女子红细胞数（×10¹²/L）",
       xlab = "红细胞计数", ylab = "频率密度")
  curve(dnorm(x, mean = m, sd = s), add = TRUE, col = "red", lwd = 2)
  qqnorm(rbc, main = "红细胞数的正态 Q-Q 图", pch = 19, cex = 0.7, col = "steelblue")
  qqline(rbc, col = "red", lwd = 2)
  par(op)
}
show_fig(fig_rbc, "2-11-3-课堂数据看正态")
cat("读图：直方图近似对称、Q-Q 图上的点大致沿直线分布——可以当作近似正态的资料描述。\n")

cat("\n===== 4. 分布函数的四种读法：d / p / q / r=====\n")
# dnorm 给密度高度（本身不是概率）；pnorm 给累积概率（曲线下面积）；
# qnorm 由概率反查分位点；rnorm 生成服从 N(μ, σ²) 的随机数。
cat(sprintf("dnorm(0)      = %.4f（标准正态在 0 处的密度高度）\n", dnorm(0)))
cat(sprintf("pnorm(1.96)   = %.4f（曲线下面积，即 P(Z ≤ 1.96)）\n", pnorm(1.96)))
cat(sprintf("qnorm(0.975)  = %.4f（由概率反查分位点）\n", qnorm(0.975)))
set.seed(20260918)
cat(sprintf("rnorm(5,0,1)  = %s\n", paste(sprintf("%.3f", rnorm(5)), collapse = ", ")))
# 教材例 2-15：血红蛋白近似 N(145.63, 12.05²)
mu15 <- 145.63; sd15 <- 12.05
cat(sprintf("例2-15  P(X<120)=%.4f；P(120≤X≤160)=%.4f；P(X>160)=%.4f（附表查得 0.0166／0.8664／0.1170）\n",
            pnorm(120, mu15, sd15),
            pnorm(160, mu15, sd15) - pnorm(120, mu15, sd15),
            1 - pnorm(160, mu15, sd15)))
cat("说明：μ±1σ、±2σ、±3σ 各覆盖多少（68-95-99.7）留到下一讲展开，本讲只用到“区间面积＝概率”这一层。\n")

cat("\n===== 5.【备查，本次课不讲】正态性检验 =====\n")
# 本次课只要求会看直方图与 Q-Q 图。检验法与 P 值留到后续课程，届时注意：
#   1) P 值大不能证明“服从正态”，只能说当前证据不足以拒绝；
#   2) 均数与标准差若由本样本估计，普通 K-S 检验的 P 值不成立，需要 Lilliefors 校正；
#   3) n 很小时检验功效低，n 很大时轻微偏离也会显著，判读始终要回到图形与专业知识。
if (interactive()) {
  print(shapiro.test(rbc))
} else {
  cat("（非交互运行，跳过检验演示；在 RStudio 中运行本段可看到结果）\n")
}

cat("\n脚本 2-11-正态分布.R 运行完毕。\n")
