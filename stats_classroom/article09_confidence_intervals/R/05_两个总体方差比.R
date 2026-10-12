# 统计课堂：从抽样分布到置信区间
# 在 RStudio 中打开本项目 .Rproj，然后打开一个脚本，点 Source 完整运行。
# 每个脚本独立，不依赖其他脚本；只使用 R 自带包。
# Console 显示计算过程和结果，Plots 显示图；output 保存 CSV 与 PNG。
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (length(file_arg) > 0) {
  script_file <- sub("^--file=", "", file_arg[1])
  root <- dirname(dirname(normalizePath(script_file)))
} else if (!is.null(source_file)) {
  root <- dirname(dirname(normalizePath(source_file)))
} else {
  root <- getwd()  # 打开本项目 .Rproj 时，这就是项目目录
}
out <- file.path(root, "output")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
options(digits = 8)

# 问题：两个小麦品种的千粒重波动，相比有多大？方向固定为品种1/品种2。
# 李春喜等《生物统计学》第5版 p80 例4.23 的汇总值。
# 本文将原检验题改成区间估计；不编造原始逐行测量值。
# 模型：两组相互独立，各来自正态总体；无需假设两个总体方差相等。
# theta是未知总体方差比。方差齐性对应theta=1；一般是R/theta服从F。
# 区间估计theta的范围；若95%等尾区间不含1，对应双侧5%等尾F检验拒绝方差齐性。
n1 <- 10
n2 <- 10
s1_sq <- 22.933
s2_sq <- 2.933
ratio <- s1_sq/s2_sq
df1 <- n1 - 1
df2 <- n2 - 1
f025 <- qf(0.025, df1, df2)
f975 <- qf(0.975, df1, df2)
ci <- c(ratio/f975, ratio/f025)
print(c(ratio=ratio, df1=df1, df2=df2, f025=f025, f975=f975))
# 检查正向区间的端点反代概率。
stopifnot(abs(pf(ratio/ci[1],df1,df2)-0.975) < 1e-10)
stopifnot(abs(pf(ratio/ci[2],df1,df2)-0.025) < 1e-10)
# 顺序反过来：同时反转样本方差比与自由度。
reverse_ci <- (1/ratio)/qf(c(0.975,0.025),df2,df1)
stopifnot(max(abs(reverse_ci - 1/rev(ci))) < 1e-10)
estimate <- ratio
unit <- "variance ratio (unitless)"
cat("解释：这是方差之比；若要比较标准差之比，应把估计值与区间端点开平方。\n")

# 保存并显示本例的区间结果
result <- data.frame(estimate=estimate, lower=ci[1], upper=ci[2], conf_level=0.95, unit=unit)
print(result)
write.csv(result, file.path(out, "05_两个总体方差比.csv"), row.names=FALSE)
# 同一绘图函数既画到Plots，也保存PNG；英文轴标签便于跨系统显示。
plot_result <- function() {
  span <- ci[2]-ci[1]
  plot(estimate, 1, xlim=c(ci[1]-0.15*span,ci[2]+0.15*span), ylim=c(0.5,1.5),
       yaxt="n", ylab="", xlab=unit, main="95% confidence interval", pch=19, col="#1769AA")
  segments(ci[1],1,ci[2],1,lwd=3,col="#1769AA")
  points(estimate,1,pch=19,col="#D65A4A")
  text(ci, c(0.85,0.85), labels=formatC(ci,digits=4,format="f"))
}
if (interactive()) plot_result()
png(file.path(out, "05_两个总体方差比.png"), width=1200, height=500, res=150)
plot_result()
dev.off()
cat("文件已保存至：", normalizePath(out), "\n")
