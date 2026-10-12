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

# 问题：农田土壤铅浓度的总体变异程度有多大？
# 李春喜等《生物统计学》第5版 p79 例4.22 的8个原始观测。
# 本文把原检验题改为95%区间估计。单位：微克/克（ug/g）。
# 条件：独立同分布正态样本；这项条件是本算例采用的模型假设。
x <- c(4.2, 4.5, 3.6, 4.7, 4.0, 3.8, 3.7, 4.2)
n <- length(x)
s_sq <- var(x)
dfree <- n - 1
q025 <- qchisq(0.025, dfree)
q975 <- qchisq(0.975, dfree)
# 大分位点作下限分母，小分位点作上限分母。
ci <- c(dfree*s_sq/q975, dfree*s_sq/q025)
print(c(n=n, s_sq=s_sq, df=dfree, q025=q025, q975=q975))
# 端点反代应恰好回到卡方分布的相应分位点。
stopifnot(abs(pchisq(dfree*s_sq/ci[1],dfree)-0.975) < 1e-10)
stopifnot(abs(pchisq(dfree*s_sq/ci[2],dfree)-0.025) < 1e-10)
# 标准差区间：对方差区间端点分别开平方，单位回到 ug/g。
print(c(SD_lower=sqrt(ci[1]), SD_upper=sqrt(ci[2])))
estimate <- s_sq
unit <- "(ug/g)^2"
cat("解释：这是总体方差区间，不是平均铅浓度区间，也不是95%土样的取值范围。\n")

# 保存并显示本例的区间结果
result <- data.frame(estimate=estimate, lower=ci[1], upper=ci[2], conf_level=0.95, unit=unit)
print(result)
write.csv(result, file.path(out, "04_一个总体方差.csv"), row.names=FALSE)
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
png(file.path(out, "04_一个总体方差.png"), width=1200, height=500, res=150)
plot_result()
dev.off()
cat("文件已保存至：", normalizePath(out), "\n")
