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

# 问题：100例患者中45例在28天达到教材规定的客观缓解标准，估计总体缓解率。
# 颜艳、王彤《医学统计学》第5版，p92-93，例6-3。
# 模型：各患者结局独立，采用同一评价标准，具有相同事件概率p。
# 1. 输入事件数和总人数；这里是汇总数据，不编造患者逐行记录。
x <- 45
n <- 100
phat <- x/n
# p未知，用样本率代入检查np与n(1-p)均>5这一经验条件。
expected_est <- c(n*phat, n*(1-phat))
print(data.frame(item=c("n*p_hat", "n*(1-p_hat)"), count=expected_est, greater_than_5=expected_est>5))
# 这只是近似的经验筛查，不是覆盖率保证。
# 2. Wald 大样本近似区间，与P38的公式一致
se <- sqrt(phat*(1 - phat)/n)
zcrit <- qnorm(0.975)
ci <- phat + c(-1, 1)*zcrit*se
print(c(events=x, nonevents=n-x, phat=phat, SE=se, z975=zcrit))
# 3. 两种其他方法作为对照，不能期待不同方法数值完全相同。
# 单样本 prop.test(correct=FALSE) 给出 Wilson 得分区间，并非 Wald。
wilson <- as.numeric(prop.test(x, n, correct=FALSE)$conf.int)
exact <- as.numeric(binom.test(x, n)$conf.int)
comparison <- data.frame(method=c("Wald", "Wilson", "Clopper-Pearson"),
                         lower=c(ci[1],wilson[1],exact[1]), upper=c(ci[2],wilson[2],exact[2]))
print(comparison)
write.csv(comparison, file.path(out, "02_三种方法对照.csv"), row.names=FALSE)
# Wilson 公式与 R 函数对照
center <- (phat + zcrit^2/(2*n))/(1 + zcrit^2/n)
half <- zcrit*sqrt(phat*(1-phat)/n + zcrit^2/(4*n^2))/(1 + zcrit^2/n)
stopifnot(max(abs(wilson - (center + c(-1,1)*half))) < 1e-8)
estimate <- phat
unit <- "proportion (0 to 1)"
cat("解释：两个端点乘100后是百分比；本例Wald区间是正态近似结果。\n")

# 保存并显示本例的区间结果
result <- data.frame(estimate=estimate, lower=ci[1], upper=ci[2], conf_level=0.95, unit=unit)
print(result)
write.csv(result, file.path(out, "02_一个总体率.csv"), row.names=FALSE)
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
png(file.path(out, "02_一个总体率.png"), width=1200, height=500, res=150)
plot_result()
dev.off()
cat("文件已保存至：", normalizePath(out), "\n")
