# 统计课堂：从抽样分布到置信区间
# 在 RStudio 中打开本项目 .Rproj，然后打开一个脚本，点 Source 完整运行。
# 各编号脚本可独立运行；需要保留项目目录及共享CSV。只使用R自带包。
# Console 显示计算过程和结果，Plots 显示图；运行结果/R 保存 CSV 与 PNG。
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
out <- file.path(root, "运行结果", "R")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
options(digits = 8)

# 问题：两块麦田的总体锈病发病率相差多少？方向为低洼地减高坡地。
# 李春喜等《生物统计学》第5版 p71-72 例4.14，p78 例4.21。
# 使用教材的原始计数；原例用99%，本文统一为95%，并用未四舍五入数据重算。
# 教学模型：两地独立取样，每地株间独立且发病概率相同。
data <- read.csv(file.path(root, "数据", "wheat_rust_counts.csv"))
print(data)
x1 <- data$events[data$group == "lowland"]
n1 <- data$n[data$group == "lowland"]
x2 <- data$events[data$group == "upland"]
n2 <- data$n[data$group == "upland"]
p1 <- x1/n1
p2 <- x2/n2
# 总体率未知，分别以两组样本率检查四项预期次数。
expected_est <- c(n1*p1, n1*(1-p1), n2*p2, n2*(1-p2))
print(data.frame(item=c("n1*p1_hat","n1*(1-p1_hat)","n2*p2_hat","n2*(1-p2_hat)"),
                 count=expected_est, at_least_5=expected_est>=5))
difference <- p1 - p2
# 分别估计两组方差；不把两个样本率合并。
se <- sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
zcrit <- qnorm(0.975)
ci <- difference + c(-1,1)*zcrit*se
print(c(p1=p1, p2=p2, difference=difference, SE=se))
print(c(events1=x1, nonevents1=n1-x1, events2=x2, nonevents2=n2-x2))
# 两样本、双侧、不校正时，prop.test 的率差区间与这里的公式核对。
check <- prop.test(c(x1,x2), c(n1,n2), correct=FALSE)$conf.int
stopifnot(max(abs(ci - as.numeric(check))) < 1e-10)
# 如需99%区间，重新取分位点，并保持原始计数计算。
print(c(CI99_lower=difference-qnorm(0.995)*se, CI99_upper=difference+qnorm(0.995)*se))
estimate <- difference
unit <- "proportion difference"
cat("解释：差值及区间乘100，单位为百分点，不是相对增加的百分比。\n")

# 保存并显示本例的区间结果
result <- data.frame(estimate=estimate, lower=ci[1], upper=ci[2], conf_level=0.95, unit=unit)
print(result)
write.csv(result, file.path(out, "03_两个总体率差.csv"), row.names=FALSE)
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
png(file.path(out, "03_两个总体率差.png"), width=1200, height=500, res=150)
plot_result()
dev.off()
cat("文件已保存至：", normalizePath(out), "\n")
