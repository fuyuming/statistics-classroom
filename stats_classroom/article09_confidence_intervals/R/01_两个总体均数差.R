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

# 问题：A、B 两批幼崽的平均出生重相差多少？单位 g。
# 马寨璞《基础生物统计学》(2015) p89，例9；原单位 kg，这里乘1000。
# 模型：两组独立，各为正态总体，总体方差相同且未知。
# 1. 输入完整原始数据
A <- c(143, 142, 143, 137)
B <- c(140, 142, 136, 138, 140)
# 2. 先计算两组均数、样本方差（var 的分母为 n-1）
n1 <- length(A)
n2 <- length(B)
m1 <- mean(A)
m2 <- mean(B)
s1_sq <- var(A)
s2_sq <- var(B)
difference <- m1 - m2
# 3. 合并方差、均数差标准误与 t 分位点
sp_sq <- ((n1 - 1)*s1_sq + (n2 - 1)*s2_sq)/(n1 + n2 - 2)
se <- sqrt(sp_sq*(1/n1 + 1/n2))
dfree <- n1 + n2 - 2
tcrit <- qt(0.975, dfree)
# 4. 下限与上限
ci <- difference + c(-1, 1)*tcrit*se
print(c(n1=n1, n2=n2, mean_A=m1, mean_B=m2, var_A=s1_sq, var_B=s2_sq,
        pooled_variance=sp_sq, SE=se, df=dfree, t975=tcrit))
# 5. 用内置函数核对。必须指定 var.equal=TRUE，默认是 Welch 方法。
check <- t.test(A, B, var.equal=TRUE, conf.level=0.95)$conf.int
stopifnot(max(abs(ci - as.numeric(check))) < 1e-10)
estimate <- difference
unit <- "g"
cat("解释：估计的是 A总体平均出生重减去B总体平均出生重。区间包含0，不能据此确定差异方向。\n")

# 保存并显示本例的区间结果
result <- data.frame(estimate=estimate, lower=ci[1], upper=ci[2], conf_level=0.95, unit=unit)
print(result)
write.csv(result, file.path(out, "01_两个总体均数差.csv"), row.names=FALSE)
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
png(file.path(out, "01_两个总体均数差.png"), width=1200, height=500, res=150)
plot_result()
dev.off()
cat("文件已保存至：", normalizePath(out), "\n")
