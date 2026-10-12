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

# 06 覆盖率模拟：重复抽样时，区间会变化，真值保持为2 g。
# 全部数据都是设定模型下的模拟，不是教材观测数据。
set.seed(20261012)
B <- 1000
n1 <- 10
n2 <- 10
mu1 <- 142
mu2 <- 140
sigma <- 2.5
truth <- mu1 - mu2
ans <- matrix(NA_real_, nrow=B, ncol=3)
colnames(ans) <- c("estimate", "lower", "upper")
for (i in seq_len(B)) {
  x1 <- rnorm(n1, mean=mu1, sd=sigma)
  x2 <- rnorm(n2, mean=mu2, sd=sigma)
  d <- mean(x1) - mean(x2)
  sp2 <- ((n1-1)*var(x1) + (n2-1)*var(x2))/(n1+n2-2)
  se <- sqrt(sp2*(1/n1 + 1/n2))
  ci <- d + c(-1,1)*qt(0.975, n1+n2-2)*se
  ans[i,] <- c(d, ci)
}
ans <- data.frame(round=seq_len(B), ans)
ans$covered <- ans$lower <= truth & ans$upper >= truth
print(c(rounds=B, coverage=mean(ans$covered), true_difference=truth))
write.csv(ans, file.path(out, "06_重复抽样.csv"), row.names=FALSE)
plot_coverage <- function() {
  show <- ans[1:40,]
  plot(NA, xlim=range(show$lower,show$upper), ylim=c(40.5,0.5),
       xlab="Population mean difference (g)", ylab="Round", main="95% CIs from 40 repeated samples")
  abline(v=truth, lty=2, col="#D65A4A")
  segments(show$lower, show$round, show$upper, show$round,
           col=ifelse(show$covered,"#1769AA","#D65A4A"), lwd=2)
}
if (interactive()) plot_coverage()
png(file.path(out,"06_重复抽样.png"),width=1000,height=1100,res=150)
plot_coverage()
dev.off()
cat("蓝线覆盖固定真值，红线未覆盖。每轮重新取样并重新计算区间。\n")
cat("表格和PNG保存在：",normalizePath(out),"\n")
