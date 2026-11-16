# 统计课堂28｜配对率差的四种区间（UTF-8）
# 数据：数据/paired_20.csv；依赖：ratesci、ExactCIdiff；运行：RStudio分节或Source；
# 输出：默认显示一幅图，SAVE_RESULTS=TRUE时写入运行结果/R。
# 完整多法入口：ratesci给Wald、MOVER-NW、AS(Tango)，ExactCIdiff另算Wang。

# ---- 准备：项目根或代码目录均可启动 ----
SAVE_RESULTS <- FALSE
args <- commandArgs(trailingOnly=FALSE)
script_arg <- grep("^--file=", args, value=TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else ""
script_dir <- if (nzchar(script_file) && file.exists(script_file)) dirname(normalizePath(script_file)) else getwd()
root_candidates <- c(getwd(), script_dir, dirname(script_dir))
root <- root_candidates[file.exists(file.path(root_candidates, "数据", "paired_20.csv"))][1]
stopifnot(file.exists(file.path(root, "数据", "paired_20.csv")))

# ---- 问题1：从共享数据计算样本率差 ----
records <- read.csv(file.path(root, "数据", "paired_20.csv"))
print(names(records))
print(head(records, 3))
changes <- records$after - records$before
n_plus <- sum(changes == 1)
n_zero <- sum(changes == 0)
n_minus <- sum(changes == -1)
n <- length(changes)
estimate <- (n_plus - n_minus) / n
cat(sprintf("(N+,T,N-)=(%d,%d,%d), n=%d, D=%.2f（后测减前测）\n", n_plus,n_zero,n_minus,n,estimate))

# ---- 问题2：四种区间如何比较？ ----
confidence <- 0.95
# ratesci的x按后测行×前测列组织；表1以前测为行，以下映射保持后测减前测。
# a、d是两次均阳性、均阴性，b是阴转阳，c是阳转阴。
a <- sum(records$before == 1 & records$after == 1)
d <- sum(records$before == 0 & records$after == 0)
x <- c(a=a, b=n_plus, c=n_minus, d=d)
# level是置信水平；std_est=TRUE统一报告样本D，cc=FALSE不加连续性校正。
rate_result <- ratesci::rdpairci(x=x, level=confidence, cc=FALSE, std_est=TRUE)
method_names <- c("Wald", "MOVER-NW", "AS(Tango)")
selected <- rate_result$estimates[method_names, ]
rate_table <- data.frame(method=method_names,
                         estimate=selected[,"est"],
                         lower=selected[,"lower"],
                         upper=selected[,"upper"], row.names=NULL)
# $estimates是多法表；PairedCI另算Wang。n12、t、n21依次为阴转阳、不变、阳转阴人数。
# conf.level指定置信水平；CItype="Two.sided"为双侧，结果$ExactCI为两端。
wang_unchanged <- x["a"] + x["d"]
wang <- ExactCIdiff::PairedCI(n12=x["b"], t=wang_unchanged,
                              n21=x["c"], conf.level=confidence,
                              CItype="Two.sided")
all_methods <- rbind(rate_table,
                     data.frame(method="Wang", estimate=unname(wang$estimate),
                                lower=wang$ExactCI[1], upper=wang$ExactCI[2]))
print(all_methods, row.names=FALSE)
cat(sprintf("confidence=%.2f；SE(Wald)=%.10f\n", confidence,
            sqrt(((x["b"]+x["c"])/sum(x)-estimate^2)/sum(x))))

# ---- 问题3：改变confidence时四种区间怎样变化？ ----
# 将问题2的confidence=.95改为.90后重跑；同一参数会同时更新四种方法。
# 横向图默认显示，不写盘；Rscript使用NULL设备，避免Rplots.pdf。
null_device <- !interactive()
if (null_device) pdf(file=NULL)
par(mar=c(4,8,3,1))
y <- rev(seq_len(nrow(all_methods)))
plot(all_methods$estimate, y, pch=19, col="#d95f02", xlim=range(c(all_methods$lower,all_methods$upper)),
     yaxt="n", xlab="Paired rate difference (after - before)",
     ylab="", main=sprintf("%.0f%% intervals for paired difference", 100*confidence))
axis(2, at=y, labels=all_methods$method, las=1)
abline(v=estimate, lty=2, col="gray50")
# i对应第i种方法；每条线连接该方法的下限和上限。
for (i in seq_len(nrow(all_methods))) {
  segments(all_methods$lower[i], y[i], all_methods$upper[i], y[i], lwd=3, col="#2f6f8f")
}
if (null_device) dev.off()

# 改参练习：将confidence=.95改为.90，R会同时重算四法且区间变窄。
if (SAVE_RESULTS) {
  output <- file.path(root,"运行结果","R")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  write.csv(all_methods,file.path(output,"四种配对率差区间.csv"),row.names=FALSE)
}
