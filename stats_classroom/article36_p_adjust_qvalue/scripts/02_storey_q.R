# 统计课堂36：同一组1000个P值，估计pi0并计算Storey q值。
# 本篇固定lambda=.5、pfdr=FALSE；qhat=pi0hat*BH是此算法的关系。
# 在项目根目录，从02.1开始；这份脚本无需先运行01。

# 02.1 准备输入和依赖 ----
if (!file.exists("data/screening_pvalues.csv")) {
  stop("请将工作目录设为项目根目录，再从02.1开始。")
}
if (!requireNamespace("qvalue",quietly=TRUE)) {
  stop("请按README安装Bioconductor的qvalue包。")
}
raw <- read.csv("data/screening_pvalues.csv")
head(raw)
p <- raw$p
m <- length(p)
lambda <- 0.5
alpha <- 0.05

# 02.2 手算真实零假设比例 ----
# 在均匀零假设模型中，P>=lambda的概率为1-lambda。
# 将右侧个数除以这段区间长度，可反推零假设个数，再除以m。
right_tail_count <- sum(p>=lambda)
pi0 <- min(1,right_tail_count/(m*(1-lambda)))
pi0_result <- data.frame(lambda=lambda,tail_count=right_tail_count,pi0=pi0)
pi0_result

# 02.3 手算本篇固定参数的q值 ----
order_p <- order(p)
sorted_p <- p[order_p]
rank_p <- seq_along(sorted_p)
bh_candidate <- m*sorted_p/rank_p
bh_sorted <- pmin(1,rev(cummin(rev(bh_candidate))))
bh <- numeric(m)
bh[order_p] <- bh_sorted
q_manual <- pi0*bh

# 02.4 调用qvalue包核对 ----
# qvalue计算pi0与qvalues；lambda明确固定，lfdr.out=FALSE不另估局部FDR。
fit_q <- qvalue::qvalue(p,lambda=lambda,pfdr=FALSE,lfdr.out=FALSE)
function_check <- data.frame(
  pi0_manual=pi0,
  pi0_function=fit_q$pi0,
  max_abs_q_error=max(abs(q_manual-fit_q$qvalues))
)
function_check
stopifnot(function_check$max_abs_q_error<1e-10)

# 02.5 查看q值与BH的发现数 ----
comparison <- raw
comparison$pi0 <- pi0
comparison$q <- fit_q$qvalues
comparison$BH <- bh
summary_result <- data.frame(
  method=c("q_Storey_lambda_0.5","BH"),
  discoveries=c(sum(comparison$q<=alpha),sum(bh<=alpha)),
  pi0=c(pi0,NA_real_)
)
summary_result

# 02.6 改变估计窗口：lambda敏感性选读 ----
lambda_grid <- c(.3,.5,.7)
pi0_grid <- numeric(length(lambda_grid))
for (i in seq_along(lambda_grid)) {
  current_lambda <- lambda_grid[i]
  pi0_grid[i] <- min(1,mean(p>=current_lambda)/(1-current_lambda))
}
lambda_sensitivity <- data.frame(lambda=lambda_grid,pi0=pi0_grid)
lambda_sensitivity
# 这一步比较估计设置的影响；正文采用预先指定的lambda=.5。

# 02.7 可选保存 ----
if (identical(Sys.getenv("ARTICLE36_SAVE_RESULTS"),"1")) {
  dir.create("results",showWarnings=FALSE)
  write.csv(comparison,"results/R_02_storey_q.csv",row.names=FALSE)
  write.csv(summary_result,"results/R_02_summary.csv",row.names=FALSE)
  write.csv(lambda_sensitivity,"results/R_02_lambda_sensitivity.csv",row.names=FALSE)
}
