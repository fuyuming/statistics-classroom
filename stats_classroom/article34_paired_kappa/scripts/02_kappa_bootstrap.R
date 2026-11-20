# 统计课堂34：配对κ的Bootstrap、百分位区间与BCa区间 ----
# 背景：第11课例4的58名疑似系统性红斑狼疮患者，各接受免疫荧光法A与乳胶凝集法B检测抗核抗体。
# 题源：颜艳、王彤《医学统计学》第5版（2020）例7-3、表7-3，教材pp110–111；四格为[11,12;2,33]。
# 目标：在给定四格经验分布下重抽样，量化κ的条件波动，并比较百分位/BCa区间。
# 计数解释：bootstrap_counts.csv的每行a,b,c,d是一次联合重抽样的四格计数，合计58，
# 不是虚构的58名病历；R节仍用boot::boot真实生成并用boot::boot.ci取区间。
# 按RStudio整文件或从第一节顺序运行；工作目录为项目根目录。

# 02.1 问题与输入 ----
# 每个患者有一对完整结果；n是58名患者，B是50000次Bootstrap重抽样，seed固定以复现共享计数。
source("helpers/kappa_bootstrap_R.R")
raw <- read.csv("data/paired_tests.csv", stringsAsFactors=FALSE)
tab <- matrix(raw$count, nrow=2, byrow=TRUE,
              dimnames=list(c("positive","negative"), c("positive","negative")))
counts <- as.integer(raw$count); n <- sum(counts); B <- 50000L; seed <- 20261120L
stopifnot(identical(counts, c(11L,12L,2L,33L)), n == 58L)

# 02.2 手动计算：原样本κ与联合重抽样对象 ----
# 四格计数顺序是A阳/B阳、A阳/B阴、A阴/B阳、A阴/B阴；先按定义展开原表κ。
a <- counts[1]; b <- counts[2]; cc <- counts[3]; d <- counts[4]
kappa_hat <- 2*(a*d-b*cc)/((a+b)*(b+d)+(cc+d)*(a+cc))
cell <- rep(seq_len(4), times=counts) # 每个元素代表一种联合结果类别
point <- data.frame(a=a, b=b, c=cc, d=d, n=n, kappa=kappa_hat)
point
# κ=2(ad−bc)/[(a+b)(b+d)+(c+d)(a+c)]；整数计数让同值比较不受浮点改写影响。

# 02.3 手动计算：delete-one jackknife、z0、acceleration与分位水平 ----
# 共享计数每行是[a,b,c,d]；先显式按整数式计算κ*，再用58个delete-one样本计算L。
shared <- as.matrix(read.csv("data/bootstrap_counts.csv"))
stopifnot(nrow(shared)==B, ncol(shared)==4L, all(rowSums(shared)==n))
sa <- shared[,1]; sb <- shared[,2]; sc <- shared[,3]; sd <- shared[,4]
shared_kappa <- 2*(sa*sd-sb*sc)/((sa+sb)*(sb+sd)+(sc+sd)*(sa+sc))
if (any(!is.finite(shared_kappa))) stop("本例出现不可定义样本，不能静默删除。")
jack <- vapply(seq_len(n), function(i) kappa_from_counts(tabulate(cell[-i], nbins=4L)), numeric(1))
if (any(!is.finite(jack))) stop("jackknife出现不可定义样本。")
L <- (n-1) * (mean(jack)-jack)
z0_manual <- qnorm(mean(shared_kappa < kappa_hat)) # 严格小于，不用midrank
acc <- sum(L^3)/(6*sum(L^2)^1.5)
zalpha <- qnorm(c(.025,.975))
alpha_bca <- pnorm(z0_manual+(z0_manual+zalpha)/(1-acc*(z0_manual+zalpha)))
alpha_bca

# 02.4 调用函数：boot::boot真实生成50000组计数 ----
# statistic每次接收58个有放回索引，返回κ及四格计数；boot::boot的$t保存B行结果。
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
statistic <- function(data, indices) {
  v <- tabulate(data[indices], nbins=4L)
  c(kappa=kappa_from_counts(v), a=v[1], b=v[2], c=v[3], d=v[4])
}
boot_fit <- boot::boot(data=cell, statistic=statistic, R=B, sim="ordinary", stype="i")
stopifnot(nrow(shared)==B, ncol(shared)==4L, all(rowSums(shared)==n),
          all(unname(boot_fit$t[,2:5]) == unname(shared)))
if (any(!is.finite(boot_fit$t[,1]))) stop("本例出现不可定义样本，不能静默删除。")
z0 <- qnorm(mean(boot_fit$t[,1] < kappa_hat))
stopifnot(abs(z0-z0_manual) < 1e-15)

# 02.5 调用函数：boot::boot.ci的percentile与BCa ----
# percentile只读取排序后的Bootstrap端点；BCa还使用L和z0，均不是“严格精确”区间。
ci <- boot::boot.ci(boot_fit, conf=.95, type=c("perc","bca"), index=1, L=L)
ci
percentile <- ci$percent[1,4:5]
bca <- ci$bca[1,4:5]
manual_percentile <- boot:::norm.inter(boot_fit$t[,1], c(.025,.975))[,2]
manual_percentile
stopifnot(max(abs(percentile-manual_percentile)) < 1e-14)

# 02.6 结果核对：列含义与预期值 ----
result <- bootstrap_summary_row(boot_fit$t, kappa_hat, jack, B, seed,
  percentile, bca, alpha_bca, shared)
result
if (identical(Sys.getenv("ARTICLE34_SAVE_RESULTS"), "1")) {
  dir.create("results", showWarnings=FALSE)
  write.csv(result, "results/R_bootstrap.csv", row.names=FALSE)
}
# bootstrap_mean/bias是给定当前样本的条件重抽样均值/偏差估计，不证明总体κhat严格无偏。
# BCa只调整分位点，不向四格加入伪计数；本脚本不承诺自动选择“最佳”区间。
