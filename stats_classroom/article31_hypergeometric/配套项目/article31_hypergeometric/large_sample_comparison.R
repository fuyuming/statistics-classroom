# 统计课堂31：两独立率的大样本近似与 Fisher 精确检验
# 题干：假设患者独立随机分组，各100人；4周用同一标准评价，达标数为新方案50、对照30。
# 目标是总体达标率差 p_new-p_control；双侧 H0: p_new=p_control，alpha=0.05。
SAVE_OUTPUT <- FALSE
args <- commandArgs(trailingOnly=FALSE)
fa <- args[grepl("^--file=", args)]
project_root <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
dat <- read.csv(file.path(project_root, "data", "large_sample_table.csv"), stringsAsFactors=FALSE)

# ---- 1. 读入数据并按组名取值 ----
new <- dat[dat$group == "new", ]
control <- dat[dat$group == "control", ]
x_new <- new$event
n_new <- new$total
x_control <- control$event
n_control <- control$total
p_new <- x_new/n_new
p_control <- x_control/n_control
RD <- p_new-p_control
observed_rates <- data.frame(p_new=p_new, p_control=p_control, rate_difference=RD)
observed_rates

# ---- 2. H0 合并率、期望表与 npq 检查 ----
# H0要求两组共享同一概率，用总达标数除总人数估计。
p_pool <- (x_new+x_control)/(n_new+n_control)
# 每行人数乘合并率，得到H0下达标/未达标的期望频数。
expected_success <- c(n_new,n_control)*p_pool
expected_failure <- c(n_new,n_control)*(1-p_pool)
# npq是二项计数方差的代入检查量，与单格期望人数不同。
npq <- c(n_new,n_control)*p_pool*(1-p_pool)
null_checks <- data.frame(group=c("new","control"), expected_success=expected_success, expected_failure=expected_failure, npq=npq)
null_checks

# ---- 3. 手算 Z、双侧正态 P ----
# 两独立样本率方差相加；H0下两组均用同一合并率。
SE0 <- sqrt(p_pool*(1-p_pool)*(1/n_new+1/n_control))
Z <- RD/SE0
# pnorm(-abs(Z))给左侧较远尾部；乘2得到双侧P。
P_normal <- 2*pnorm(-abs(Z))
manual_z <- data.frame(p_pool=p_pool, SE0=SE0, Z=Z, P_normal=P_normal, alpha=0.05)
manual_z

# ---- 4. 完整近似函数核对：prop.test 不加连续性校正 ----
prop_fit <- prop.test(c(x_new,x_control), c(n_new,n_control), alternative="two.sided", correct=FALSE)
prop_check <- data.frame(P_normal=P_normal, P_prop_test=prop_fit$p.value)
prop_check

# ---- 5. 固定边际下超几何分布与概率排序双侧 P ----
N <- n_new+n_control
K <- x_new+x_control
draw <- n_new
A <- 0:K
# m/n为总体达标/未达标数，k为新方案人数；这里求每张表的概率。
q <- dhyper(A, m=K, n=N-K, k=draw)
q_observed <- dhyper(x_new,m=K,n=N-K,k=draw)
# 1e-7仅是浮点相等容差，纳入概率不大于观察表的所有表。
P_hyper_two_sided <- sum(q[q <= q_observed*(1+1e-7)])
E_A <- sum(A*q)
Var_A <- sum((A-E_A)^2*q)
hypergeom_check <- data.frame(sum_probability=sum(q), E_A=E_A, Var_A=Var_A, P_ordered_two_sided=P_hyper_two_sided)
hypergeom_check

# ---- 6. Fisher 精确检验核对（只报告 P） ----
tab <- matrix(c(x_new,new$no_event,x_control,control$no_event), nrow=2, byrow=TRUE)
fisher_fit <- fisher.test(tab, alternative="two.sided", conf.int=FALSE)
fisher_check <- data.frame(P_hyper=P_hyper_two_sided, P_fisher=fisher_fit$p.value)
# P约.00593734；近似P约.00389242，两者都小于.05，支持两总体率不同。
fisher_check

# ---- 7. 可选保存结果（不绘图） ----
if (SAVE_OUTPUT) {
  out <- file.path(project_root,"results","r")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write.csv(data.frame(A=A,q=q),file.path(out,"large_sample_probability.csv"),row.names=FALSE)
  write.csv(data.frame(RD=RD,SE0=SE0,Z=Z,P_normal=P_normal,P_prop_test=prop_fit$p.value,P_hyper=P_hyper_two_sided,P_fisher=fisher_fit$p.value,E_A=E_A,Var_A=Var_A),file.path(out,"large_sample_comparison.csv"),row.names=FALSE)
}
