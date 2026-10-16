# 问题：先算100尾的功效，再计算至少需要多少尾。
# 教材数值另拟右侧设计；只用R基础包。打开项目Rproj，Source本文件。
args <- commandArgs(trailingOnly=FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(file_arg)>0) dirname(dirname(normalizePath(file_arg[1]))) else getwd()
out <- file.path(root, "运行结果", "R")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
parameters <- read.csv(file.path(root, "数据", "设计参数.csv"))
print(parameters) # 教学设计参数，不是原始观测数据；体长单位cm
mu0 <- parameters$mu0
mu1 <- parameters$mu1
sigma <- parameters$sigma
alpha <- parameters$alpha
target <- parameters$target_power
delta <- mu1-mu0
d <- delta/sigma
# 1. qnorm输入左侧累计概率，返回标准正态分位数。
n <- 100
se <- sigma/sqrt(n)
cutoff <- mu0+qnorm(1-alpha)*se
# pnorm计算指定真实均数下未达到拒绝门槛的概率。
beta <- pnorm((cutoff-mu1)/se)
power <- 1-beta
# 2. ceiling向上取整；检查前一个整数尚未达到目标。
n_continuous <- ((qnorm(1-alpha)+qnorm(target))/d)^2
n_min <- ceiling(n_continuous)
power_previous <- pnorm(qnorm(1-alpha)-d*sqrt(n_min-1), lower.tail=FALSE)
power_minimum <- pnorm(qnorm(1-alpha)-d*sqrt(n_min), lower.tail=FALSE)
result <- data.frame(quantity=c("se","cutoff","beta","power","d","n_continuous","n_min","power_previous","power_minimum"), value=c(se,cutoff,beta,power,d,n_continuous,n_min,power_previous,power_minimum))
print(result)
write.csv(result,file.path(out,"结果.csv"),row.names=FALSE)
stopifnot(power_previous < target, power_minimum >= target)
cat("功效是指定真实差异下的检出概率，不是H1为真的概率。\n")
# 3. 以下改为双侧观察结果，三种情境的差值均为0.40 cm。
observed_n <- c(25,100,400)
observed_se <- sigma/sqrt(observed_n)
p_two_sided <- 2*pnorm(abs(0.40/observed_se),lower.tail=FALSE)
observed <- data.frame(n=observed_n,difference=0.40,p_two_sided=p_two_sided,ci_low=0.40-qnorm(.975)*observed_se,ci_high=0.40+qnorm(.975)*observed_se)
print(observed)
write.csv(observed,file.path(out,"观察结果情境.csv"),row.names=FALSE)
