# 精读03 v5：同一个作者模型，网格与 quap 两种计算。
# 打开本课 Rproj 后 Source；或 Rscript 代码/howell_v5.R。
# Howell1 的共享 CSV 已随项目提供；依赖 rethinking。
library(rethinking)
args_all <- commandArgs(FALSE)
file_arg <- grep('^--file=', args_all, value=TRUE)
root <- if(length(file_arg)) dirname(dirname(normalizePath(sub('^--file=', '', file_arg[1])))) else getwd()
out <- file.path(root,'运行结果/v5')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- read.csv(file.path(root,'运行结果/00_教材数据_Howell1成年人.csv'))
print(head(d)) # height 厘米，weight 千克；每行是一位成年人。
dat <- list(W=d$weight,H=d$height,Hbar=mean(d$height))
# 1. 作者 03_howell_new_weight_model.r 中 m_adults 的模型。
# dnorm 的第二个参数是标准差；dlnorm 的 0、1 是 log(b) 的均值和标准差。
fit <- quap(alist(
  W ~ dnorm(mu,sigma),
  mu <- a+b*(H-Hbar),
  a ~ dnorm(60,10),
  b ~ dlnorm(0,1),
  sigma ~ dunif(0,10)
),data=dat,start=list(a=45,b=.63,sigma=4.2))
print(precis(fit,prob=.89))
write.csv(data.frame(parameter=names(coef(fit)),mean=coef(fit),sd=sqrt(diag(vcov(fit)))),file.path(out,'R_quap_parameters.csv'),row.names=FALSE)
write.csv(vcov(fit),file.path(out,'R_quap_covariance.csv'),row.names=FALSE)
# 2. 网格：三参数都未知，不把 sigma 固定。
# 这些范围是计算窗口，不是新的先验。Python 入口另做加密和扩大窗口检查。
a <- seq(43,47,length.out=101)
b <- seq(.45,.8,length.out=101)
s <- seq(3.3,5.5,length.out=101)
g <- expand.grid(a=a,b=b,sigma=s)
x <- d$height-mean(d$height)
y <- d$weight
n <- length(y)
# 用汇总量计算残差平方和，与逐个累加 (y-mu)^2 完全相同。
sse <- sum(y*y)-2*g$a*sum(y)-2*g$b*sum(x*y)+n*g$a^2+2*g$a*g$b*sum(x)+g$b^2*sum(x*x)
logw <- -n*log(g$sigma)-sse/(2*g$sigma^2)+dnorm(g$a,60,10,log=TRUE)+dlnorm(g$b,0,1,log=TRUE)
# 相同网格体积与 Uniform(0,10) 的常数在归一化时约掉。
w <- exp(logw-max(logw))
w <- w/sum(w)
summary <- data.frame(parameter=c('a','b','sigma'),mean=c(sum(w*g$a),sum(w*g$b),sum(w*g$sigma)))
print(summary)
write.csv(summary,file.path(out,'R_grid_means.csv'),row.names=FALSE)
# 3. 从 quap 联合后验抽样，再模拟新个体（作者 link/sim 的同一思路）。
set.seed(20261005)
p <- extract.samples(fit,n=100000)
mu160 <- p$a+p$b*(160-dat$Hbar)
w160 <- rnorm(length(mu160),mu160,p$sigma)
intervals <- data.frame(quantity=c('mean_at_160','new_person_at_160'),mean=c(mean(mu160),mean(w160)),lower=c(quantile(mu160,.055),quantile(w160,.055)),upper=c(quantile(mu160,.945),quantile(w160,.945)))
write.csv(intervals,file.path(out,'R_quap_simulation_160.csv'),row.names=FALSE)
cat('斜率描述群体关联；新人区间比平均体重区间宽。随机模拟区间与正文数值可能略有差异。\n')
writeLines(c(R.version.string,paste('rethinking',packageVersion('rethinking'))),file.path(out,'R_versions.txt'))
