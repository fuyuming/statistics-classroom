# 正文跟读：按文章顺序保留核心计算，RStudio Source即可。
# 完整批注、中文图与导出结果仍见lesson05_v2.R或tidyverse版。
source_path <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if(length(arg)) source_path <- sub("^--file=", "", arg[1])
root <- if(!is.null(source_path)) dirname(dirname(normalizePath(source_path))) else getwd()
old_dir <- getwd()
setwd(root)
# 片段 1 --------------------
n_x1 <- c(10, 90)   # 两个Z组内，X=1的总人数
n_x0 <- c(90, 10)   # 两个Z组内，X=0的总人数
y1_x1 <- c(1, 81)   # 其中Y=1的人数
y1_x0 <- c(9, 9)
c(全体X1 = sum(y1_x1) / sum(n_x1),
  全体X0 = sum(y1_x0) / sum(n_x0)) # 0.82与0.18
cbind(Z = 0:1, X1组 = y1_x1/n_x1, X0组 = y1_x0/n_x0)

# 片段 2 --------------------
d <- read.csv("数据/v2/WaffleDivorce.csv")
head(d, 3)  # 先核对原始数值
zscore <- function(x) as.numeric(scale(x))
dat <- list(D = zscore(d$Divorce),
            M = zscore(d$Marriage),
            A = zscore(d$MedianAgeMarriage))

# 片段 3 --------------------
library(rethinking)
fit <- quap(
  alist(D ~ dnorm(mu, sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0, 0.2),
        bM ~ dnorm(0, 0.5), bA ~ dnorm(0, 0.5),
        sigma ~ dexp(1)),
  data = dat,
  start = list(a=0, bM=0, bA=-0.5, sigma=0.8)
)
precis(fit, prob = 0.89) # 查看均值、SD和89%区间，重点看bM、bA

# 片段 4 --------------------
set.seed(20261009)
post <- extract.samples(fit, n = 10000) # 从quap近似后验抽样
post <- post[post$sigma > 0, ]          # sigma必须为正
age_same <- sample(dat$A, nrow(post), replace = TRUE)
mu0 <- with(post, a + bM * 0 + bA * age_same)
mu1 <- with(post, a + bM * 1 + bA * age_same)
delta_mean <- mu1 - mu0
all.equal(delta_mean, post$bM)         # TRUE：核对斜率含义
quantile(delta_mean, c(0.055, 0.5, 0.945))

# 片段 5 --------------------
D0 <- rnorm(nrow(post), mean = mu0, sd = post$sigma)
D1 <- rnorm(nrow(post), mean = mu1, sd = post$sigma)
delta_observed <- D1 - D0
c(均值差的后验SD = sd(delta_mean),
  两次观测之差的SD = sd(delta_observed))

# 片段 6 --------------------
p_fungus <- c(不处理 = 0.5, 处理 = 0.1)
mean_growth <- (1-p_fungus)*5 + p_fungus*2
mean_growth                       # 3.5、4.7
unname(mean_growth[2]-mean_growth[1]) # 总效应1.2
c(无真菌时 = 5-5, 有真菌时 = 2-2)     # 两组内的差都为0

# 片段 7 --------------------
# 四项依次为：低低、低高、高低、高高（新颖性在前）
all_n <- c(100, 100, 100, 100)
funded_n <- all_n * c(0.2, 0.9, 0.9, 0.9)
trust_rate <- function(n) {
  c(新颖性低 = n[2]/sum(n[1:2]),
    新颖性高 = n[4]/sum(n[3:4]))
}
rbind(全部申请 = trust_rate(all_n),
      获资助者 = trust_rate(funded_n))

# 片段 8 --------------------
set.seed(20261009)
n_sim <- 100000
X <- rbinom(n_sim, 1, 0.5)
Z <- rbinom(n_sim, 1, 0.1+0.8*X)
Y <- rbinom(n_sim, 1, 0.1+0.8*Z)
A <- rbinom(n_sim, 1, 0.1+0.8*Z)
c(全体 = cor(X,Y),
  固定Z = cor(X[Z==0],Y[Z==0]),
  固定A = cor(X[A==0],Y[A==0]))

setwd(old_dir)
