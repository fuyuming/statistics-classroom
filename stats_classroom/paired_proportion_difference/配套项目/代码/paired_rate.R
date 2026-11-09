# 配对率差的三项模型：R主入口（UTF-8，base R）。
# 打开本篇Rproj后Source，结果分别保存在运行结果/R。
root <- if (basename(getwd())=="代码") normalizePath("..") else getwd()
data_dir <- file.path(root,"数据")
out_dir <- file.path(root,"运行结果","R")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

# 一、逐人相加得到精确分布，kernel依次对应-1、0、+1。
distribution <- function(n,qp,qm) {
  q0 <- 1-qp-qm
  probability <- 1
  for (i in seq_len(n)) {
    # 直接按三种变化移动概率，避免相关/卷积函数对顺序的不同约定。
    next_probability <- numeric(length(probability)+2)
    index <- seq_along(probability)
    next_probability[index] <- next_probability[index]+qm*probability
    next_probability[index+1] <- next_probability[index+1]+q0*probability
    next_probability[index+2] <- next_probability[index+2]+qp*probability
    probability <- next_probability
  }
  k <- -n:n
  d <- k/n
  theta <- qp-qm
  sd <- sqrt((qp+qm-theta^2)/n)
  edges <- (d[-length(d)]+d[-1])/2
  # 相邻取值中点作为边界，两个端点格合并正态外尾。
  normal <- diff(c(0,pnorm(edges,theta,sd),1))
  data.frame(k=k,D=d,exact_probability=probability,normal_cell_probability=normal)
}
spec <- data.frame(scenario=c("normal_n2","normal_n20","normal_n1000","rare_n1000"),n=c(2,20,1000,1000),qp=c(.05,.05,.05,.001),qm=c(.25,.25,.25,.002))
moments <- list()
for (i in seq_len(nrow(spec))) {
  setting <- spec[i,]
  result <- distribution(setting$n,setting$qp,setting$qm)
  write.csv(result,file.path(out_dir,paste0("distribution_",setting$scenario,".csv")),row.names=FALSE)
  theta <- setting$qp-setting$qm
  variance <- setting$qp+setting$qm-theta^2
  error <- max(abs(cumsum(result$exact_probability)-cumsum(result$normal_cell_probability)))
  moments[[i]] <- data.frame(scenario=setting$scenario,n=setting$n,qp=setting$qp,q0=1-setting$qp-setting$qm,qm=setting$qm,theta=theta,mean_D=theta,var_D=variance/setting$n,max_abs_cdf_error=error)
}
write.csv(do.call(rbind,moments),file.path(out_dir,"moments.csv"),row.names=FALSE)

# 二、读个人记录，计算率差、插件标准误与Wald近似区间。
people <- read.csv(file.path(data_dir,"paired_20.csv"))
changes <- people$after-people$before
stopifnot(all(people$before %in% c(0,1)),all(people$after %in% c(0,1)),all(changes==people$z))
scenarios <- read.csv(file.path(data_dir,"scenarios.csv"))
examples <- list()
for (i in seq_len(nrow(scenarios))) {
  setting <- scenarios[i,]
  if (setting$scenario=="main_20") {
    n <- nrow(people)
    plus <- sum(changes==1)
    zero <- sum(changes==0)
    minus <- sum(changes== -1)
    both_positive <- sum(people$before==1 & people$after==1)
  } else {
    n <- setting$n
    plus <- setting$nplus
    zero <- setting$nzero
    minus <- setting$nminus
    both_positive <- setting$both_positive
  }
  difference <- (plus-minus)/n
  se <- sqrt(((plus+minus)/n-difference^2)/n)
  before <- (both_positive+minus)/n
  after <- (both_positive+plus)/n
  h <- 2*asin(sqrt(after))-2*asin(sqrt(before))
  interval <- difference+c(-1,1)*qnorm(.975)*se
  examples[[i]] <- data.frame(scenario=setting$scenario,n=n,nplus=plus,nzero=zero,nminus=minus,D=difference,SE_Wald=se,CI_low_Wald=interval[1],CI_high_Wald=interval[2],before_rate=before,after_rate=after,Cohen_h_marginal=h)
}
result <- do.call(rbind,examples)
print(result)
write.csv(result,file.path(out_dir,"paired_examples.csv"),row.names=FALSE)

# 三、用R重现正文三张统计图 ----
# 分布与尾部试算均由R计算，PNG与PDF保存在文章配图目录。
source(file.path(root, "代码", "plot_figures.R"))

# 四、逐组列出20人三类人数的231种结果及其概率 ----
source(file.path(root, "代码", "three_category_probabilities.R"))
