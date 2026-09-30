# 明哥的微生物世界 · 均值与方差的独立性 · R / RStudio
# 打开本课Rproj，再Source本脚本。仅需基础R，UTF-8。
# 输入均为共享人工/模型模拟数据，非真实实验。模拟不证明独立性。
# 输出：运行结果/R；英文标签配图：文章配图/R。

# 0. 项目与输出 ------------------------------------------------------------
stopifnot(dir.exists('数据'))
out <- file.path('运行结果', 'R')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
messages <- character()
explain <- function(x) {
  lines <- capture.output(print(x, row.names=FALSE))
  cat(paste(lines, collapse='\n'), '\n')
  messages <<- c(messages, lines)
}
save_table <- function(x, name) {
  write.csv(x, file.path(out, name), row.names=FALSE, fileEncoding='UTF-8')
}

# 1. 均值相同，方差可以不同吗？ --------------------------------------------
examples <- read.csv('数据/position_examples.csv')
explain('问题1：输入的每行包含一组人工教学数据。')
explain(examples)
values <- as.matrix(examples[c('x1', 'x2', 'x3')])
stopifnot(all(is.finite(values)))
examples$mean <- rowMeans(values)
# apply的第二个参数为1：按行；var默认分母n-1。
examples$variance <- apply(values, 1, var)
save_table(examples, 'position_results.csv')
explain(examples)
stopifnot(isTRUE(all.equal(examples$variance, c(4,4,16))))
explain('均值10、20、10；样本方差4、4、16。此例不证明独立性。')

# 2. Cov(M,D)的展开和n=10的重复抽样 -----------------------------------------
n <- 10
mu <- 10
sigma <- 2
# qt输入左侧累计概率；0.975对应双侧5%的上界。
t_limit <- qt(0.975, df=n-1)
summaries <- list()
conditions <- list()
identities <- list()
for (model in c('normal', 'shifted_exponential')) {
  data <- read.csv(file.path('数据', paste0(model, '_samples.csv')))
  explain(paste('问题2：', model, '共享输入前三行'))
  explain(head(data, 3))
  x <- as.matrix(data[paste0('x', 1:n)])
  stopifnot(all(is.finite(x)))
  means <- rowMeans(x)
  variances <- apply(x, 1, var)
  # R从1开始索引；x[,1]取每行的第一个观测。
  x1 <- x[,1]
  x2 <- x[,2]
  midpoint <- (x1+x2)/2
  difference <- x1-x2
  cov_md <- cov(midpoint, difference)
  rhs <- (var(x1)-var(x2))/2
  pair_variance <- apply(x[,1:2], 1, var)
  pair_error <- max(abs(pair_variance-difference^2/2))
  identities[[model]] <- data.frame(model=model, cov_MD=cov_md,
    half_var_difference=rhs, identity_error=abs(cov_md-rhs), max_pair_variance_error=pair_error)
  stopifnot(abs(cov_md-rhs)<1e-12, pair_error<1e-12)
  z <- (means-mu)/(sigma/sqrt(n))
  q <- (n-1)*variances/sigma^2
  t_constructed <- z/sqrt(q/(n-1))
  t_direct <- (means-mu)/sqrt(variances/n)
  stopifnot(max(abs(t_constructed-t_direct))<1e-12)
  save_table(data.frame(sample_id=data$sample_id, mean=means, variance=variances,
    M_n2=midpoint, D_n2=difference, S2_n2=pair_variance, Z=z, Q=q, T=t_direct),
    paste0(model, '_statistics.csv'))
  theoretical_cov <- if (model=='normal') 0 else 1.6
  rejection <- mean(abs(t_direct)>t_limit)
  summaries[[model]] <- data.frame(model=model, repetitions=nrow(x), n=n,
    mean_of_means=mean(means), mean_of_variances=mean(variances),
    cov_mean_variance=cov(means, variances), theoretical_cov=theoretical_cov,
    corr_mean_variance=cor(means, variances), t_rejection_fraction=rejection,
    rejection_mcse=sqrt(rejection*(1-rejection)/nrow(x)), t_975=t_limit)
  selections <- list(all=rep(TRUE,nrow(x)), mean_lt_9_5=means<9.5, mean_gt_10_5=means>10.5)
  for (label in names(selections)) {
    selected <- selections[[label]]
    v <- variances[selected]
    # type=7与NumPy默认线性分位数一致。
    quantiles <- quantile(v, c(.25,.5,.75), type=7, names=FALSE)
    conditions[[paste(model,label)]] <- data.frame(model=model, selection=label,
      count=sum(selected), mean_variance=mean(v), q25_variance=quantiles[1],
      median_variance=quantiles[2], q75_variance=quantiles[3])
  }
}
summary_table <- do.call(rbind,summaries)
conditional_table <- do.call(rbind,conditions)
save_table(do.call(rbind,identities), 'covariance_identity.csv')
save_table(summary_table, 'sampling_summary.csv')
save_table(conditional_table, 'conditional_summary.csv')
explain(summary_table)
explain(conditional_table)
explain('经验协方差不必等于理论值。条件分位数相近不能替代独立性的数学证明。')
explain('只有正态模型在此设定下保证精确t(9)分布和5%的理论拒绝率。')

# 3. 三种理论协方差为0的模型 -----------------------------------------------
pairs <- read.csv('数据/joint_pairs.csv')
explain('问题3：pair_id为模拟编号，sign为独立公平硬币决定的正负1。')
explain(head(pairs,3))
stopifnot(all(is.finite(as.matrix(pairs))), all(pairs$sign %in% c(-1,1)))
pairs$independent_y_calculated <- pairs$independent_y
pairs$square_y_calculated <- pairs$x^2
pairs$random_sign_y_calculated <- pairs$sign*pairs$x
pair_summary <- list()
for (label in c('independent','square','random_sign')) {
  y <- pairs[[paste0(label,'_y_calculated')]]
  pair_summary[[label]] <- data.frame(model=label, cov_XY=cov(pairs$x,y),
    corr_XY=cor(pairs$x,y), theoretical_cov=0)
}
save_table(pairs,'figure2_coordinates.csv')
save_table(do.call(rbind,pair_summary),'pair_summary.csv')
explain(do.call(rbind,pair_summary))
explain('后两种模型存在非线性依赖；第三种虽各自正态，却不联合正态。')

# 4. 用同一数据生成三张图；英文标签避免字体依赖 ------------------------------
fig_dir <- file.path('文章配图','R')
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
png(file.path(fig_dir,'01_position.png'),width=1260,height=1080,res=150)
par(mfrow=c(3,1),mar=c(3,4,3,1))
for (i in 1:3) {
  xs <- as.numeric(examples[i,c('x1','x2','x3')])
  plot(xs,rep(0,3),xlim=c(4,24),ylim=c(-.5,.7),yaxt='n',xlab='Value (dimensionless)',ylab='',
    pch=19,col='#167d80',cex=2,main=paste('Group',examples$group[i],
    ': mean =',examples$mean[i],'; sample variance =',examples$variance[i]))
  abline(v=examples$mean[i],lty=2,col='#cc773e')
  text(xs,.3,labels=xs)
}
dev.off()
png(file.path(fig_dir,'02_dependence.png'),width=1260,height=1950,res=150)
par(mfrow=c(3,1),mar=c(4,4,3,1))
labels <- c('Independent normals','Y = X squared: dependent','Y = random sign * X: dependent')
keys <- c('independent','square','random_sign')
for (i in 1:3) {
  limits <- if (i==2) c(-.3,10) else c(-3.4,3.4)
  plot(pairs$x,pairs[[paste0(keys[i],'_y_calculated')]],xlim=c(-3.4,3.4),ylim=limits,
    pch=16,cex=.5,col=adjustcolor('#167d80',alpha.f=.35),xlab='X',ylab='Y',main=labels[i])
}
dev.off()
source('代码/plot_parabola.R', encoding='UTF-8')
explain(R.version.string)
explain('练习：改变筛选均值的门槛，观察入选样本数和结果波动；不能用模拟证明独立。')
writeLines(messages,file.path(out,'运行记录.txt'),useBytes=TRUE)
