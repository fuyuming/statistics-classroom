# 统计课堂09：期望、方差、协方差。只用基础R，无需安装包。
# 打开统计课堂09.Rproj，Source完整脚本；UTF-8保存。
# 输入是课堂教学例子；随机数据是模型模拟，不是真实实验。
args <- commandArgs(FALSE)
file_arg <- sub('^--file=', '', args[grepl('^--file=', args)])
source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
root <- getwd()
if (length(file_arg)) root <- dirname(dirname(normalizePath(file_arg[1])))
if (!is.null(source_file)) root <- dirname(dirname(normalizePath(source_file)))
out <- file.path(root, '运行结果', 'R')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sink(file.path(out,'逐题结果.txt'), split=TRUE)
tryCatch({
groups <- read.csv(file.path(root,'数据','groups.csv'))
categories <- read.csv(file.path(root,'数据','categories.csv'))
cat('共享输入：n为样本量，p为给定总体概率，observed为教学计数。\n')
print(groups)
print(categories)
stopifnot(all(groups$n>0), all(groups$p>0 & groups$p<1),
          nrow(categories)==4, all(abs(categories$p-.25)<1e-12))

# 1. 两独立样本率之差：先求各自方差，再相加 ----
n1 <- groups$n[1]
n2 <- groups$n[2]
p1 <- groups$p[1]
p2 <- groups$p[2]
mean_d <- p1-p2
var1 <- p1*(1-p1)/n1
var2 <- p2*(1-p2)/n2
var_d <- var1+var2
se_d <- sqrt(var_d)
cat(sprintf('\n问题1：率差期望=%.6f；方差=%.6f；标准误=%.6f。\n',mean_d,var_d,se_d))
cat('解释：两组独立，因此协方差为0。标准误为9.19个百分点。\n')

# 2. 总数固定的分类计数：不同类别有负协方差 ----
observed <- categories$observed
p <- categories$p
n <- sum(observed)
expected <- n*p
# diag(p)把概率放到对角线；outer(p,p)给出所有p_j*p_l。
cov_counts <- n*(diag(p)-outer(p,p))
cat('\n问题2：各类计数协方差表，对角线为方差。\n')
print(cov_counts)
cat('所有元素相加=计数之和的方差：',sum(cov_counts),'\n')
contribution <- (observed-expected)^2/expected
Q <- sum(contribution)
z1 <- (observed[1]-observed[2])/sqrt(n/2)
z2 <- (observed[3]-observed[4])/sqrt(n/2)
z3 <- (observed[1]+observed[2]-observed[3]-observed[4])/sqrt(n)
z <- c(z1,z2,z3)
cat('\n问题3：每类对Q的贡献：',contribution,'；Q=',Q,'\n')
cat('三个比较的平方：',z^2,'；合计=',sum(z^2),'\n')
cat('解释：Q=4是一次数据的值；给定且正确的概率模型下E(Q)=3。\n')
# 可选矩阵核对：三个比较方差均为1，互相协方差为0。
a <- rbind(c(1,-1,0,0),c(0,0,1,-1),c(1,1,-1,-1))
a <- sweep(a,1,sqrt(c(n/2,n/2,n)),'/')
cov_z <- a %*% cov_counts %*% t(a)
stopifnot(max(abs(cov_z-diag(3)))<1e-12, abs(Q-sum(z^2))<1e-12,
          abs(sum(cov_counts))<1e-12)
metrics <- data.frame(name=c('mean_d','var_d','se_d','var_count','cov_count','Q','z1_sq','z2_sq','z3_sq','EQ'),
                      value=c(mean_d,var_d,se_d,cov_counts[1,1],cov_counts[1,2],Q,z^2,length(p)-1))
write.csv(metrics,file.path(out,'theory.csv'),row.names=FALSE)

# 3. 重复抽样，观察统计量的波动 ----
set.seed(20261011)
repeats <- 20000
c1 <- rbinom(repeats,size=n1,prob=p1)
c2 <- rbinom(repeats,size=n2,prob=p2)
differences <- c1/n1-c2/n2
counts <- t(rmultinom(repeats,size=n,prob=p)) # 每行一轮，四列为四类。
centered <- sweep(counts,2,expected,'-')
scaled_squares <- sweep(centered^2,2,expected,'/')
q_sim <- rowSums(scaled_squares)
write.csv(data.frame(difference=differences),file.path(out,'rate_simulation.csv'),row.names=FALSE)
write.csv(data.frame(counts,Q=q_sim),file.path(out,'category_simulation.csv'),row.names=FALSE)
cat(sprintf('\n模拟：率差均值=%.6f；方差=%.6f；Q均值=%.6f。\n',mean(differences),var(differences),mean(q_sim)))
cat('var使用n-1分母估计模拟分布方差；模拟不是证明，三语言随机结果不必逐行相同。\n')
cat('练习：四类计数改为10、20、25、25，Q应为7.5。概率保持各0.25。\n')
# 英文标签避免不同系统中文字体缺失；公众号中文图由Python绘图脚本生成。
png(file.path(out,'Q_check.png'),width=1200,height=800,res=140)
hist(q_sim,breaks=seq(0,max(q_sim)+.5,.5),probability=TRUE,
     col='#a7cdc7',border='white',xlim=c(0,20),main='Pearson Q: multinomial simulation',xlab='Q')
curve(dchisq(x,df=3),from=0,to=20,add=TRUE,col='#5264a7',lwd=2)
dev.off()
}, finally=sink())
