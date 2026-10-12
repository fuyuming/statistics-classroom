# 统计课堂08：精确分布、正态近似、率差与Q。只用基础R。
# 打开统计课堂09.Rproj，再Source完整文件；所有文本保存为UTF-8。
args <- commandArgs(FALSE)
file_arg <- sub('^--file=', '', args[grepl('^--file=',args)])
source_file <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
root <- getwd()
if(length(file_arg)) root <- dirname(dirname(normalizePath(file_arg[1])))
if(!is.null(source_file)) root <- dirname(dirname(normalizePath(source_file)))
out <- file.path(root,'运行结果','v2','R')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,'逐题结果.txt'),split=TRUE)
tryCatch({
save_csv <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,fileEncoding='UTF-8')

# 1. 单样本率：概率质量与正态小区间概率比较 ----
cases <- read.csv(file.path(root,'数据','binomial_cases_v2.csv'))
curves <- list()
for(i in seq_len(nrow(cases))) {
  n <- cases$n[i]; p <- cases$p[i]; x <- 0:n
  se <- sqrt(p*(1-p)/n)
  curves[[i]] <- data.frame(case=cases$case[i],n=n,p=p,x=x,proportion=x/n,
    exact=dbinom(x,n,p),normal_mass=pnorm((x+.5)/n,p,se)-pnorm((x-.5)/n,p,se))
}
save_csv(do.call(rbind,curves),'single_distribution.csv')
n <- 10; p <- .1; k <- 3
se <- sqrt(p*(1-p)/n)
# lower.tail=FALSE是P(X>q)，故“至少3”要传q=2。
metrics <- c(single_mean=p,single_var=se^2,single_se=se,
  single_exact=pbinom(k-1,n,p,lower.tail=FALSE),
  single_normal=pnorm(k/n,p,se,lower.tail=FALSE),
  single_cc=pnorm((k-.5)/n,p,se,lower.tail=FALSE))
cat('问题1：n=10、p=0.1，样本率至少0.3的概率。\n')
print(metrics)
cat('解释：计数边界从3校正到2.5，比例边界从0.3校正到0.25。np=1，校正仍有误差。\n')

# 2. 两组独立：枚举所有联合组合，合并同一个率差 ----
G <- read.csv(file.path(root,'数据','groups.csv'))
n1 <- G$n[1]; n2 <- G$n[2]; p1 <- G$p[1]; p2 <- G$p[2]
stopifnot(n1>0,n2>0,p1>0,p1<1,p2>0,p2<1)
mean_d <- p1-p2
var_d <- p1*(1-p1)/n1+p2*(1-p2)/n2
se_d <- sqrt(var_d)
# 最大公约数用于建立整数格点键，避免浮点数把相同率差拆开。
gcd <- function(a,b) {while(b!=0) {temp <- b; b <- a%%b; a <- temp}; a}
L <- n1*n2/gcd(n1,n2)
pairs <- expand.grid(x1=0:n1,x2=0:n2)
pairs$key <- pairs$x1*(L/n1)-pairs$x2*(L/n2)
pairs$difference <- pairs$key/L
pairs$probability <- dbinom(pairs$x1,n1,p1)*dbinom(pairs$x2,n2,p2)
save_csv(pairs,'joint_pairs.csv')
pmf <- aggregate(probability~key,pairs,sum)
pmf$difference <- pmf$key/L
save_csv(pmf,'difference_pmf.csv')
cut_key <- ceiling(.2*L-1e-12)
stopifnot(all(c(cut_key-1,cut_key,cut_key+1) %in% pmf$key))
cc_boundary <- (cut_key-.5)/L
exact_tail <- sum(pmf$probability[pmf$key>=cut_key])
metrics <- c(metrics,mean_d=mean_d,var_d=var_d,se_d=se_d,pair_count=nrow(pairs),
  difference_exact=exact_tail,difference_normal=pnorm(.2,mean_d,se_d,lower.tail=FALSE),
  difference_cc=pnorm(cc_boundary,mean_d,se_d,lower.tail=FALSE),cc_boundary=cc_boundary)
cat(sprintf('\n问题2：率差中心=%.9f，方差=%.9f，标准误=%.9f。\n',mean_d,var_d,se_d))
cat('解释：组间独立，协方差为0，差的方差才等于两份方差之和。\n')
cat(sprintf('问题3：%d种联合组合；精确尾概率=%.9f，未校正=%.9f，校正=%.9f。\n',
  nrow(pairs),exact_tail,metrics['difference_normal'],metrics['difference_cc']))
cat('给定总体率差0.1，计算样本率差至少0.2的概率；这不是检验两总体率相等的P值。\n')
stopifnot(abs(sum(pmf$probability)-1)<1e-12,
  abs(sum(pmf$difference*pmf$probability)-mean_d)<1e-12,
  abs(sum((pmf$difference-mean_d)^2*pmf$probability)-var_d)<1e-12)

# 3. 分类计数：负协方差、Q与三个标准化比较 ----
C <- read.csv(file.path(root,'数据','categories.csv'))
observed <- C$observed; p <- C$p; n <- sum(observed)
stopifnot(length(p)==4,all(abs(p-.25)<1e-12))
expected <- n*p
cov_counts <- n*(diag(p)-outer(p,p))
Q <- sum((observed-expected)^2/expected)
a <- rbind(c(1,-1,0,0),c(0,0,1,-1),c(1,1,-1,-1))
a <- sweep(a,1,sqrt(c(n/2,n/2,n)),'/')
z <- as.vector(a%*%observed)
stopifnot(max(abs(a%*%cov_counts%*%t(a)-diag(3)))<1e-12,
  abs(Q-sum(z^2))<1e-12)
metrics <- c(metrics,var_count=cov_counts[1,1],cov_count=cov_counts[1,2],Q=Q,
  z1_sq=z[1]^2,z2_sq=z[2]^2,z3_sq=z[3]^2,EQ=length(p)-1)
save_csv(as.data.frame(cov_counts),'count_covariance.csv')
cat('\n问题4：四分类计数协方差矩阵：\n'); print(cov_counts)
cat('计数总和方差=',sum(cov_counts),'；Q=',Q,'；三个平方=',z^2,'；E(Q)=3。\n')
cat('解释：三个方向只有渐近独立正态，Q的卡方参照是大样本近似。\n')
save_csv(data.frame(name=names(metrics),value=unname(metrics)),'theory.csv')
# 英文核对图避免中文字体依赖；中文公众号配图由Python入口生成。
png(file.path(out,'difference_check.png'),width=1200,height=800,res=140)
plot(pmf$difference,pmf$probability,type='h',xlim=c(-.25,.45),
     xlab='Difference of sample proportions',ylab='Probability',main='Exact and normal interval probabilities')
normal_mass <- pnorm(pmf$difference+.5/L,mean_d,se_d)-pnorm(pmf$difference-.5/L,mean_d,se_d)
lines(pmf$difference,normal_mass,col='#D67D38',lwd=2)
dev.off()
},finally=sink())
