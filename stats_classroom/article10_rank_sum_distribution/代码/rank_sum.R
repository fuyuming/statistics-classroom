# 统计课堂10：秩和分布。基础R即可；UTF-8保存。
# 打开统计课堂10.Rproj后Source完整文件。输入为人工教学评分。
args <- commandArgs(FALSE)
farg <- sub('^--file=', '', args[grepl('^--file=',args)])
src <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
root <- getwd()
if (length(farg)) root <- dirname(dirname(normalizePath(farg[1])))
if (!is.null(src)) root <- dirname(dirname(normalizePath(src)))
out <- file.path(root,'运行结果','R')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sink(file.path(out,'逐题结果.txt'),split=TRUE)
tryCatch({
data <- read.csv(file.path(root,'数据','scores.csv'))
cat('人工教学输入：group组别、value无并列评分、value_tied并列评分。\n')
print(data)
is_a <- data$group=='A'
m <- sum(is_a)
n <- sum(!is_a)
N <- m+n
stopifnot(m==3,n==3,length(unique(data$value))==N)
# 1. 选观测记录的下标，每列是A组的一种分配；并列时也不删除记录。
allocations <- combn(seq_len(N),m)
mean_w <- m*(N+1)/2
results <- list()
for (column in c('value','value_tied')) {
    ranks <- rank(data[[column]],ties.method='average')
    w_obs <- sum(ranks[is_a])
    u_obs <- w_obs-m*(m+1)/2
    all_w <- apply(allocations,2,function(index) sum(ranks[index]))
    # 完整等可能分布的方差用概率平均；不能用var(all_w)的19分母。
    var_exact <- mean((all_w-mean(all_w))^2)
    p_two <- mean(abs(all_w-mean_w)>=abs(w_obs-mean_w)-1e-12)
    tied_sizes <- as.numeric(table(data[[column]]))
    correction <- sum(tied_sizes^3-tied_sizes)
    var_formula <- m*n/12*(N+1-correction/(N*(N-1)))
    cat('\n问题：',column,'；合并秩：',ranks,'\n')
    cat(sprintf('W=%.6f；U=%.6f；均值=%.6f；方差=%.6f；双侧P=%.6f。\n',w_obs,u_obs,mean(all_w),var_exact,p_two))
    cat('解释：P按离零假设中心至少同样远定义；20种记录分配各占1/20。\n')
    allocation_table <- data.frame(t(allocations),W=all_w)
    names(allocation_table) <- c('record1','record2','record3','W')
    write.csv(allocation_table,file.path(out,paste0(column,'_allocations.csv')),row.names=FALSE)
    freq <- table(all_w)
    pmf <- data.frame(W=as.numeric(names(freq)),count=as.numeric(freq),probability=as.numeric(freq)/length(all_w))
    print(pmf)
    write.csv(pmf,file.path(out,paste0(column,'_pmf.csv')),row.names=FALSE)
    stopifnot(abs(var_exact-var_formula)<1e-12,abs(mean(all_w)-mean_w)<1e-12)
    results[[column]] <- c(W=w_obs,U=u_obs,variance=var_exact,p=p_two)
}
# 2. 注意R输出名称W，实际减掉了第一组最小秩和。
test <- wilcox.test(data$value[is_a],data$value[!is_a],exact=TRUE)
cat('\nR两独立样本无并列精确检验：\n')
print(test)
stopifnot(abs(unname(test$statistic)-results$value['U'])<1e-12,
          abs(test$p.value-results$value['p'])<1e-12)
cat('此处软件W=1，对应本文U；本文原始秩和W=7。\n')
p_normal <- 2*pnorm((results$value['W']+.5-mean_w)/sqrt(results$value['variance']))
cat('本例连续性校正双侧正态近似=',p_normal,'；精确值=0.2。\n')
name <- c('W','U','mean_W','var_W','p_exact','var_tied','p_tied','p_normal_cc')
value <- c(results$value['W'],results$value['U'],mean_w,results$value['variance'],results$value['p'],results$value_tied['variance'],results$value_tied['p'],p_normal)
write.csv(data.frame(name=name,value=as.numeric(value)),file.path(out,'theory.csv'),row.names=FALSE)

# 3. 大样本模拟：不放回抽名次，不要求模拟原始数值。
set.seed(20261011)
m_large <- 20
n_large <- 20
N_large <- m_large+n_large
repeats <- 20000
w_sim <- replicate(repeats,sum(sample.int(N_large,m_large,replace=FALSE)))
center <- m_large*(N_large+1)/2
sd_w <- sqrt(m_large*n_large*(N_large+1)/12)
z_sim <- (w_sim-center)/sd_w
write.csv(data.frame(W=w_sim,Z=z_sim),file.path(out,'large_simulation.csv'),row.names=FALSE)
cat(sprintf('\n模拟：W理论中心=%.6f；SD=%.6f；Z模拟均值=%.6f、方差=%.6f。\n',center,sd_w,mean(z_sim),var(z_sim)))
cat('模拟不是证明；同样随机种子不保证三语言抽样序列相同。\n')
# 英文标签便于跨系统显示；公众号中文图由Python绘图文件生成。
png(file.path(out,'normal_check.png'),width=1200,height=800,res=140)
hist(z_sim,breaks=seq(-8,8,.2),probability=TRUE,xlim=c(-4.5,4.5),col='#a6cdc6',border='white',main='Rank-sum null simulation: m=n=20',xlab='Standardized W')
curve(dnorm(x),from=-4.5,to=4.5,add=TRUE,col='#5363aa',lwd=2)
dev.off()
},finally=sink())
