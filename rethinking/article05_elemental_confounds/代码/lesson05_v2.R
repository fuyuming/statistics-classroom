# 精读05 v2：RStudio打开精读05_v2.Rproj，再Source本文件。
# 作者PPT31页的quap模型；输入已附，无须先跑其他语言。
# 二元结构精确枚举，区别于PPT中1000次随机模拟的频数。
# 第一部分对应正文二至四：真实婚姻数据的贝叶斯回归。
# 第二部分对应四种因果结构：已知生成规则下的概率计算，不需要拟合后验。
library(rethinking)
# Source时使用脚本位置；项目根目录中运行也可。
source_path <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
args <- commandArgs(FALSE)
arg <- grep("^--file=", args, value=TRUE)
if (length(arg)) source_path <- sub("^--file=", "", arg[1])
root <- if (!is.null(source_path)) dirname(dirname(normalizePath(source_path))) else getwd()
out <- file.path(root,"运行结果/v2/R")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(file.path(root,"数据/v2/WaffleDivorce.csv"))
# 先读数据：每行一个州或特区；Divorce、Marriage单位为每千名成年人，
# MedianAgeMarriage是结婚年龄中位数（岁）。本项目CSV列名应保持一致。
# 换数据时先核对列名、单位与缺失值，不能把个体记录直接套入本例解释。
print(head(d))
# standardize减均值除样本标准差，与PPT一致。
dat <- list(D=standardize(d$Divorce), M=standardize(d$Marriage), A=standardize(d$MedianAgeMarriage))
# D：离婚率，M：结婚率，A：结婚年龄，三者均已标准化。
# 似然 D ~ Normal(mu,sigma)：描述同一条件均值周围的数据波动。
# mu中的bM比较结婚年龄相同时的结婚率差异，bA反过来比较。
# a、bM、bA的正态先验约束合理尺度；sigma的指数先验只支持正数。
# quap寻找后验峰值，再用峰顶曲率给出联合正态近似；不是MCMC抽样。
# 参数怎么填：dnorm(均值,标准差)，不是填方差；0.2和0.5对应标准化尺度。
# dexp(1)的1是速率，均值为1/速率；不能把它理解为固定sigma=1。
# data=dat把公式中的D/M/A与数据对应。start只是优化起点，
# 依次给a、bM、bA、sigma初值；sigma必须为正，初值不决定先验均值。
# 调整先验应改dnorm/dexp中的参数，并重新检查先验预测；不要只改start。
fit <- quap(alist(D ~ dnorm(mu,sigma), mu <- a+bM*M+bA*A,
                 a ~ dnorm(0,0.2), bM ~ dnorm(0,0.5), bA ~ dnorm(0,0.5),
                 sigma ~ dexp(1)), data=dat, start=list(a=0,bM=0,bA=-0.5,sigma=0.8))
# 近似分布以峰值为中心：此处mean是近似后验均值。
# vcov保留参数间协方差，区间同时反映先验与数据的信息。
means <- coef(fit)[c("a","bM","bA","sigma")]
cv <- vcov(fit)[names(means),names(means)]
sds <- sqrt(diag(cv))
# 区间怎么填：中央概率为c时，上分位点是(1+c)/2。
# 本例c=0.89，所以填0.945；若改95%，填0.975，并同步改列名和图注。
q <- qnorm(.945) # 中间89%，两端各5.5%
result <- data.frame(mean=means,sd=sds,lower89=means-q*sds,upper89=means+q*sds)
write.csv(result,file.path(out,"posterior.csv"),row.names=FALSE)
# 结果怎么看：posterior.csv四行依次为a、bM、bA、sigma，sd是参数的后验标准差。
# 其中sigma本身是数据的残差标准差；sigma那一行的sd是对sigma估计的不确定性。
# bM约-0.065，89%区间约[-0.306,0.176]：年龄相同后，方向仍不确定；
# 区间跨0不证明作用为0。bA约-0.614，表示每增加1个年龄标准差，
# 平均离婚率降低约0.614个离婚率标准差，前提是结婚率相同。
print(result)
cat("bM是在年龄相同条件下的斜率；因果解释仍依赖图和模型假设。\n")
# 给出所有0/1组合，而不是用随机次数估计概率。
# 下面切换到已知概率的教学模型，X/Y/Z/A是0或1状态；这里的A不是结婚年龄。
# bern(value,p)的p填“该变量为1的概率”，应在0到1之间。
# 0.1+0.8*z：z=0时概率0.1，z=1时概率0.9；乘法来自生成规则的条件独立。
# 这只是枚举所有数据状态，不是在参数网格上近似后验。
bern <- function(value,p) ifelse(value==1,p,1-p)
# model：1叉、2管、3对撞、4中间变量的后代。
# 前三种的A只是独立占位变量，乘0.5后求和即可消去。
# 每一行的概率沿箭头相乘，全部行加总为1。
rows <- expand.grid(model=1:4,X=0:1,Z=0:1,Y=0:1,A=0:1)
rows$probability <- 0
for(i in seq_len(nrow(rows))) {
 x <- rows$X[i]; z <- rows$Z[i]; y <- rows$Y[i]; a <- rows$A[i]
 m <- rows$model[i]
 w <- switch(m,
  .5*bern(x,.1+.8*z)*bern(y,.1+.8*z)*.5,
  .5*bern(z,.1+.8*x)*bern(y,.1+.8*z)*.5,
  .25*bern(z,ifelse(x+y>0,.9,.2))*.5,
  .5*bern(z,.1+.8*x)*bern(y,.1+.8*z)*bern(a,.1+.8*z))
 rows$probability[i] <- w
}
answers <- data.frame()
for(m in 1:4) for(g in 0:4) {
 sub <- rows[rows$model==m,]
 if(g>0) {
  col <- if(g<3) "Z" else "A"
  sub <- sub[sub[[col]]==(g-1)%%2,]
 }
# 分组后除以组内总概率，得到条件概率；据此计算E(XY)-E(X)E(Y)。
 w <- sub$probability/sum(sub$probability)
 ex <- sum(w*sub$X); ey <- sum(w*sub$Y)
 corr <- (sum(w*sub$X*sub$Y)-ex*ey)/sqrt(ex*(1-ex)*ey*(1-ey))
 answers <- rbind(answers,data.frame(model=m,group=g,correlation=corr))
}
write.csv(rows,file.path(out,"joint.csv"),row.names=FALSE)
# 枚举输出怎么看：model=1/2/3/4分别为叉/管/对撞/后代。
# group=0全体、1固定Z=0、2固定Z=1、3固定A=0、4固定A=1。
# 后代模型全体相关约0.64，固定Z后为0，固定A后约0.390。
# 练习：只把后代模型最后的bern(A,0.1+0.8*Z)概率改成0.5
# （原版变量名是小写a、z）；固定A后的相关应回到0.64，其余机制保持不变。
write.csv(answers,file.path(out,"associations.csv"),row.names=FALSE)
# 原书171页生成规则的期望计算，不是一次样本拟合。
# 植物参数：treatment的0/1分别是不处理/处理；0.5和0.1是出现真菌的概率。
# 5与2是无真菌/有真菌时的平均增长。按两种状态的概率加权，
# 得到3.5与4.7，两者之差1.2是生成规则给出的理论总效应。
plant <- data.frame(treatment=0:1,fungus_probability=c(.5,.1))
plant$mean_growth <- (1-plant$fungus_probability)*5+plant$fungus_probability*2
write.csv(plant,file.path(out,"plant_expectation.csv"),row.names=FALSE)
print(plant)

# 两张中文图在RStudio Plots面板显示，同时用ragg保存。
library(systemfonts)
library(ragg)
fonts <- unique(system_fonts()$family)
font <- intersect(c("PingFang SC","Microsoft YaHei","Noto Sans CJK SC","Heiti SC"),fonts)[1]
if(is.na(font)) stop("请安装Noto Sans CJK SC中文字体。")
# 图怎么读：点是近似后验均值，横线是89%后验区间，竖虚线为零。
# 区间宽度反映参数不确定性，不是地区之间的离婚率波动。
plot_model <- function() {
 par(family=font,mar=c(5,7,3,1))
 plot(means[2:3],c(2,1),xlim=c(-1,.4),ylim=c(.5,2.5),yaxt="n",pch=19,col="#167e83",
 xlab="标准化斜率及89%后验区间",ylab="",main="年龄相同以后，结婚率的斜率接近零")
 axis(2,at=c(2,1),labels=c("结婚率","结婚年龄"),las=1)
 segments(means[2:3]-q*sds[2:3],c(2,1),means[2:3]+q*sds[2:3],c(2,1),col="#167e83",lwd=3)
 abline(v=0,lty=2,col="gray")
}
plot_descendant <- function() {
 par(family=font)
 vals <- answers$correlation[answers$model==4 & answers$group %in% c(0,1,3)]
 barplot(vals,names.arg=c("全体","固定中间变量Z","固定后代A"),ylim=c(0,.75),col="#167e83",ylab="理论相关",main="后代携带中间变量的信息")
}
for(name in c("model","descendant")) {
 fun <- get(paste0("plot_",name))
 if(interactive()) fun()
 agg_png(file.path(out,paste0(name,".png")),width=1500,height=900,res=180)
 fun()
 dev.off()
}
