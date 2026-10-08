# 精读05 v2：RStudio打开精读05_v2.Rproj，再Source本文件。
# 作者PPT31页的quap模型；输入已附，无须先跑其他语言。
# 二元结构精确枚举，区别于PPT中1000次随机模拟的频数。
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
print(head(d))
# standardize减均值除样本标准差，与PPT一致。
dat <- list(D=standardize(d$Divorce), M=standardize(d$Marriage), A=standardize(d$MedianAgeMarriage))
fit <- quap(alist(D ~ dnorm(mu,sigma), mu <- a+bM*M+bA*A,
                 a ~ dnorm(0,0.2), bM ~ dnorm(0,0.5), bA ~ dnorm(0,0.5),
                 sigma ~ dexp(1)), data=dat, start=list(a=0,bM=0,bA=-0.5,sigma=0.8))
means <- coef(fit)[c("a","bM","bA","sigma")]
cv <- vcov(fit)[names(means),names(means)]
sds <- sqrt(diag(cv))
q <- qnorm(.945) # 中间89%，两端各5.5%
result <- data.frame(mean=means,sd=sds,lower89=means-q*sds,upper89=means+q*sds)
write.csv(result,file.path(out,"posterior.csv"),row.names=FALSE)
print(result)
cat("bM是在年龄相同条件下的斜率；因果解释仍依赖图和模型假设。\n")
# 给出所有0/1组合，而不是用随机次数估计概率。
bern <- function(value,p) ifelse(value==1,p,1-p)
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
 w <- sub$probability/sum(sub$probability)
 ex <- sum(w*sub$X); ey <- sum(w*sub$Y)
 corr <- (sum(w*sub$X*sub$Y)-ex*ey)/sqrt(ex*(1-ex)*ey*(1-ey))
 answers <- rbind(answers,data.frame(model=m,group=g,correlation=corr))
}
write.csv(rows,file.path(out,"joint.csv"),row.names=FALSE)
write.csv(answers,file.path(out,"associations.csv"),row.names=FALSE)
# 原书171页生成规则的期望计算，不是一次样本拟合。
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
 barplot(vals,names.arg=c("全体","固定中介Z","固定后代A"),ylim=c(0,.75),col="#167e83",ylab="理论相关",main="后代携带中介的信息")
}
for(name in c("model","descendant")) {
 fun <- get(paste0("plot_",name))
 if(interactive()) fun()
 agg_png(file.path(out,paste0(name,".png")),width=1500,height=900,res=180)
 fun()
 dev.off()
}
