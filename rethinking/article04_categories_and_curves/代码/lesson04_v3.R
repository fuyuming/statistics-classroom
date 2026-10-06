# 精读04：原作者分类模型与年龄—身高样条。RStudio 打开本课Rproj后 Source。
# 依赖 rethinking、systemfonts、ragg；splines 随R提供。
# 自动计算并显示六张中文图，也保存到运行结果/v3/R_figures。
library(rethinking)
library(splines)
z <- grep('^--file=',commandArgs(FALSE),value=TRUE)
root <- if(length(z)) dirname(dirname(normalizePath(sub('^--file=','',z[1])))) else getwd()
out <- file.path(root,'运行结果/v3'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
# 数据取自作者R包。共享CSV随项目提供，其他语言无需先运行R。
data(Howell1); d <- Howell1[Howell1$age>=18,]
write.csv(d,file.path(root,'数据/v3/Howell1_adults.csv'),row.names=FALSE)
dat <- list(W=d$weight,S=d$male+1,H=d$height,Hbar=mean(d$height))
print(head(d)); cat('1表示女性，2表示男性；编号是标签，不是数量。\n')
# 模型一：每组一个均值，sigma仍未知。与作者m_SW相同。
m1 <- quap(alist(W~dnorm(mu,sigma),mu<-a[S],a~dnorm(60,10),sigma~dunif(0,10)),data=dat,start=list(a=c(42,49),sigma=5))
# 模型二：每组一条线，比较同身高处的平均体重。按PPT第59页：斜率Uniform(0,1)；公开脚本另有LogNormal版本。
m2 <- quap(alist(W~dnorm(mu,sigma),mu<-a[S]+b[S]*(H-Hbar),a~dnorm(60,10),b~dunif(0,1),sigma~dunif(0,10)),data=dat,start=list(a=c(45,45),b=c(.5,.6),sigma=4))
# 曲线：PPT第84–92页的年龄—身高例子，数值设定见作者04_prior_pred_spline.r后半段。
data(Howell1)
write.csv(Howell1,file.path(root,'数据/v3/Howell1_all.csv'),row.names=FALSE)
ch <- Howell1[,c('age','height')]
knots <- seq(min(ch$age),max(ch$age),length.out=20)
B <- bs(ch$age,knots=knots,degree=3,intercept=FALSE)
ms <- quap(alist(Y~dnorm(mu,exp(log_sigma)),mu<-a0+as.vector(a%*%B),a0~dnorm(120,1),a~dnorm(0,25),log_sigma~dnorm(0,.5)),data=list(Y=ch$height,B=t(B)),start=list(a0=120,a=rep(0,ncol(B)),log_sigma=2),control=list(maxit=20000,reltol=1e-12))
write.csv(ch,file.path(root,'数据/v3/age_height.csv'),row.names=FALSE)
write.csv(B,file.path(root,'数据/v3/age_basis.csv'),row.names=FALSE)
curve_year <- seq(min(ch$age),max(ch$age),length.out=300)
Bp <- predict(B,curve_year)
write.csv(data.frame(age=curve_year,Bp),file.path(root,'数据/v3/age_basis_prediction.csv'),row.names=FALSE)
write.csv(data.frame(knot=knots),file.path(root,'数据/v3/age_knots.csv'),row.names=FALSE)
# 系数和协方差可完整核对；quap均值就是其正态近似的中心。
for(nm in c('m1','m2','ms')) {
 f <- get(nm)
 write.csv(data.frame(parameter=names(coef(f)),mean=coef(f),sd=sqrt(diag(vcov(f)))),file.path(out,paste0('R_',nm,'.csv')),row.names=FALSE)
 write.csv(vcov(f),file.path(out,paste0('R_',nm,'_cov.csv')),row.names=FALSE)
}
# PPT第97页Bonus：把身高和体重一起写入联合模型。
full <- quap(alist(W~dnorm(mu,sigma),mu<-a[S]+b[S]*(H-Hbar),a~dnorm(60,10),b~dunif(0,1),sigma~dunif(0,10),H~dnorm(nu,tau),nu<-h[S],h~dnorm(160,10),tau~dunif(0,10)),data=dat,start=list(a=c(45,45),b=c(.65,.61),sigma=4,h=c(150,160),tau=6),control=list(maxit=10000,reltol=1e-12))
write.csv(data.frame(parameter=names(coef(full)),mean=coef(full),sd=sqrt(diag(vcov(full)))),file.path(out,'R_full.csv'),row.names=FALSE)
set.seed(20261005);post_full <- extract.samples(full,n=50000)
meanS1 <- post_full$a[,1]+post_full$b[,1]*(post_full$h[,1]-dat$Hbar)
meanS2 <- post_full$a[,2]+post_full$b[,2]*(post_full$h[,2]-dat$Hbar)
# 对每套参数积分掉身高和体重噪声，得到总的平均效应。
total_mean <- meanS2-meanS1
# 两个独立模拟个体之差，不能当成同一个人的两种反事实之差。
H1 <- rnorm(50000,post_full$h[,1],post_full$tau)
H2 <- rnorm(50000,post_full$h[,2],post_full$tau)
person1 <- rnorm(50000,post_full$a[,1]+post_full$b[,1]*(H1-dat$Hbar),post_full$sigma)
person2 <- rnorm(50000,post_full$a[,2]+post_full$b[,2]*(H2-dat$Hbar),post_full$sigma)
write.csv(data.frame(total_mean=total_mean,independent_person_difference=person2-person1),file.path(out,'R_full_contrasts.csv'),row.names=FALSE)
source(file.path(root,'代码/lesson04_plots_v3.R'),local=TRUE)
