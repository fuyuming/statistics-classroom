# 六张原生R图，当前RStudio设备显示，同时保存中文PNG。
if(!requireNamespace('systemfonts',quietly=TRUE)||!requireNamespace('ragg',quietly=TRUE)) stop('请安装 systemfonts、ragg')
families <- unique(systemfonts::system_fonts()$family)
font <- intersect(c('PingFang SC','Microsoft YaHei','SimHei','Noto Sans CJK SC','Heiti SC'),families)[1]
if(is.na(font)) stop('请安装 Noto Sans CJK SC 中文字体后重启。')
figdir <- file.path(out,'R_figures');dir.create(figdir,showWarnings=FALSE)
render <- function(name,draw) {
 ragg::agg_png(file.path(figdir,paste0(name,'.png')),width=1200,height=1000,res=140)
 par(family=font);tryCatch(draw(),finally=dev.off())
 if(interactive()) {par(family=font);draw()}
}
color <- c('#167e83','#ca6841');z89 <- qnorm(.945)
p1 <- coef(m1);v1 <- vcov(m1);p2 <- coef(m2);v3 <- vcov(m2)
render('01-个体与组平均',function() {
 set.seed(4);plot(jitter(dat$S,amount=.13),dat$W,col=adjustcolor(color[dat$S],.5),pch=16,xlim=c(.6,2.6),xaxt='n',xlab='组别',ylab='体重（千克）',main='01 个体体重与组平均')
 axis(1,1:2,c('女性','男性'))
 points((1:2)+.25,p1[1:2],pch=16,cex=1.5)
 arrows((1:2)+.25,p1[1:2]-z89*sqrt(diag(v1))[1:2],(1:2)+.25,p1[1:2]+z89*sqrt(diag(v1))[1:2],angle=90,code=3,length=.08,lwd=2)
})
set.seed(20261006);post <- extract.samples(m1,n=50000)
diffmean <- post$a[,2]-post$a[,1]
diffperson <- diffmean+rnorm(50000,0,post$sigma)-rnorm(50000,0,post$sigma)
render('02-均值差与个体差',function() {
 op<-par(mfrow=c(2,1),mar=c(4,4,3,1));on.exit(par(op))
 for(j in 1:2) {
  v<-if(j==1) diffmean else diffperson
  plot(density(v),xlim=c(-20,35),main=c('02 组平均体重之差','两个独立新个体的体重之差')[j],xlab='男性 − 女性（千克）',ylab='概率密度',col=color[1],lwd=2)
  abline(v=0,col=color[2],lty=2)
 }
})
render('03-每组一条线',function() {
 plot(dat$H,dat$W,col=adjustcolor(color[dat$S],.4),pch=16,xlab='身高（厘米）',ylab='体重（千克）',main='03 比较相同身高处的平均体重')
 for(j in 1:2) {
 hh<-seq(min(dat$H[dat$S==j]),max(dat$H[dat$S==j]),length.out=100)
 lines(hh,p2[j]+p2[j+2]*(hh-dat$Hbar),col=color[j],lwd=3)
 }
 legend('topleft',c('女性','男性'),col=color,lty=1,bty='n')
})
render('04-同身高平均差',function() {
 hh<-seq(max(tapply(dat$H,dat$S,min)),min(tapply(dat$H,dat$S,max)),length.out=100)
 C<-cbind(1,-1,hh-dat$Hbar,-(hh-dat$Hbar),0)
 mm<-as.vector(C%*%p2);ss<-sqrt(rowSums((C%*%v3)*C))
 plot(hh,mm,type='n',ylim=range(mm-z89*ss,mm+z89*ss),xlab='共同观测范围内身高（厘米）',ylab='女性 − 男性（千克）',main='04 同身高处均值差与89%可信区间')
 polygon(c(hh,rev(hh)),c(mm-z89*ss,rev(mm+z89*ss)),col=adjustcolor(color[1],.2),border=NA)
 lines(hh,mm,col=color[1],lwd=2);abline(h=0,lty=2,col=color[2])
})
render('05-样条积木',function() {
 op<-par(mfrow=c(2,1),mar=c(4,4,3,1));on.exit(par(op))
 matplot(curve_year,Bp[,9:11],type='l',lty=1,xlim=c(20,55),xlab='年龄（岁）',ylab='基函数的值',main='05 三块局部积木（示意）')
 plot(curve_year,as.vector(Bp[,9:11]%*%c(6,-3,5)),type='l',xlim=c(20,55),xlab='年龄（岁）',ylab='加权和',main='教学权重6、−3、5；不是拟合结果',col=color[1],lwd=2)
})
XP<-cbind(1,Bp);b<-coef(ms)[seq_len(ncol(XP))];V<-vcov(ms)[seq_len(ncol(XP)),seq_len(ncol(XP))]
mu<-as.vector(XP%*%b);se<-sqrt(rowSums((XP%*%V)*XP))
render('06-年龄身高曲线',function() {
 plot(ch$age,ch$height,pch=16,cex=.5,col='#91a5a670',xlab='年龄（岁）',ylab='身高（厘米）',main='06 年龄身高曲线：均值及89%可信区间')
 polygon(c(curve_year,rev(curve_year)),c(mu-z89*se,rev(mu+z89*se)),col=adjustcolor(color[1],.25),border=NA)
 lines(curve_year,mu,col=color[1],lwd=2)
})
cat('已生成六张中文R图；RStudio中可在Plots面板前后浏览。\n')
