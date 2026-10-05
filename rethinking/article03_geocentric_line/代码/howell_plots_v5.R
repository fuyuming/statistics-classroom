# 由 howell_v5.R 调用。所有图使用刚算出的 R 结果，不读取 Python 图片。
# 自动寻找实际已安装的中文字体；缺少字体时明确提示，不静默生成方框。
if (!requireNamespace('systemfonts',quietly=TRUE) || !requireNamespace('ragg',quietly=TRUE)) {
  stop('中文绘图需要 systemfonts 和 ragg；请先运行 install.packages(c("systemfonts", "ragg"))')
}
font_candidates <- c('PingFang SC','Microsoft YaHei','SimHei','Noto Sans CJK SC','Source Han Sans SC','WenQuanYi Micro Hei','Heiti SC')
installed <- unique(systemfonts::system_fonts()$family)
matched <- font_candidates[font_candidates %in% installed]
if (!length(matched)) stop('未找到中文字体，请安装 Noto Sans CJK SC，重启 RStudio 后重跑。')
cn_font <- matched[1]
cat('中文绘图字体：',cn_font,'\n')
plot_dir <- file.path(out,'R_figures')
dir.create(plot_dir,showWarnings=FALSE)
teal <- '#167e83'; orange <- '#ca6841'
hseq <- seq(min(d$height),max(d$height),length.out=60)
xseq <- hseq-dat$Hbar
center <- summary$mean
avg <- center[1]+center[2]*xseq
# 给每张图同样的保存/显示行为：先存PNG，回到原有设备再画一次。
# RStudio 下使用现有 RStudioGD；命令行无GUI时仅保存文件。
render_plot <- function(name, draw) {
  ragg::agg_png(file.path(plot_dir,paste0(name,'.png')),width=1200,height=1000,res=130)
  par(family=cn_font)
  tryCatch(draw(),finally=dev.off())
  if (interactive()) {
    par(family=cn_font)
    draw()
  }
}
render_plot('01_relationship',function() {
  plot(d$height,y,pch=16,col='#91a5a680',xlab='身高（厘米）',ylab='体重（千克）',main='01  身高与体重的平均关系')
  lines(hseq,avg,col=teal,lwd=3)
})
ix <- which.min(abs(d$height-160)+.03*abs(y-(center[1]+center[2]*x)-6))
mx <- center[1]+center[2]*x[ix]
render_plot('02_residual',function() {
  plot(d$height,y,pch=16,col='#91a5a650',xlab='身高（厘米）',ylab='体重（千克）',main='02  一个人的残差')
  lines(hseq,avg,col=teal,lwd=2)
  segments(d$height[ix],mx,d$height[ix],y[ix],col=orange,lwd=3)
  points(rep(d$height[ix],2),c(mx,y[ix]),col=c(teal,orange),pch=16,cex=1.5)
  legend('topleft',legend=sprintf('实际体重 %.2f；模型均值 %.2f 千克',y[ix],mx),bty='n')
})
render_plot('03_conditional_normal',function() {
  xx <- seq(30,70,length.out=400)
  mm <- center[1]+center[2]*(160-dat$Hbar)
  plot(xx,dnorm(xx,mm,center[3]),type='l',lwd=3,col=teal,xlab='体重（千克）',ylab='概率密度',main='03  身高160厘米：固定在后验平均参数处')
  abline(v=mm,lty=2,col=orange)
})
set.seed(20261005)
prior_a <- rnorm(40,60,10); prior_b <- rlnorm(40,0,1)
selected <- sample(seq_len(nrow(g)),40,replace=TRUE,prob=w)
render_plot('04_prior_posterior',function() {
  oldpar <- par(mfrow=c(2,1),mar=c(4,4,3,1))
  on.exit(par(oldpar))
  for(k in 1:2) {
    plot(d$height,y,pch=16,col='#ca684150',xlim=range(hseq),ylim=c(-150,300),xlab='身高（厘米）',ylab='体重（千克）',main=c('04  看数据前的均值线','看数据后的均值线')[k])
    aa <- if(k==1) prior_a else g$a[selected]
    bb <- if(k==1) prior_b else g$b[selected]
    for(j in 1:40) lines(hseq,aa[j]+bb[j]*xseq,col='#167e8340')
  }
})
render_plot('05_grid_quap',function() {
  oldpar <- par(mfrow=c(3,1),mar=c(4,4,2,1),mgp=c(2.2,.7,0))
  on.exit(par(oldpar))
  for(j in 1:3) {
    values <- sort(unique(g[[j]]))
    mass <- as.numeric(tapply(w,g[[j]],sum))
    plot(values,mass/(values[2]-values[1]),type='l',col=teal,lwd=2,xlab=c('截距 α（千克）','斜率 β（千克/厘米）','标准差 σ（千克）')[j],ylab='概率密度',main=c('05  网格与 quap 后验对照','','')[j])
    lines(values,dnorm(values,coef(fit)[j],sqrt(vcov(fit)[j,j])),col=orange,lwd=2,lty=2)
    if(j==1) legend('topright',c('网格','quap'),col=c(teal,orange),lty=c(1,2),bty='n')
  }
})
# 使用联合 quap 抽样，保留参数相关性；R模拟端点允许抽样波动。
set.seed(20261005)
pp <- extract.samples(fit,n=20000)
mu <- outer(pp$b,xseq,'*')+pp$a
new_weight <- matrix(rnorm(length(mu),as.vector(mu),rep(pp$sigma,ncol(mu))),nrow=nrow(mu))
ci <- apply(mu,2,quantile,probs=c(.055,.945))
pi <- apply(new_weight,2,quantile,probs=c(.055,.945))
render_plot('06_prediction',function() {
  plot(hseq,colMeans(mu),type='n',ylim=range(pi),xlab='身高（厘米）',ylab='体重（千克）',main='06  quap：平均与新个体的89%区间')
  polygon(c(hseq,rev(hseq)),c(pi[1,],rev(pi[2,])),col='#cce3df',border=NA)
  polygon(c(hseq,rev(hseq)),c(ci[1,],rev(ci[2,])),col=teal,border=NA)
  lines(hseq,colMeans(mu),lwd=2)
  legend('topleft',c('新个体','平均体重'),fill=c('#cce3df',teal),bty='n')
})
cat('六张原生 R 图已保存到运行结果/v5/R_figures；交互运行时也显示在当前绘图设备。\n')
