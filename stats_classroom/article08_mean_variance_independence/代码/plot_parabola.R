# X~N(0,1)、Y=X²：理论函数网格，不是抽样数据。
curve_x <- seq(-3,3,length.out=601)
curve_y <- curve_x^2
curve_product <- curve_x*(curve_y-1)
curve_density <- exp(-curve_x^2/2)/sqrt(2*pi)
parabola <- data.frame(x=curve_x,y=curve_y,centered_product=curve_product,
  normal_density=curve_density,weighted_integrand=curve_product*curve_density)
stopifnot(max(abs(curve_product+rev(curve_product)))<1e-12)
dir.create('运行结果/R',recursive=TRUE,showWarnings=FALSE)
write.csv(parabola,'运行结果/R/parabola_coordinates.csv',row.names=FALSE)
dir.create('文章配图/R',recursive=TRUE,showWarnings=FALSE)
png('文章配图/R/04_parabola.png',width=1260,height=1305,res=150)
par(mar=c(5,5,4,1),oma=c(3,0,0,0))
plot(curve_x,curve_y,type='l',lwd=3,col='#167d80',xlim=c(-3.2,3.2),ylim=c(0,10),
  xlab='X',ylab='Y = X squared',main='X ~ N(0,1), Y = X squared: dependent, covariance = 0')
abline(h=1,lty=2,col='gray')
points(c(-2,2),c(4,4),pch=19,col=c('#cc773e','#167d80'),cex=1.5)
text(c(-1.5,1.5),c(7,7),c('X = -2; Y = 4\nX(Y - 1) = -6','X = 2; Y = 4\nX(Y - 1) = 6'))
arrows(c(-1.5,1.5),c(6,6),c(-2,2),c(4.2,4.2),length=.1)
mtext('Each x and -x: opposite centered products, equal normal densities.',side=1,outer=TRUE,line=.5)
mtext('Plot limited to [-3,3]; covariance expectation is over the full real line.',side=1,outer=TRUE,line=2)
dev.off()
