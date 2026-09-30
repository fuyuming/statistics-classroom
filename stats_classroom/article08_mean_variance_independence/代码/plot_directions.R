# 人工方向示例，不是随机抽样；箭头不能单独证明协方差为0。
# 在项目根目录运行；只需基础R。
direction_data <- read.csv('数据/direction_examples.csv')
direction_data$M <- (direction_data$x1+direction_data$x2)/2
direction_data$D <- direction_data$x1-direction_data$x2
stopifnot(all(direction_data$M==c(10,11,10,11)),all(direction_data$D==c(0,2,0,-2)))
dir.create('运行结果/R',recursive=TRUE,showWarnings=FALSE)
write.csv(direction_data,'运行结果/R/direction_coordinates.csv',row.names=FALSE)
print(direction_data)
dir.create('文章配图/R',recursive=TRUE,showWarnings=FALSE)
png('文章配图/R/03_directions.png',width=1260,height=1650,res=150)
par(mfrow=c(2,1),mar=c(4,5,4,2),oma=c(4,0,1,0))
for (i in 1:2) {
  rows <- if (i==1) 1:2 else 3:4
  a <- direction_data[rows,]
  color <- c('#167d80','#cc773e')[i]
  plot(a$M,a$D,xlim=c(9.7,11.6),ylim=c(-2.8,2.8),pch=19,col=color,
       xlab='Mean M',ylab='Difference D',main=c('Only X1: 10 to 12; X2 stays at 10',
       'Only X2: 10 to 12; X1 stays at 10')[i])
  arrows(a$M[1],a$D[1],a$M[2],a$D[2],length=.15,lwd=3,col=color)
  text(10,.4,'Start (10, 0)')
  text(11.05,a$D[2],paste0('(11, ',a$D[2],')'),pos=4,col=color)
  grid()
}
mtext('Under independence: Cov(M,D) = Var(X1)/2 - Var(X2)/2',outer=TRUE,side=1,line=1)
mtext('Equal variances imply zero covariance. Arrows are illustrative.',outer=TRUE,side=1,line=2.4)
dev.off()
