# 配套绘图文件：由 normal_chisq_t_f.R 在最后调用。
# 初学者先读主脚本的统计计算，再按需要阅读这里的绘图参数。
# 使用主脚本已经设置好的 out、figs、dfs、q；不单独运行。
# png() 打开图片文件；plot()/lines() 画图；dev.off() 完成文件写入。
# 横轴是统计量取值，纵轴是概率密度；曲线高度不是概率。
# 下方每一块先生成坐标并保存 CSV，再生成对应图片。

# 二、导出配图所用的全部曲线坐标，供其他软件复画。
x <- seq(-5, 5, length.out = 1001)
curves <- data.frame(x = x, normal = dnorm(x), t_df1 = dt(x, 1),
                     t_df5 = dt(x, 5), t_df30 = dt(x, 30))
write.csv(curves, file.path(out, 't_normal_curves.csv'), row.names = FALSE)
png(file.path(figs, '02-t分布与标准正态比较.png'), width = 1600, height = 1000, res = 160)
par(mar = c(4.5, 4.5, 3, 1), las = 1)
matplot(x, curves[, -1], type = 'l', lty = c(2, 1, 1, 1),
        col = c('#222222', '#178C88', '#4144A5', '#E36854'), lwd = 2,
        xlab = 'Statistic value', ylab = 'Density',
        main = 'Student t and standard normal distributions')
legend('topright', c('N(0, 1)', 't: df = 1', 't: df = 5', 't: df = 30'),
       col = c('#222222', '#178C88', '#4144A5', '#E36854'),
       lty = c(2, 1, 1, 1), lwd = 2, bty = 'n', cex = .8)
dev.off()

# 二补充、卡方曲线：分面避免 df=1 在零附近的高密度压扁其他曲线。
chi_x <- seq(.001, 55, length.out = 3001)
chi_df <- c(1, 2, 3, 5, 10, 30)
chi_curves <- data.frame(x = chi_x)
for (k in chi_df) chi_curves[[paste0('df_', k)]] <- dchisq(chi_x, k)
write.csv(chi_curves, file.path(out, 'chisq_curves.csv'), row.names = FALSE)
png(file.path(figs, '03-卡方分布是一族曲线.png'), width = 1600, height = 1500, res = 160)
par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1), las = 1)
plot(chi_x, dchisq(chi_x, 1), type = 'l', col = '#178C88', lwd = 2,
     xlim = c(0, 8), ylim = c(0, 1.5), xlab = 'Chi-square value', ylab = 'Density',
     main = 'Chi-square: df = 1 and 2')
lines(chi_x, dchisq(chi_x, 2), col = '#4144A5', lwd = 2, lty = 2)
legend('topright', c('df = 1', 'df = 2'), col = c('#178C88', '#4144A5'), lty = c(1, 2), lwd = 2, bty = 'n')
mtext('df = 1: density tends to infinity as x approaches 0; upper part clipped', side = 3, cex = .7, line = .2)
cols <- c('#178C88', '#4144A5', '#E36854', '#333333')
matplot(chi_x, chi_curves[, c('df_3', 'df_5', 'df_10', 'df_30')], type = 'l',
        col = cols, lty = 1:4, lwd = 2, xlab = 'Chi-square value', ylab = 'Density',
        main = 'Chi-square: increasing degrees of freedom')
legend('topright', paste('df =', c(3, 5, 10, 30)), col = cols, lty = 1:4, lwd = 2, bty = 'n')
dev.off()

# F 曲线：分别固定分母和分子自由度，观察另一个自由度的影响。
f_x <- seq(.001, 6, length.out = 2001)
pairs <- data.frame(d1 = c(3, 5, 10, 5, 5), d2 = c(10, 10, 10, 5, 30))
f_curves <- data.frame(x = f_x)
for (j in seq_len(nrow(pairs))) {
  name <- paste0('F_', pairs$d1[j], '_', pairs$d2[j])
  f_curves[[name]] <- df(f_x, pairs$d1[j], pairs$d2[j])
}
write.csv(f_curves, file.path(out, 'f_curves.csv'), row.names = FALSE)
png(file.path(figs, '04-F分布的两个自由度.png'), width = 1600, height = 1500, res = 160)
par(mfrow = c(2, 1), mar = c(4, 4.5, 3, 1), las = 1)
for (panel in 1:2) {
  keys <- if (panel == 1) c('F_3_10', 'F_5_10', 'F_10_10') else c('F_5_5', 'F_5_10', 'F_5_30')
  labs <- if (panel == 1) c('F(3, 10)', 'F(5, 10)', 'F(10, 10)') else c('F(5, 5)', 'F(5, 10)', 'F(5, 30)')
  matplot(f_x, f_curves[, keys], type = 'l', col = cols[1:3], lty = 1:3, lwd = 2,
          xlab = 'F value', ylab = 'Density',
          main = if (panel == 1) 'Fixed denominator df = 10' else 'Fixed numerator df = 5')
  abline(v = 1, col = '#999999', lty = 3)
  legend('topright', labs, col = cols[1:3], lty = 1:3, lwd = 2, bty = 'n')
}
dev.off()

# 构造关系图：每一行是独立的一条构造；不是对原始资料的正态性要求。
png(file.path(figs, '01-三个分布的构造关系.png'), width=1800, height=1300, res=180)
par(mar=c(1,1,3,1));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1))
title('Normal, chi-square, t and F: construction', cex.main=1.2)
yrows <- c(.80,.51,.22)
inputs <- c('Independent Z[1], ..., Z[k] ~ N(0,1)',
            'Z ~ N(0,1), U ~ chi-square(nu); independent',
            'U ~ chi-square(d1), V ~ chi-square(d2); independent')
forms <- expression(U == sum(Z[i]^2, i==1,k) %~% chi^2 * plain("(") * k * plain(")"),
                    T == frac(Z,sqrt(U/nu)) %~% t(nu),
                    F == frac(U/d[1],V/d[2]) %~% F(d[1],d[2]))
for (j in 1:3) {
 rect(.025,yrows[j]-.10,.975,yrows[j]+.13,col='#f0f7f6',border='#c3ddda')
 text(.5,yrows[j]+.075,inputs[j],cex=.90,col='#29424d')
 text(.5,yrows[j]-.035,forms[j],cex=1.4,col='#167d80')
}
text(.5,.025,expression(T %~% t(nu) ~~ implies ~~ T^2 %~% F(1,nu)),cex=1.1)
dev.off()

# 理论第一类错误率，不是模拟比例。t 界值一列应为 0.05。
error_rates <- data.frame(df=dfs, n_one_sample=dfs+1,
 correct_t=2*pt(qt(.975,dfs),dfs,lower.tail=FALSE),
 wrong_1_96=2*pt(1.96,dfs,lower.tail=FALSE))
stopifnot(max(abs(error_rates$correct_t-.05))<1e-10)
write.csv(error_rates,file.path(out,'error_rates.csv'),row.names=FALSE)
png(file.path(figs,'05-错用正态界值的实际错误率.png'),width=1600,height=1000,res=160)
par(mar=c(4.5,4.5,4,1),las=1)
bp <- barplot(100*t(as.matrix(error_rates[,c('correct_t','wrong_1_96')])),
 beside=TRUE,names.arg=paste0('df=',dfs),col=c('#178C88','#E36854'),
 ylim=c(0,38),ylab='Type I error (%)',xlab='Degrees of freedom',
 main='When a t statistic is compared with 1.96')
abline(h=5,lty=3,col='#555555')
text(bp,100*t(as.matrix(error_rates[,c('correct_t','wrong_1_96')]))+1.1,
 labels=sprintf('%.2f',100*t(as.matrix(error_rates[,c('correct_t','wrong_1_96')]))),cex=.7)
legend('topright',c('Correct t cutoff','Wrong cutoff: 1.96'),fill=c('#178C88','#E36854'),bty='n')
mtext('Exact probabilities under H0; independent normal sample, unknown sigma',side=3,line=.3,cex=.75)
dev.off()

