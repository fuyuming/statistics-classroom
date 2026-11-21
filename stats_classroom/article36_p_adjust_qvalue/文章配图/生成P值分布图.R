# 教师绘图入口：命令行依次传入共享P值CSV和输出PNG路径。
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dat <- read.csv(args[1])
p <- dat$p
lambda <- .5
pi0_hat <- min(1, mean(p >= lambda)/(1-lambda))

# 先独立统计每个宽0.1区间的人数，再把结果交给ggplot2绘图。
bin <- cut(p, breaks = seq(0, 1, by=.1), right=FALSE, include.lowest=TRUE)
counts <- as.integer(table(bin))
plot_data <- data.frame(midpoint=seq(.05,.95,by=.1), count=counts)
plot_data$tail <- plot_data$midpoint >= lambda
expected_per_bin <- length(p)*pi0_hat*.1

library(ggplot2)
g <- ggplot(plot_data, aes(x=midpoint,y=count,fill=tail)) +
  geom_col(width=.096, show.legend=FALSE) +
  scale_fill_manual(values=c('FALSE'='#df9b48','TRUE'='#348dc4')) +
  geom_hline(yintercept=expected_per_bin, linetype='dashed', linewidth=.65, color='#27343a') +
  geom_vline(xintercept=lambda, linetype='dotted', linewidth=.6, color='#27343a') +
  annotate('text', x=.76, y=max(counts)*.90,
           label=paste0('P >= 0.5: ',sum(p>=lambda)), size=5, color='#27343a') +
  scale_x_continuous(breaks=seq(0,1,by=.2),limits=c(0,1),expand=c(0,0)) +
  scale_y_continuous(expand=expansion(mult=c(0,.06))) +
  labs(x='P value', y='Count') +
  theme_classic(base_size=17) +
  theme(text=element_text(color='#27343a'),axis.text=element_text(color='#27343a'),
        plot.margin=margin(12,12,8,8))
ggsave(args[2],g,width=5.2,height=3.7,dpi=240,bg='white')
write.csv(plot_data,sub('\\.png$','_坐标.csv',args[2]),row.names=FALSE)
