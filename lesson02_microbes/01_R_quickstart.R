# 模拟教学数据：同一非致病微生物在两类培养基中的固定终点OD600。
# 每行假设为一个独立培养物；不是同一培养物的技术重复，也不是时间序列。
# OD600已按教学设定扣除空白，是无量纲浊度读数，不直接等于活菌数或生长速率。
# 明哥的微生物世界｜R 快速入门｜2026-09-17
# 自拟教学数据，不是真实实验，不进行假设检验。
# 打开同目录 course_demo.Rproj，然后 source("01_R_quickstart.R")。
# 本脚本只用基础R，无需安装额外包。重复运行只覆盖 output/R 中的练习结果。

# 1. 定位文件：命令行、source()自动识别脚本；逐行运行时使用项目目录。
source_files <- unlist(lapply(sys.frames(), function(x) x$ofile))
cli <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])
file_path <- if(length(source_files)) tail(source_files,1) else cli
root <- if(length(file_path)) dirname(normalizePath(file_path)) else getwd()
out <- file.path(root,"output","R")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# 2. 导入、检查：空白或NA作为缺失；字符列不自动转因子。
microbes <- read.csv(file.path(root,"data","microbes_od600_demo.csv"),
                   na.strings=c("","NA"),stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
print(head(microbes)); str(microbes)
# 检查数据结构，避免误把文件读错还继续分析。
stopifnot(all(c("culture_id","group","od600") %in% names(microbes)),
          is.numeric(microbes$od600), !anyDuplicated(microbes$culture_id),
          !anyNA(microbes$group), all(microbes$group %in% c("medium_A","medium_B")),
          all(is.finite(microbes$od600) | is.na(microbes$od600)))
microbes$group <- factor(microbes$group,levels=c("medium_A","medium_B"))
# 缺失不等于0；先显示缺失，再在每组计算时暂时去掉缺失值。
print(colSums(is.na(microbes)))
print(microbes[is.na(microbes$od600),])

# 3. 定义一个汇总函数：输入一组数值，输出该组统计量。
# n_total为记录数，n为有效观测数，n_missing为缺失数。
# 样本标准差分母为n-1；四分位数使用type=7，三种语言统一此定义。
summarise_one <- function(raw) {
  x <- raw[!is.na(raw)]
  if(length(x)<2) stop("每组至少需要2个有效观测，才能计算样本标准差。")
  q <- quantile(x,c(.25,.75),type=7,names=FALSE)
  data.frame(n_total=length(raw),n=length(x),n_missing=sum(is.na(raw)),
    mean=mean(x),sd=sd(x),median=median(x),q1=q[1],q3=q[2],
    iqr=q[2]-q[1],min=min(x),max=max(x))
}
# lapply对每组调用同一个函数；do.call把结果按行拼接。
groups <- levels(microbes$group)
summary_table <- do.call(rbind,lapply(groups,function(g) {
  cbind(group=g,summarise_one(microbes$od600[microbes$group==g]))
}))
rownames(summary_table) <- NULL
print(summary_table)
write.csv(summary_table,file.path(out,"summary.csv"),row.names=FALSE)
# 保留数据的类型；CSV只保存表格值，RDS还保留因子信息。
saveRDS(microbes,file.path(out,"microbes.rds"))

# 4. 散点图：每一点是一个有效观测，横线是样本均值。
# 为避免完全重叠，用固定的微小水平偏移；纵坐标从未改变。
cols <- c("#267D87","#D69645")
png(file.path(out,"group_plot.png"),width=1600,height=1000,res=180)
par(mar=c(5,5,4,2))
plot(NA,xlim=c(.6,2.4),ylim=c(.3,.65),xaxt="n",xlab="Group",
     ylab="OD600 (blank-corrected)",main="Microbial OD600: synthetic teaching data")
axis(1,at=1:2,labels=c("Medium A","Medium B"))
for(idx in seq_along(groups)) {
  x <- microbes$od600[microbes$group==groups[idx]]
  x <- x[!is.na(x)]
  points(idx+seq(-.06,.06,length.out=length(x)),x,pch=19,col=cols[idx],cex=1.4)
  segments(idx-.13,mean(x),idx+.13,mean(x),lwd=3,col=cols[idx])
}
mtext("Synthetic classroom data; missing values omitted from plotted measurements",side=1,line=3.5,cex=.7)
dev.off()
# 5. 直方图：按组绘制，使用相同分箱边界，便于比较。
png(file.path(out,"histogram.png"),width=1600,height=900,res=180)
par(mfrow=c(1,2),mar=c(5,4,4,1))
for(idx in seq_along(groups)) {
  x <- microbes$od600[microbes$group==groups[idx]]
  hist(x[!is.na(x)],breaks=seq(.3,.65,by=.05),right=FALSE,
       col=cols[idx],border="white",main=groups[idx],xlab="OD600 (blank-corrected)",ylim=c(0,4))
}
dev.off()
# 6. 自检与版本记录：这里只检查练习是否按预期处理了数据。
stopifnot(sum(summary_table$n_total)==16,sum(summary_table$n)==14,
          sum(summary_table$n_missing)==2,all(summary_table$n==7))
capture.output(sessionInfo(),file=file.path(out,"versions.txt"))
cat("完成：统计表和两张图在",out,"\n")
