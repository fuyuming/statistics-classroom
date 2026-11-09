# 配对率差区间：正文R函数与三语言共同算法（UTF-8）。
# 在本篇R项目中Source本文件，先安装ExactCIdiff与PropCIs。
# install.packages(c("ExactCIdiff","PropCIs"))
root <- if (basename(getwd())=="代码") normalizePath("..") else getwd()
output <- file.path(root,"运行结果","R")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
required <- c("ExactCIdiff","PropCIs")
if (!all(vapply(required,requireNamespace,logical(1),quietly=TRUE))) stop("请先安装ExactCIdiff与PropCIs。")

# 一、小样本20人：Wang归纳构造的精确区间。
# 1个+1、14个0、5个-1，方向后测减前测。
wang <- ExactCIdiff::PairedCI(n12=1,t=14,n21=5,conf.level=.95,CItype="Two.sided",grid.one=60,grid.two=40)
print(wang)
capture.output(print(wang),file=file.path(output,"PairedCI_实际打印.txt"))
write.csv(data.frame(method="Wang_PairedCI",estimate=wang$estimate,lower=wang$ExactCI[1],upper=wang$ExactCI[2]),file.path(output,"Wang精确区间.csv"),row.names=FALSE)

# 二、大样本1000人：配对Wald近似，明确插件标准误与计算方向。
n <- 1000
n_plus <- 50
n_minus <- 250
difference <- (n_plus-n_minus)/n
se <- sqrt(((n_plus+n_minus)/n-difference^2)/n)
wald <- difference+c(-1,1)*qnorm(.975)*se
print(wald)
# scoreci.mp的源代码计算(c-b)/n，故c填+1人数，b填-1人数。
score <- PropCIs::scoreci.mp(b=250,c=50,n=1000,conf.level=.95)
print(score)
capture.output(print(score),file=file.path(output,"Tango_scoreCI_实际打印.txt"))
write.csv(data.frame(method=c("paired_Wald","Tango_score"),estimate=difference,lower=c(wald[1],score$conf.int[1]),upper=c(wald[2],score$conf.int[2])),file.path(output,"大样本区间.csv"),row.names=FALSE)

# 三、三语言共同演示：固定D排序的三项无条件精确区间。
# 与PairedCI的排序/区间构造不同，不能把两个输出混称。
source(file.path(root,"代码","paired_exact.R"))
common <- paired_exact(1,14,5)
print(common)
write.csv(common,file.path(output,"共同精确区间.csv"),row.names=FALSE)
print(c(ExactCIdiff=as.character(packageVersion("ExactCIdiff")),PropCIs=as.character(packageVersion("PropCIs"))))
