# 统计课堂29｜独立样本率差 A-B（UTF-8）
# 问题：假设A组9/20阳性，B组13/20阳性，估计总体率差的双侧区间。
# RStudio打开本课Rproj，从准备节开始逐节运行；也可Source整文件。
# 依赖DescTools；默认控制台显示结果、Plots显示图，保存关闭。
# SAVE_RESULTS=TRUE才写运行结果/R。三种区间均调用R包现成函数。

# ---- 准备与参数 ----
confidence <- 0.95
SAVE_RESULTS <- FALSE
if (!requireNamespace("DescTools", quietly=TRUE)) {
  stop('请先在Console运行 install.packages("DescTools")')
}
# Rproj启动后工作目录为项目根；若当前在代码目录，则向上一层查找。
root <- getwd()
if (!file.exists(file.path(root, "数据", "independent_counts.csv"))) {
  root <- dirname(root)
}
data_path <- file.path(root, "数据", "independent_counts.csv")

# ---- 问题1：两批独立的人，A率减B率是多少？ ----
# 数据为假设的两行汇总计数；先查看字段与实际行，不补造个体记录。
data <- read.csv(data_path, fileEncoding="UTF-8", stringsAsFactors=FALSE)
cat("字段：", names(data), "\n")
print(data)
stopifnot(all(c("group", "positive", "total") %in% names(data)))
stopifnot(nrow(data)==2, setequal(data$group, c("A", "B")))
counts <- as.matrix(data[c("positive", "total")])
stopifnot(all(is.finite(counts)), all(counts==floor(counts)))
stopifnot(all(data$total>0), all(data$positive>=0), all(data$positive<=data$total))
stopifnot(length(confidence)==1, is.finite(confidence), confidence>0, confidence<1)
# [条件, ]按行选择A/B，保留全部列；每组恰有一行。
A <- data[data$group=="A", ]
B <- data[data$group=="B", ]
xA <- A$positive
nA <- A$total
xB <- B$positive
nB <- B$total
pA <- xA/nA
pB <- xB/nB
d <- pA-pB
cat(sprintf("A阳性/阴性=%d/%d；B阳性/阴性=%d/%d。\n", xA,nA-xA,xB,nB-xB))
cat(sprintf("A=%.2f%%，B=%.2f%%；A-B=%.2f个百分点。\n",100*pA,100*pB,100*d))

# ---- 问题2：Wald怎样从标准误得到区间？ ----
# 思路：独立两组的率差方差为各自方差相加；用样本率估计未知率。
alpha <- 1-confidence
# qnorm输入左侧累计概率；95%双侧区间取.975分位点，输出约1.96。
z <- qnorm(1-alpha/2)
se <- sqrt(pA*(1-pA)/nA+pB*(1-pB)/nB)
# c(-1,1)分别取“减去误差限”和“加上误差限”。
wald_by_hand <- d+c(-1,1)*z*se
cat(sprintf("Wald标准误=%.8f，正态临界值=%.6f。\n",se,z))
cat("Wald手算端点（百分点）：",100*wald_by_hand,"\n")

# ---- 问题3：用三个现成方法估计同一总体率差 ----
# 参数顺序A阳性、A总人数、B阳性、B总人数；估计第一组减第二组。
# 返回矩阵三列est、lwr.ci、upr.ci分别是原始差、下限、上限。
wald <- DescTools::BinomDiffCI(xA,nA,xB,nB,conf.level=confidence,method="wald")
# 此包的score指Newcombe无连续性校正；Python同构造名为newcomb。
newcombe <- DescTools::BinomDiffCI(xA,nA,xB,nB,conf.level=confidence,method="score")
ac <- DescTools::BinomDiffCI(xA,nA,xB,nB,conf.level=confidence,method="ac")
# rbind将每种方法的一行上下限排在一起；按列名取值，避免混淆est与中心。
intervals <- rbind(wald,newcombe,ac)
method_names <- c("Wald","Newcombe","Agresti-Caffo")
results <- data.frame(method=method_names, estimate=intervals[,"est"],
                      lower=intervals[,"lwr.ci"], upper=intervals[,"upr.ci"])
results$lower_pp <- 100*results$lower
results$upper_pp <- 100*results$upper
results$width_pp <- results$upper_pp-results$lower_pp
print(results, digits=9, row.names=FALSE)
# 每组各加1阳性1阴性，用于计算区间；CSV原始计数保持原值。
pA_adjusted <- (xA+1)/(nA+2)
pB_adjusted <- (xB+1)/(nB+2)
ac_center <- pA_adjusted-pB_adjusted
cat(sprintf("AC计算中心=%.4f个百分点；est仍为原始率差%.2f。\n",100*ac_center,100*d))
cat(sprintf("Newcombe的%.0f%%区间为%.2f至%.2f个百分点。\n",
            100*confidence,100*newcombe[1,"lwr.ci"],100*newcombe[1,"upr.ci"]))
if (newcombe[1,"lwr.ci"]<=0 && newcombe[1,"upr.ci"]>=0) {
  cat("该区间包含0，尚不能确定总体率差方向；不等于证明两总体率相等。\n")
} else {
  cat("该区间不含0，方向由区间正负判断；差异大小仍需结合研究背景。\n")
}
cat("这些是近似区间，不保证所有真实率下都恰好达到所标置信水平。\n")

# ---- 问题4：把区间与原始差值放在一张图中 ----
# 英文图签便于跨平台显示；横轴乘100，单位为百分点。
plot_intervals <- function() {
  par(mar=c(5,9,4,1))
  plot(NA,xlim=range(c(0,results$lower_pp,results$upper_pp)),ylim=c(.5,3.5),
       yaxt="n",ylab="",xlab="Population difference A-B (percentage points)",
       main=sprintf("%.0f%% two-sided CIs | A %d/%d, B %d/%d",100*confidence,xA,nA,xB,nB))
  axis(2,at=3:1,labels=method_names,las=1)
  abline(v=0,lty=2,col="gray")
  segments(results$lower_pp,3:1,results$upper_pp,3:1,col="#0000F5",lwd=3)
  points(rep(100*d,3),3:1,pch=16)
  points(100*ac_center,1,pch=18,col="#E89623",cex=1.4)
  legend("topright",legend=c("Raw A-B","AC center"),pch=c(16,18),
         col=c("black","#E89623"),bty="n",cex=.8)
}
plot_intervals()
if (SAVE_RESULTS) {
  output <- file.path(root,"运行结果","R")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  write.csv(results,file.path(output,"intervals.csv"),row.names=FALSE,fileEncoding="UTF-8")
  png(file.path(output,"intervals.png"),width=1280,height=768,res=160)
  plot_intervals()
  dev.off()
}
# 练习：①CSV改45/100、65/100；②恢复原数据，confidence改.90；③两组均0/20。
# 每次只改一项，从准备/读取开始按顺序运行，观察点估计、区间、图标题。
