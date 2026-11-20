# 统计课堂34：Cohen κ的逐步W表推导 ----
# 2020年颜艳、王彤《医学统计学》第5版例7-3、表7-3，pp.110–111：抗核抗体检测。
# A为免疫荧光法，B为乳胶凝集法；四格计数依次为11、12、2、33。
# 目标：从三个0/1均数计算κ的方差、SE与1.96 Wald置信区间。
# 同一患者的X/Y/C有关联，但不同患者相互独立；联合四格表保留这种联系。
# 运行前将工作目录设为本项目根目录；默认不保存，设置ARTICLE34_SAVE_RESULTS=1才写结果。

# 03.1 输入与完整病例 ----
# X=A阳性，Y=B阳性，C=两法结果一致；每位患者只有一条完整配对记录。
tab_raw <- read.csv("data/paired_tests.csv", stringsAsFactors=FALSE)
tab_raw$A <- factor(tab_raw$A, levels=c("positive", "negative"))
tab_raw$B <- factor(tab_raw$B, levels=c("positive", "negative"))
tab <- xtabs(count ~ A + B, data=tab_raw)
tab
n <- sum(tab)                         # 58名患者，每位恰好落入四格之一

# 03.2 按四格建立X/Y/C均数 ----
# 理论p=EX、q=EY、r=EC未知；本题以三个0/1记录的样本均数作帽号代入。
p_hat <- (tab["positive", "positive"] + tab["positive", "negative"]) / n
q_hat <- (tab["positive", "positive"] + tab["negative", "positive"]) / n
r_hat <- (tab["positive", "positive"] + tab["negative", "negative"]) / n
means <- data.frame(p_hat=p_hat, q_hat=q_hat, r_hat=r_hat)
means
# p_hat、q_hat、r_hat分别是A阳性、B阳性、两法一致的样本均数，均来自同一批患者。

# 03.3 计算三个系数 ----
# 机会一致率e=pq+(1-p)(1-q)，用帽号值代入；κ=(r-e)/(1-e)。
e_hat <- p_hat*q_hat + (1-p_hat)*(1-q_hat)
kappa <- (r_hat-e_hat)/(1-e_hat)
# cr、cp、cq是κ对C、X、Y的局部一阶系数；按样本代入而非假定三变量独立。
cr <- 1/(1-e_hat)
cp <- -(1-kappa)*(2*q_hat-1)/(1-e_hat)
cq <- -(1-kappa)*(2*p_hat-1)/(1-e_hat)
coefficients <- data.frame(e_hat=e_hat, kappa=kappa, cr=cr, cp=cp, cq=cq)
coefficients

# 03.4 联合W表和两种矩 ----
# 四格顺序是(A+,B+)、(A+,B−)、(A−,B+)、(A−,B−)，对应W=cr*C+cp*X+cq*Y。
W <- matrix(c(cr+cp+cq, cp, cq, cr), nrow=2, byrow=TRUE,
            dimnames=dimnames(tab))
W
counts <- tab
m1 <- sum(counts * W) / n               # E(W)的经验概率矩，分母是全部患者n
m2 <- sum(counts * W^2) / n             # E(W^2)，先平方W再按四格加权
moments <- data.frame(m1=m1, m2=m2)
moments
# 四格W表保留X/Y/C在同一人的联系；不要把三个边际当成三份独立样本。

# 03.5 方差、SE与1.96 Wald区间 ----
# 经验概率矩给出Var(W)=E(W²)-[E(W)]²；按独立患者数n得到κ的一阶方差。
variance <- (m2-m1^2) / n
se <- sqrt(variance)
ci <- kappa + c(-1, 1) * 1.96 * se
stepwise_result <- data.frame(p_hat=p_hat, q_hat=q_hat, r_hat=r_hat, e_hat=e_hat, kappa=kappa,
  cr=cr, cp=cp, cq=cq, Wa=W[1,1], Wb=W[1,2], Wc=W[2,1], Wd=W[2,2],
  m1=m1, m2=m2, variance=variance, se=se, ci_low=ci[1], ci_high=ci[2])
stepwise_result
# 本例variance≈0.0133022442、SE≈0.1153353554；Wald区间是κ±1.96SE。
# 这是一阶代入近似，不宣称样本κ估计量严格无偏。

# 03.6 完整函数核对与可选保存 ----
# vcd::Kappa以完整2×2计数表计算未加权κ及渐近SE；其实现用于独立核对点估计和SE。
vk <- vcd::Kappa(tab)
vcd_check <- data.frame(kappa_manual=kappa,
  kappa_vcd=vk$Unweighted[["value"]], se_manual=se,
  se_vcd=vk$Unweighted[["ASE"]])
vcd_check
# 可选保存：默认不写文件；只有显式设置环境变量才输出
if (identical(Sys.getenv("ARTICLE34_SAVE_RESULTS"), "1")) {
  dir.create("results", showWarnings=FALSE)
  write.csv(stepwise_result, "results/R_stepwise.csv", row.names=FALSE)
}
