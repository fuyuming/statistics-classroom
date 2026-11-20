# 统计课堂34：同一张配对表的率差与Cohen κ ----
# 例题：58名可疑系统性红斑狼疮患者，各接受免疫荧光法A与乳胶法B检测抗核抗体。
# 问：B−A阳性率差多少？判定有多一致？总体κ的95%区间怎样算？
# 数据：颜艳、王彤《医学统计学》第5版例7-3、表7-3，教材联合汇总计数。
# 第11课PDF P16及P19–25为本例课堂入口；一般κ区间及固定边际界限为本文扩展。
# 条件：不同患者联合记录独立，每人一对完整的阳/阴结果；本表按A行、B列。
# RStudio打开本项目Rproj，再从第一节按顺序运行；扩展包vcd提供完整κ/区间接口。

# 输入：显式规定行列类别顺序 ----
tab_raw <- read.csv("data/paired_tests.csv", stringsAsFactors=FALSE)
tab_raw
# factor的levels保证两轴都按阳性、阴性排列，避免按字母排序改变表的读法。
tab_raw$A <- factor(tab_raw$A, levels=c("positive", "negative"))
tab_raw$B <- factor(tab_raw$B, levels=c("positive", "negative"))
tab <- xtabs(count ~ A + B, data=tab_raw)
tab
n <- sum(tab)                         # 全部58位患者，每位恰好落入一格
a <- tab["positive", "positive"]
b <- tab["positive", "negative"]
c_count <- tab["negative", "positive"]
d <- tab["negative", "negative"]

# 手动计算：率差、观察一致率、机会一致率及κ ----
p <- tab / n                         # 四格的经验联合概率
r <- rowSums(p)                       # A的两类边际比例
cprob <- colSums(p)                   # B的两类边际比例
pA <- r["positive"]
pB <- cprob["positive"]
rd <- pB - pA                        # B−A的概率差；乘100变为百分点
po <- sum(diag(p))                    # 对角线上为判定相同者
pe <- sum(r * cprob)                  # 两法各类边际比例乘积相加，构造机会基准
kappa <- (po - pe) / (1 - pe)
point_estimates <- data.frame(n=n, pA=pA, pB=pB, rd=rd, po=po, pe=pe, kappa=kappa)
point_estimates
# rd≈−0.1724，表示B阳性率低17.24个百分点；Po≈.7586，κ≈.4550均无量纲。

# 手动计算：一般多项delta标准误与95%Wald区间 ----
# g[i,j]是κ对第i行第j列比例的局部变化系数；按四格循环计算完整公式。
g <- matrix(0, nrow=2, ncol=2)
for (i in 1:2) {
  for (j in 1:2) {
    diagonal <- as.numeric(i == j)     # 同类别的对角格使观察一致率增加
    pe_sensitivity <- r[j] + cprob[i] # 机会一致率对该格比例的局部变化系数
    g[i,j] <- (diagonal - (1-kappa)*pe_sensitivity) / (1-pe)
  }
}
g
# 四格比例合计1，格子之间存在协方差。不是四份互相独立的样本。
g_mean <- sum(p * g)
g_second <- sum(p * g^2)
variance <- (g_second - g_mean^2) / n  # 四值辅助量方差除以独立患者数
se <- sqrt(variance)
# 双侧95%区间两端各留2.5%，qnorm把累计概率换成标准正态分位点。
ci <- kappa + qnorm(c(.025, .975)) * se
interval_result <- data.frame(g_mean=g_mean, g_second=g_second, variance=variance,
                              se=se, ci_low=ci[1], ci_high=ci[2])
interval_result
# 一般SE≈.115335，95%Wald近似区间≈[.228980,.681087]。

# 调用函数：vcd::Kappa与stats::confint ----
# vcd是扩展包，首次使用先按README安装；完整函数以整张方形计数表为输入。
vk <- vcd::Kappa(tab)
vcd_estimate <- vk$Unweighted          # value为未加权κ，ASE为一般渐近标准误
vcd_interval <- confint(vk, level=.95)["Unweighted", ]
vcd_estimate
vcd_interval

# 结果核对：手算κ、一般SE、区间对应完整函数字段 ----
function_comparison <- data.frame(
  quantity=c("kappa", "general_SE", "CI_low", "CI_high"),
  manual=c(kappa, se, ci[1], ci[2]),
  vcd=c(vcd_estimate[["value"]], vcd_estimate[["ASE"]], vcd_interval[["lwr"]], vcd_interval[["upr"]]))
function_comparison
# 本节取value/ASE与区间；κ=0的专用SE及右尾P在下一节按Rosner计算。

# 手动计算：κ=0对预先指定κ>0的右尾近似检验 ----
# SE0是在零假设下计算的渐近标准误，和区间使用的一般SE服务于不同任务。
null_numerator <- pe + pe^2 - sum(r * cprob * (r + cprob))
se0 <- sqrt(null_numerator / (n * (1-pe)^2))
z <- kappa / se0
p_greater <- pnorm(z, lower.tail=FALSE) # 事先问“高于”，所以取标准正态右尾
null_test <- data.frame(se0=se0, z=z, p_greater=p_greater)
null_test
# SE0≈.120950，Z≈3.7622，单侧P≈.00008422，支持高于机会基准。

# 手动计算：McNemar条件二项双侧精确P ----
m <- b + c_count                      # 14个不一致对子；率差的分母仍是n=58
# H0下两个分歧方向等可能；本题二项分布对称，较小尾概率乘2，最大不超过1。
mcnemar_exact_p <- min(1, 2 * pbinom(min(b, c_count), size=m, prob=.5))
mcnemar_exact_p

# 调用函数：stats::binom.test的条件精确二项接口 ----
exact_fit <- binom.test(x=c_count, n=m, p=.5, alternative="two.sided")
exact_fit
# p.value对应精确McNemar的P≈.01294；此函数自带conf.int针对方向概率，不是总体率差。

# 结果核对：双侧精确P与函数字段 ----
exact_comparison <- data.frame(manual=mcnemar_exact_p, function_P=exact_fit$p.value)
exact_comparison

# 手动枚举：固定样本边际，改变共同阳性人数 ----
a_pos <- sum(tab["positive", ])
b_pos <- sum(tab[, "positive"])
a_min <- max(0, a_pos+b_pos-n)         # 保证共同阴性人数也非负
a_max <- min(a_pos, b_pos)
# 每次循环建立一张满足相同边际的构造表；教材原表在a=11这一行。
bounds_list <- list()
for (aa in seq.int(a_min, a_max)) {
  bb <- a_pos-aa
  cc <- b_pos-aa
  dd <- n-aa-bb-cc
  poo <- (aa+dd)/n
  # 边际保持不变，故pe及B−A率差保持不变。
  one_row <- data.frame(a=aa,b=bb,c=cc,d=dd,rd=(cc-bb)/n,po=poo,kappa=(poo-pe)/(1-pe))
  bounds_list[[length(bounds_list)+1]] <- one_row
}
bounds <- do.call(rbind, bounds_list)
bounds
kappa_min <- min(bounds$kappa)
kappa_max <- max(bounds$kappa)
po_min <- min(bounds$po)
po_max <- max(bounds$po)
# 该范围固定当前观察边际；总体κ区间还考虑边际估计的不确定性，所以不按此范围截断。

# 汇总：查看方向、单位与两类推断结果 ----
summary <- data.frame(n=n,pA=pA,pB=pB,rd=rd,po=po,pe=pe,kappa=kappa,
  se0=se0,z=z,p_greater=p_greater,se=se,ci_low=ci[1],ci_high=ci[2],
  mcnemar_exact_p=mcnemar_exact_p,kappa_min=kappa_min,kappa_max=kappa_max,po_min=po_min,po_max=po_max)
summary

# 思考：先预测，再改动副本中的联合计数 ----
# 将a改为13，同时按固定边际重算b/c/d；预测率差、Po和κ各怎样变化。
# κ≈.455按Landis–Koch六级参考为moderate，按Rosner三级为good；分级不改变数值或区间。

# 可选保存：首次练习默认不写文件 ----
# 若需要保存，在Console先运行Sys.setenv(ARTICLE34_SAVE_RESULTS="1")，再运行本节。
if (identical(Sys.getenv("ARTICLE34_SAVE_RESULTS"), "1")) {
  dir.create("results", showWarnings=FALSE)
  write.csv(summary,"results/R_summary.csv",row.names=FALSE)
  write.csv(bounds,"results/R_bounds.csv",row.names=FALSE)
}
