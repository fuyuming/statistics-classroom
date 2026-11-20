# 问题与输入 ----
# 假设100名不同患者各有一段心电。两算法是否标记至少一次室早，构成一对结果。
# 目标比较B−A总体阳性判定率；无金标准，不比较准确性。患者间独立、无缺失。
tab <- as.matrix(read.csv("data/ecg_2x2.csv", row.names=1))
tab
# 两行是A阳性/阴性，两列是B阳性/阴性；每格是联合类别的人数，不是一位患者。
n <- sum(tab)                       # 全部配对数，率差的分母
b <- tab[1,2]                       # A阳性、B阴性：个人差值−1
c_count <- tab[2,1]                 # A阴性、B阳性：个人差值＋1
m <- b + c_count                    # 不一致对子数；本例30，供条件二项法使用

# 手动计算：率差及普通配对Wald区间 ----
z <- c(rep(-1,b), rep(0,n-m), rep(1,c_count))
RD <- (c_count-b)/n
se_RD <- sd(z)/sqrt(n)              # 样本差值SD分母n−1；SE用全部配对数
RD_CI <- RD + qnorm(c(.025,.975))*se_RD
# 使用标准正态临界值、无连续性校正；不是下述校正检验的区间反演。

# 手动计算：连续性校正McNemar及条件二项精确P ----
Q <- max(abs(c_count-b)-1,0)^2/m
p_chi <- pchisq(Q, 1, lower.tail=FALSE) # χ²自由度1的右尾概率
# H0下，不一致对中的两个方向各占1/2；此对称二项情形可取较小一侧尾概率的2倍。
p_exact <- min(1, 2*pbinom(min(b,c_count),m,.5))
result <- data.frame(n=n,b=b,c=c_count,m=m,Q=Q,p=p_chi,exact=p_exact,
                     RD=RD,se=se_RD,low=RD_CI[1],high=RD_CI[2])
result

# 调用函数：基础R的McNemar和二项检验 ----
fit <- mcnemar.test(tab, correct=TRUE)
fit
fit_exact <- binom.test(min(b,c_count), m, p=.5, alternative="two.sided")
fit_exact$p.value
# 这里仅提取二项P；该函数的条件比例区间不是总体率差RD的区间。
# 结果：RD=.20，校正Q≈12.0333、P≈.0005226；精确P≈.0003249。
# RD区间约[.09956,.30044]，乘100后为百分点；不一致m不是率差分母。
write.csv(result, "results/R_binary.csv", row.names=FALSE)
