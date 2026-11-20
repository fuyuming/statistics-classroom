# 问题与输入 ----
# 假设8名不同患者各提供一段30秒心电；同一行的A/B读数来自同一记录。
# 单位为次/分，目标为总体平均差B−A；不同对子独立，按差值近似正态演算。
# 从项目根目录运行；本例8个ID唯一、无缺失，n是完整配对数。
rates <- read.csv("data/ecg_rates.csv")
rates
A <- rates$A_rate
B <- rates$B_rate
d <- B - A

# 手动计算：差值、标准误与配对t ----
n <- length(d)
dbar <- mean(d)
sd_d <- sd(d)                       # 样本SD使用n−1作分母
se <- sd_d / sqrt(n)                # 平均差的标准误，仍为次/分
df <- n - 1
t_value <- dbar / se                # H0：总体平均差为0
p_two <- 2 * pt(-abs(t_value), df)   # 对称t分布：左尾乘2得到双侧P
ci <- dbar + qt(c(.025,.975), df) * se
result <- data.frame(n=n, mean_diff=dbar, sd_diff=sd_d, se=se,
                     t=t_value, df=df, p=p_two, low=ci[1], high=ci[2])
result

# 调用函数：基础R的配对t检验 ----
fit <- t.test(B, A, paired=TRUE, conf.level=.95)
fit
# 对应检查：均值差−2、t≈−6.1101、df=7、双侧P≈.0004862，区间单位次/分。
# 这里只比较平均读数，没有金标准，不能判断哪种算法更准确。

# 练习：对应关系影响差值的波动 ----
# 仅在练习副本倒序B：B列均值不变，差值SD变大；实测数据必须按真实ID配对。
data.frame(mean_B=mean(B), mean_B_reordered=mean(rev(B)),
           sd_diff=sd(d), sd_diff_reordered=sd(rev(B)-A))
write.csv(result, "results/R_paired.csv", row.names=FALSE)
