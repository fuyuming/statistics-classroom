# 第10课例8-3：肺癌组RD是否倾向于高于硅沉着病0期工人组？
# 本脚本以本目录为工作目录运行；不使用cat或print解释结果。
options(width = 42)
rd <- read.csv("../data/independent_RD.csv")
A <- rd$RD[rd$group == 1]
B <- rd$RD[rd$group == 2]

# 先将22个观测合并编秩，再按原组别分别求和；3.23并列，取平均秩2.5。
pooled_rank <- rank(rd$RD, ties.method = "average")
rank_table <- data.frame(group = rd$group, RD = rd$RD, rank = pooled_rank)
rank_table
R_A <- sum(pooled_rank[rd$group == 1])
R_B <- sum(pooled_rank[rd$group == 2])
U_A <- R_A - length(A) * (length(A) + 1) / 2
U_B <- R_B - length(B) * (length(B) + 1) / 2
c(R_A = R_A, R_B = R_B, U_A = U_A, U_B = U_B)
# R_A=141.5，R_B=111.5；U_A=86.5，U_B=33.5。

# 第一组是肺癌组；greater对应预先提出的较大方向。
# 沿用课堂的正态近似；软件修正并列值方差，correct=TRUE作连续性校正。
rd_result <- wilcox.test(A, B, alternative = "greater",
                        exact = FALSE, correct = TRUE)
rd_result
