# 第一部分：看数据→明问题→立假设；显式导入，支持项目根或 code/ 工作目录。
options(stringsAsFactors = FALSE)
find_data <- function(name) {
  candidates <- c(file.path("data", name), file.path("..", "data", name))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(paste("找不到共享数据：", name))
  }
  hit[[1]]
}
textbook <- read.csv(find_data("textbook_example.csv"))
rd_data <- read.csv(find_data("independent_RD.csv"))
# 教材例5-8：甲、乙各5个独立小池，单位万尾；预定方向是甲高于乙。
# H0：两组总体分布相同；前提：两组独立、数据可排序；无并列值。
A <- textbook$value[textbook$group == "A"]
B <- textbook$value[textbook$group == "B"]
m <- length(A)
n <- length(B)
# 第二部分：选方法→查条件；合并编秩，解释原始秩和与U。
# 本例原始秩和R甲=37、R乙=18，扣除组内秩和15后U甲=22、U乙=3。
pooled_rank <- rank(c(A, B), ties.method = "average")
RA <- sum(pooled_rank[seq_len(m)])
RB <- sum(pooled_rank[m + seq_len(n)])
UA <- RA - m * (m + 1) / 2
UB <- RB - n * (n + 1) / 2
win_A <- sum(outer(A, B, ">"))
win_B <- sum(outer(A, B, "<"))
win_tie <- sum(outer(A, B, "=="))
effect_A <- UA / (m * n)
textbook_summary <- data.frame(n_A=m,n_B=n,R_A=RA,R_B=RB,U_A=UA,U_B=UB,wins_A=win_A,wins_B=win_B,ties=win_tie,effect_A=effect_A)
textbook_summary
# 第三部分：做检验→判显著→看效应；枚举10选5的252个标签分配，U甲≥22有7种。
all_values <- c(A, B)
N <- length(all_values)
assignment <- combn(seq_len(N), m)
u_all <- apply(assignment, 2, function(index) {
  sum(index) - m * (m + 1) / 2
})
count_assignments <- length(u_all)
p_upper <- mean(u_all >= UA)
p_lower <- mean(u_all <= UA)
p_two_sided <- min(1, 2 * min(p_upper, p_lower))
u_critical_lower <- 4
p_critical_lower <- mean(u_all <= u_critical_lower)
lower_count <- sum(u_all <= u_critical_lower)
# HL=5；共同形状位置平移模型下，U临界值4对应第(u临界值+1)=5个有序差值2，实际覆盖20/21=.952381。
cross_diff <- sort(as.vector(outer(A, B, "-")))
HL <- median(cross_diff)
lower_index <- u_critical_lower + 1
lower_95_nominal <- cross_diff[lower_index]
actual_confidence <- 1 - p_critical_lower
exact_summary <- data.frame(assignments=count_assignments,upper_count=sum(u_all >= UA),p_upper=p_upper,p_lower=p_lower,p_two_sided=p_two_sided,u_critical_lower=u_critical_lower,lower_count=lower_count,p_critical_lower=p_critical_lower,HL=HL,ci_lower=lower_95_nominal,ci_upper=Inf,actual_confidence=actual_confidence)
exact_summary
# R的W=22是甲组U型统计量；conf.int给位置差区间，需共同形状位置模型。
r_result <- wilcox.test(A,B,alternative="greater",exact=TRUE,conf.int=TRUE)
r_result
r_ci <- unname(r_result$conf.int)
# 交换组回归检查：交换输入后改用B<A的等价方向，单侧P应为0.0277778。
swap_result <- wilcox.test(B,A,alternative="less",exact=TRUE)
swap_u <- unname(swap_result$statistic)
swap_p_upper <- mean(u_all >= swap_u)
swap_p_lower <- mean(u_all <= swap_u)
swap_p_two_sided <- min(1, 2 * min(swap_p_upper, swap_p_lower))
swap_summary <- data.frame(swapped_U=swap_u,swapped_p_less=swap_result$p.value,swapped_p_two_sided=swap_p_two_sided)
swap_summary
# 第四部分：RD课件对照；两组n=10、12，有并列值，秩取平均，使用近似+连续性校正。
RD1 <- rd_data$RD[rd_data$group == 1]
RD2 <- rd_data$RD[rd_data$group == 2]
rd_rank <- rank(c(RD1, RD2), ties.method = "average")
rd_n1 <- length(RD1)
rd_n2 <- length(RD2)
RA_RD <- sum(rd_rank[seq_len(rd_n1)])
RB_RD <- sum(rd_rank[rd_n1 + seq_len(rd_n2)])
UA_RD <- RA_RD - rd_n1 * (rd_n1 + 1) / 2
UB_RD <- RB_RD - rd_n2 * (rd_n2 + 1) / 2
rd_win_tie <- sum(outer(RD1, RD2, "=="))
rd_effect <- UA_RD / (rd_n1 * rd_n2)
rd_result <- wilcox.test(RD1,RD2,alternative="greater",exact=FALSE,correct=TRUE)
rd_summary <- data.frame(n_1=rd_n1,n_2=rd_n2,R_1=RA_RD,U_1=UA_RD,U_2=UB_RD,ties=rd_win_tie,effect_1=rd_effect,p_approx_continuity=rd_result$p.value)
rd_summary
# 第五部分：配对资料；同一问题可写成差值的单样本符号秩或paired调用。
# 来源为颜艳、王彤教材课件的8对大鼠维生素A数据，X与Y按pair对应。
# H0：差值分布关于0对称；d=round(X-Y,2)，零差剔除；本例n=8、无绝对差并列。
paired_data <- read.csv(find_data("paired_vitaminA.csv"))
X <- paired_data$X
Y <- paired_data$Y
d <- round(X - Y, 2)
d_nonzero <- d[d != 0]
d_rank <- rank(abs(d_nonzero), ties.method = "average")
positive_rank_sum <- sum(d_rank[d_nonzero > 0])
negative_rank_sum <- sum(d_rank[d_nonzero < 0])
signed_rank_summary <- data.frame(n=length(d_nonzero),positive_rank_sum=positive_rank_sum,negative_rank_sum=negative_rank_sum)
signed_rank_summary
one_sample_result <- wilcox.test(d,mu=0,alternative="two.sided",exact=TRUE,conf.int=TRUE)
paired_result <- wilcox.test(X,Y,paired=TRUE,alternative="two.sided",exact=TRUE,conf.int=TRUE)
one_sample_result
paired_result
# 可选第六部分：保存可分享的关键结果CSV；默认不写文件。
result_export <- data.frame(language="R",R_A=RA,R_B=RB,U_A=UA,U_B=UB,p_upper=p_upper,p_two_sided=p_two_sided,HL=HL,ci_lower=r_ci[1],ci_upper=r_ci[2],actual_confidence=actual_confidence,RD_R1=RA_RD,RD_U1=UA_RD,RD_U2=UB_RD,RD_effect=rd_effect,RD_p=rd_result$p.value)
save_csv <- FALSE
if (save_csv) {
  project_root <- dirname(dirname(normalizePath(find_data("textbook_example.csv"))))
  output_dir <- file.path(project_root, "results")
  dir.create(output_dir, showWarnings = FALSE)
  write.csv(result_export, file.path(output_dir, "R_results.csv"), row.names = FALSE)
}
