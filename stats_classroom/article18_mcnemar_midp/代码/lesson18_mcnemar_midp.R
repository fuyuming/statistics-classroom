# 目标：沿用第11课58人例题，比较四P并计算条件第一类错误率。
# 编码UTF-8；从项目根目录运行；共享输入在数据/，结果在运行结果/对应语言。
# 依赖、安装、对象解释和预期结果见README；按四部分顺序执行。
# 统计课堂18：McNemar精确检验与mid-P；基础R；按四部分从上到下运行。
project_dir <- normalizePath(getwd(), mustWork=TRUE)
if (!dir.exists(file.path(project_dir, "数据"))) {
  project_dir <- normalizePath(file.path(project_dir, ".."), mustWork=TRUE)
}
data <- read.csv(file.path(project_dir,"数据","sle_mcnemar.csv"),stringsAsFactors=FALSE)
data
# a、b、c、d分别是荧光行与乳胶列的四个格。
a <- data$latex_positive[1]
b <- data$latex_negative[1]
c <- data$latex_positive[2]
d <- data$latex_negative[2]
N <- a+b+c+d
m <- b+c
data_summary <- data.frame(N=N,a=a,b=b,c=c,d=d,m=m)
data_summary

# 第一部分：H0下固定m后，b服从Binomial(m,.5)。
if (m == 0) {
  stop("m=0时停止检验，不生成新的检验结果。")
}
k_values <- 0:m
# dbinom给出P(B=k)：m是试验次数，0.5是H0下的方向概率。
probabilities <- dbinom(k_values,m,.5)
# T是距中心m/2的距离乘2，不改变极端程度的排序。
t_values <- abs(2*k_values-m)
t_observed <- abs(b-c)
# 方括号是逻辑筛选，sum把被选中的概率相加。
strict_probability <- sum(probabilities[t_values>t_observed])
boundary_probability <- sum(probabilities[t_values==t_observed])
exact_P <- strict_probability+boundary_probability
midP <- strict_probability+.5*boundary_probability

# 第二部分：普通Q、连续性校正Q及两种离散P。
ordinary_Q <- (b-c)^2/m
continuity_Q <- max(0,abs(b-c)-1)^2/m
# lower.tail=FALSE取卡方分布右尾；df=1指定自由度。
ordinary_P <- pchisq(ordinary_Q,df=1,lower.tail=FALSE)
continuity_P <- pchisq(continuity_Q,df=1,lower.tail=FALSE)
observed_four_methods <- data.frame(method=c("ordinary","continuity_correction",
    "exact","midP"),Q=c(ordinary_Q,continuity_Q,NA,NA),df=c(1,1,NA,
    NA),P=c(ordinary_P,continuity_P,exact_P,midP),strict_probability=c(NA,
    NA,strict_probability,strict_probability),boundary_probability=c(NA,
    NA,boundary_probability,boundary_probability))
boundary_components <- data.frame(b=b,c=c,m=m,t_observed=t_observed,
    strict_probability=strict_probability,boundary_probability=boundary_probability,
    exact_P=exact_P,midP=midP,binomtest_P=binom.test(b,m,.5)$p.value)

# 第三部分：函数输入四格数，逐步返回四种P值；m=0明确返回NA。
four_methods <- function(a_value,b_value,c_value,d_value,label) {
  N_value <- a_value+b_value+c_value+d_value
  m_value <- b_value+c_value
  if (m_value==0) { return(data.frame(exercise=label,a=a_value,b=b_value,
      c=c_value,d=d_value,N=N_value,ordinary_P=NA,continuity_P=NA,
      exact_P=NA,midP=NA)) }
  k <- 0:m_value
  p <- dbinom(k,m_value,.5)
  T <- abs(2*k-m_value)
  t <- abs(b_value-c_value)
  strict <- sum(p[T>t])
  boundary <- sum(p[T==t])
  Q <- (b_value-c_value)^2/m_value
  Qc <- max(0,abs(b_value-c_value)-1)^2/m_value
  ordinary <- pchisq(Q,1,lower.tail=FALSE)
  corrected <- pchisq(Qc,1,lower.tail=FALSE)
  return(data.frame(exercise=label,a=a_value,b=b_value,c=c_value,d=d_value,
      N=N_value,ordinary_P=ordinary,continuity_P=corrected,exact_P=strict+boundary,
      midP=strict+.5*boundary))
}
exercise_11_3 <- rbind(four_methods(a,b,c,d,"original_12_2"),four_methods(11,
    11,3,33,"teaching_11_3"))

# 第四部分：每个k都列出，按prob加权而不是除以行数；拒绝规则为P<=alpha。
all_outcomes <- data.frame(m=m,k=k_values,prob=probabilities,ordinary_P=NA,
    continuity_P=NA,exact_P=NA,midP=NA,ordinary_reject_alpha05=FALSE,
    continuity_reject_alpha05=FALSE,exact_reject_alpha05=FALSE,midP_reject_alpha05=FALSE)
for (i in seq_along(k_values)) {
  Q <- (2*k_values[i]-m)^2/m
  Qc <- max(0,abs(2*k_values[i]-m)-1)^2/m
  strict <- sum(probabilities[t_values>t_values[i]])
  boundary <- sum(probabilities[t_values==t_values[i]])
  all_outcomes$ordinary_P[i] <- pchisq(Q,1,lower.tail=FALSE)
  all_outcomes$continuity_P[i] <- pchisq(Qc,1,lower.tail=FALSE)
  all_outcomes$exact_P[i] <- strict+boundary
  all_outcomes$midP[i] <- strict+.5*boundary
  all_outcomes$ordinary_reject_alpha05[i] <- all_outcomes$ordinary_P[i]<=.05
  all_outcomes$continuity_reject_alpha05[i] <- all_outcomes$continuity_P[i]<=.05
  all_outcomes$exact_reject_alpha05[i] <- all_outcomes$exact_P[i]<=.05
  all_outcomes$midP_reject_alpha05[i] <- all_outcomes$midP[i]<=.05
}
conditional_size <- data.frame(m=m,alpha=.05,ordinary_size=sum(probabilities*all_outcomes$ordinary_reject_alpha05),
    continuity_size=sum(probabilities*all_outcomes$continuity_reject_alpha05),
    exact_size=sum(probabilities*all_outcomes$exact_reject_alpha05),
    midP_size=sum(probabilities*all_outcomes$midP_reject_alpha05))

# 人工练习：改alpha或m后重新枚举，并同时给出拒绝列。
alpha_exercise <- data.frame()
for (m_exercise in c(14,10)) {
  alpha_value <- if (m_exercise==14) .01 else .05
  k <- 0:m_exercise
  p <- dbinom(k,m_exercise,.5)
  T <- abs(2*k-m_exercise)
  for (i in seq_along(k)) {
    strict <- sum(p[T>T[i]])
    boundary <- sum(p[T==T[i]])
    Q <- (2*k[i]-m_exercise)^2/m_exercise
    Qc <- max(0,abs(2*k[i]-m_exercise)-1)^2/m_exercise
    exact_value <- strict+boundary
    mid_value <- strict+.5*boundary
    alpha_exercise <- rbind(alpha_exercise,data.frame(m=m_exercise,
        alpha=alpha_value,k=k[i],prob=p[i],ordinary_P=pchisq(Q,1,lower.tail=FALSE),
        continuity_P=pchisq(Qc,1,lower.tail=FALSE),exact_P=exact_value,
        midP=mid_value,ordinary_reject=pchisq(Q,1,lower.tail=FALSE)<=alpha_value,
        continuity_reject=pchisq(Qc,1,lower.tail=FALSE)<=alpha_value,
        exact_reject=exact_value<=alpha_value,midP_reject=mid_value<=alpha_value))
  }
}
alpha_exercise_summary <- data.frame()
for (m_value in unique(alpha_exercise$m)) {
  part <- alpha_exercise[alpha_exercise$m==m_value,]
  alpha_value <- unique(part$alpha)
  alpha_exercise_summary <- rbind(alpha_exercise_summary,data.frame(m=m_value,
      alpha=alpha_value,ordinary_size=sum(part$prob*part$ordinary_reject),
      continuity_size=sum(part$prob*part$continuity_reject),exact_size=sum(part$prob*part$exact_reject),
      midP_size=sum(part$prob*part$midP_reject)))
}
# 每个汇总行同时保存拒绝域中的k。
for (i in seq_len(nrow(alpha_exercise_summary))) {
  selected_m <- alpha_exercise_summary$m[i]
  part <- alpha_exercise[alpha_exercise$m == selected_m, ]
  for (method in c("ordinary", "continuity", "exact", "midP")) {
    reject_column <- paste0(method, "_reject")
    rejected_k <- part$k[part[[reject_column]]]
    alpha_exercise_summary[i, paste0(method, "_reject_k")] <- paste(rejected_k, collapse=" ")
  }
}
output_dir <- file.path(project_dir,"运行结果","R")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
write.csv(observed_four_methods,file.path(output_dir,"observed_four_methods.csv"),row.names=FALSE)
write.csv(boundary_components,file.path(output_dir,"boundary_components.csv"),row.names=FALSE)
write.csv(all_outcomes,file.path(output_dir,"all_outcomes.csv"),row.names=FALSE)
write.csv(conditional_size,file.path(output_dir,"conditional_size.csv"),row.names=FALSE)
write.csv(exercise_11_3,file.path(output_dir,"exercise_11_3.csv"),row.names=FALSE)
write.csv(alpha_exercise,file.path(output_dir,"alpha_exercise.csv"),row.names=FALSE)
write.csv(alpha_exercise_summary,file.path(output_dir,"alpha_exercise_summary.csv"),
    row.names=FALSE)
