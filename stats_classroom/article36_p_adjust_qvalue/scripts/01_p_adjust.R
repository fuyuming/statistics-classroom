# 统计课堂36：从入选门槛理解Bonferroni、Holm、BH和BY。
# 共享输入含1000个教学P值：800项无差异、200项有差异；标签只用于核对。
# is_null只用于说明模拟真相，调整算法只读取p列。
# 在RStudio打开本项目Rproj，从01.1依次运行；默认不保存文件。

# 01.1 读取完整比较族 ----
if (!file.exists("data/screening_pvalues.csv")) {
  stop("请先将工作目录设为本项目根目录，再从01.1开始。")
}
raw <- read.csv("data/screening_pvalues.csv")
head(raw)
p <- raw$p
m <- length(p)
alpha <- 0.05

# 01.2 手算Bonferroni与Holm ----
# 保存排序索引，最后用它将调整值放回原指标顺序。
order_p <- order(p)
sorted_p <- p[order_p]
rank_p <- seq_along(sorted_p)
bonf_manual <- pmin(1, m * p)
# Holm走到某项时，前面的每一步都必须通过，所以取前向累计最大值。
holm_candidate <- (m - rank_p + 1) * sorted_p
holm_sorted <- pmin(1, cummax(holm_candidate))
holm_manual <- numeric(m)
holm_manual[order_p] <- holm_sorted
head(data.frame(rank=rank_p, p=sorted_p, candidate=holm_candidate,
                adjusted=holm_sorted))

# 01.3 手算BH与BY ----
# 后面的某个位置若通过，前面的项也会入选，所以BH查看整个后段的最小值。
bh_candidate <- m * sorted_p / rank_p
# 先倒序，再cummin，最后倒回原来的升序位置。
bh_sorted <- pmin(1, rev(cummin(rev(bh_candidate))))
bh_manual <- numeric(m)
bh_manual[order_p] <- bh_sorted
harmonic_m <- sum(1 / rank_p)
by_manual <- pmin(1, harmonic_m * bh_manual)
head(data.frame(rank=rank_p, p=sorted_p, candidate=bh_candidate,
                adjusted=bh_sorted))

# 01.4 完整函数核对 ----
# stats::p.adjust以完整P值向量为输入，输出顺序与输入相同。
function_results <- data.frame(
  bonf=stats::p.adjust(p,"bonferroni"),
  holm=stats::p.adjust(p,"holm"),
  BH=stats::p.adjust(p,"BH"),
  BY=stats::p.adjust(p,"BY")
)
manual_results <- data.frame(bonf=bonf_manual, holm=holm_manual,
                             BH=bh_manual, BY=by_manual)
function_check <- data.frame(
  method=names(manual_results),
  # vapply逐方法返回一个最大绝对误差，结果组成数值列。
  max_abs_error=vapply(names(manual_results), function(name) {
    max(abs(manual_results[[name]]-function_results[[name]]))
  }, numeric(1))
)
function_check
stopifnot(all(function_check$max_abs_error < 1e-10))

# 01.5 两个各自单独构成比较族的三项Holm小例：显示停止与累计最大值 ----
# 每组三项单独构成m=3的比较族，不改变前面1000项的m。
mini_p_stop <- c(.020,.023,.040)
mini_p_increasing <- c(.020,.040,.100)
mini_rank <- 1:3
mini_stop_candidate <- (3-mini_rank+1) * mini_p_stop  # 乘数依次为3、2、1。
mini_stop_holm <- cummax(mini_stop_candidate)         # 前向累计最大值。
mini_increasing_candidate <- (3-mini_rank+1) * mini_p_increasing
mini_increasing_holm <- cummax(mini_increasing_candidate)
mini_example <- rbind(
  data.frame(example="stop", p=mini_p_stop,
             holm_candidate=mini_stop_candidate,
             holm=mini_stop_holm),
  data.frame(example="increasing", p=mini_p_increasing,
             holm_candidate=mini_increasing_candidate,
             holm=mini_increasing_holm)
)
print(mini_example)
cat("stop例第一关未通过，后续步骤停止；调整值仍返回完整三项。\n")
mini_expected <- c(.060,.060,.060,.060,.080,.100) # 仅作核对，不作为计算列。
stopifnot(identical(mini_example$example, c(rep("stop",3),rep("increasing",3))),
          max(abs(mini_example$holm - mini_expected)) < 1e-10,
          max(abs(mini_example$holm - c(p.adjust(mini_p_stop,"holm"),
                                         p.adjust(mini_p_increasing,"holm")))) < 1e-10)

# 01.6 原题三疗法复算：每个2x2比较使用Yates连续性校正 ----
# therapy前缀避免污染前面1000项的p、m和调整结果。
therapy_raw <- read.csv("data/therapies.csv")
therapy_table <- therapy_raw[, c("effective", "ineffective")]
therapy_pairs <- matrix(c(1, 2, 1, 3, 2, 3), ncol=2, byrow=TRUE)
therapy_names <- c("physical-drug", "physical-topical", "drug-topical")
therapy_raw_p <- numeric(nrow(therapy_pairs))
therapy_expected_min <- numeric(nrow(therapy_pairs))
for (therapy_i in seq_len(nrow(therapy_pairs))) {
  # correct=TRUE：明确开启2×2卡方检验的Yates连续性校正。
  therapy_fit <- chisq.test(therapy_table[therapy_pairs[therapy_i, ], , drop=FALSE],
                            correct=TRUE)
  therapy_raw_p[therapy_i] <- therapy_fit$p.value
  therapy_expected_min[therapy_i] <- min(therapy_fit$expected)
}
therapy_bonf <- pmin(1, length(therapy_raw_p) * therapy_raw_p)
therapy_order <- order(therapy_raw_p)
therapy_sorted <- therapy_raw_p[therapy_order]
therapy_holm_sorted <- pmin(1, cummax((length(therapy_raw_p) - seq_along(therapy_sorted) + 1) * therapy_sorted))
therapy_holm <- numeric(length(therapy_raw_p))
therapy_holm[therapy_order] <- therapy_holm_sorted
therapy_results <- data.frame(pair=therapy_names, raw_p_yates=therapy_raw_p,
                              bonferroni_p=therapy_bonf, holm_p=therapy_holm,
                              expected_min=therapy_expected_min)
therapy_results
stopifnot(max(abs(therapy_raw_p - c(.016754306570385496, .000009344615515985061,
                                     .047763774036015214))) < 1e-10,
          sum(therapy_bonf <= .05) == 1,
          all(therapy_holm <= .05))

# 01.7 查看完整结果和发现数 ----
results <- cbind(raw,manual_results)
summary_result <- data.frame(
  method=c("raw","bonf","holm","BH","BY"),
  discoveries=c(sum(p<=alpha),colSums(manual_results<=alpha))
)
summary_result
# 当前一批的发现数不是重复研究的FWER/FDR；这些错误率的定义见正文。

# 01.8 可选保存 ----
# 需要保存时先设置Sys.setenv(ARTICLE36_SAVE_RESULTS="1")，再顺序运行。
if (identical(Sys.getenv("ARTICLE36_SAVE_RESULTS"),"1")) {
  dir.create("results",showWarnings=FALSE)
  write.csv(results,"results/R_01_p_adjust.csv",row.names=FALSE)
  write.csv(summary_result,"results/R_01_summary.csv",row.names=FALSE)
  write.csv(mini_example,"results/R_01_mini_example.csv",row.names=FALSE)
  write.csv(therapy_results,"results/R_01_therapies.csv",row.names=FALSE)
}
