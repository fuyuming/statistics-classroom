# 统计课堂37：区分完整构成、指定类别比例和稀疏表检验（第11课教学例题计数，非真实调查）。
# faults为每组80台已失效记录器的故障构成；sparse为每组10台的合格/不合格计数。
# 依赖 effectsize/rcompanion/rstatix/report；UTF-8保存。打开Rproj后按九节运行，默认显示结果。
# 设 ARTICLE37_SAVE_RESULTS=1 后第9节才保存四张表到 results/R。
# 1. 数据与参数 ----
needed <- c("effectsize", "rcompanion", "rstatix", "report")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("请先安装依赖包：install.packages('", pkg, "')", sep = ""))
  }
}
if (!file.exists("data/faults.csv")) {
  stop("请先打开本项目的 Rproj，确认工作目录包含 data 文件夹。")
}
faults <- read.csv("data/faults.csv")
sparse <- read.csv("data/sparse.csv")
# 先看字段和全部三行：faults每行一组、分母为80台已失效记录器；sparse每行一组、分母为10台。
print(names(faults))
print(faults)
print(names(sparse))
print(sparse)
# 后续矩阵保留计数列，分别回答完整构成、单类比例和稀疏表检验问题。
tab <- as.matrix(faults[c("open", "electrode", "board")])
rownames(tab) <- faults$row
names(dimnames(tab)) <- c("packaging", "fault")
sparse_tab <- as.matrix(sparse[c("pass", "fail")])
rownames(sparse_tab) <- sparse$row
conf_level <- 0.95
alpha <- 0.05
# 每列为一个组对：第1列F-G，第2列F-H，第3列G-H。
pair_index <- combn(seq_len(nrow(tab)), 2)

# 2. 整体检验与关联效应量 ----
# 问题：三组的完整故障构成是否相同？同时报告整体Q、P和V。
fit <- stats::chisq.test(tab, correct = FALSE)
# fit是原生chisq.test返回的检验对象；adjust=FALSE保留未修正课堂V。
v_all <- as.data.frame(effectsize::cramers_v(
  tab, adjust = FALSE, ci = conf_level, alternative = "two.sided"
))
main_overall <- data.frame(
  q = unname(fit$statistic), df = unname(fit$parameter),
  p = fit$p.value, n = sum(tab), V = v_all$Cramers_v,
  V_low = v_all$CI_low, V_high = v_all$CI_high
)
print(main_overall)
cat("Q/P回答整体差异，V及区间概括关联强度；这里尚未进行组对Holm调整。\n")
# report用于组织本次整体检验的文字结果。
overall_report <- report::report(
  fit, adjust = FALSE, ci = conf_level, alternative = "two.sided"
)
print(overall_report)
# 实跑：Q=44.79167，P=4.392904e-09，V=0.3054766。

# 3. 完整分类构成的三项组对 ----
# 思路：每次保留三类故障，重算2×3子表，再对三对原始P作Holm调整。
full_list <- vector("list", ncol(pair_index))
for (k in seq_len(ncol(pair_index))) {
  index <- pair_index[, k]
  sub_tab <- tab[index, , drop = FALSE] # index按行取一对组，保留三类故障
  sub_fit <- stats::chisq.test(sub_tab, correct = FALSE)
  v <- as.data.frame(effectsize::cramers_v(
    sub_tab, adjust = FALSE, ci = conf_level, alternative = "two.sided"
  ))
  full_list[[k]] <- data.frame(
    pair = paste(rownames(tab)[index], collapse = "-"),
    q = unname(sub_fit$statistic), df = unname(sub_fit$parameter),
    p_raw = sub_fit$p.value, min_expected = min(sub_fit$expected),
    V = v$Cramers_v[1], V_low = v$CI_low[1], V_high = v$CI_high[1]
  )
}
full_pairs <- do.call(rbind, full_list)
full_pairs$p_holm <- stats::p.adjust(full_pairs$p_raw, method = "holm")
detected_full <- full_pairs$pair[full_pairs$p_holm < alpha]
cat("完整构成比较（alpha=", alpha, "）：\n", sep = "")
if (length(detected_full) > 0) {
  cat("检出差异的组对：", paste(detected_full, collapse = ", "), "\n")
} else {
  cat("无检出差异的组对。\n")
}
full_pairs <- full_pairs[c("pair", "q", "df", "p_raw", "p_holm",
                           "min_expected", "V", "V_low", "V_high")]
# 正文使用的包入口；与逐对子表计算核对，stats=FALSE只保留所选P值。
package_full <- rcompanion::pairwiseNominalIndependence(
  tab, compare = "row", fisher = FALSE, gtest = FALSE, chisq = TRUE,
  yates = FALSE, method = "holm", digits = 15, stats = FALSE
)
stopifnot(max(abs(package_full$p.adj.Chisq - full_pairs$p_holm)) < 1e-10)
print(full_pairs)
# p_raw是单个2×3检验，p_holm是F-G、F-H、G-H这一比较族的调整P。
# 实跑：三个Holm P约0.002365、4.298e-09、0.004714。

# 4. 单类比例：电极脱落与其他 ----
# 问题：只比较电极脱落比例时，把其余两类合并为other。
electrode_tab <- cbind(
  electrode = tab[, "electrode"],
  # rowSums按组求合计；总数减电极脱落得到other。
  other = rowSums(tab) - tab[, "electrode"]
)
category_list <- vector("list", ncol(pair_index))
for (k in seq_len(ncol(pair_index))) {
  index <- pair_index[, k]
  x <- electrode_tab[index, "electrode"]
  n_group <- rowSums(electrode_tab[index, , drop = FALSE])
  # 原生prop.test返回P值和区间；estimate/conf.int按输入顺序，差为前组减后组。
  rate_fit <- stats::prop.test(x, n_group, correct = TRUE, conf.level = conf_level)
  category_list[[k]] <- data.frame(
    pair = paste(rownames(tab)[index], collapse = "-"),
    q_yates = unname(rate_fit$statistic), df = unname(rate_fit$parameter),
    p_raw = rate_fit$p.value, p1 = x[1] / n_group[1], p2 = x[2] / n_group[2],
    diff = unname(rate_fit$estimate[1] - rate_fit$estimate[2]),
    diff_low = rate_fit$conf.int[1], diff_high = rate_fit$conf.int[2]
  )
}
category_pairs <- do.call(rbind, category_list)
category_pairs$p_holm <- stats::p.adjust(category_pairs$p_raw, method = "holm")
detected_category <- category_pairs$pair[category_pairs$p_holm < alpha]
cat("电极脱落比例比较（alpha=", alpha, "）：\n", sep = "")
if (length(detected_category) > 0) {
  cat("检出差异的组对：", paste(detected_category, collapse = ", "), "\n")
} else {
  cat("无检出差异的组对。\n")
}
category_pairs <- category_pairs[c("pair", "q_yates", "df", "p_raw", "p_holm",
                                   "p1", "p2", "diff", "diff_low", "diff_high")]
package_category <- rstatix::pairwise_prop_test(
  electrode_tab, p.adjust.method = "holm", correct = TRUE
)
stopifnot(max(abs(package_category$p.adj - category_pairs$p_holm)) < 1e-10)
print(category_pairs)
cat("p_holm仍按三对调整；p1/p2为前后组比例，diff及区间均为前组减后组。\n")
# 实跑：三个Holm P为1、1、0.4249341；G-H差0.125，CI约[-0.03627,0.28627]。

# 5. 稀疏表的Fisher组对 ----
# 思路：每对保留合格/不合格两列；小预期频数时用双侧Fisher，再调整三对P。
sparse_list <- vector("list", ncol(pair_index))
for (k in seq_len(ncol(pair_index))) {
  index <- pair_index[, k]
  sub_tab <- sparse_tab[index, , drop = FALSE]
  # 原生fisher.test返回条件OR estimate及conf.int；alternative指定双侧。
  exact_fit <- stats::fisher.test(sub_tab, alternative = "two.sided",
                                  conf.level = conf_level)
  # 用行列合计另算预期数作检查，不用卡方近似计算本节P值。
  expected <- outer(rowSums(sub_tab), colSums(sub_tab)) / sum(sub_tab)
  sparse_list[[k]] <- data.frame(
    pair = paste(rownames(sparse_tab)[index], collapse = "-"),
    p_raw = exact_fit$p.value, min_expected = min(expected),
    OR = unname(exact_fit$estimate), OR_low = exact_fit$conf.int[1],
    OR_high = exact_fit$conf.int[2]
  )
}
sparse_pairs <- do.call(rbind, sparse_list)
sparse_pairs$p_holm <- stats::p.adjust(sparse_pairs$p_raw, method = "holm")
detected_sparse <- sparse_pairs$pair[sparse_pairs$p_holm < alpha]
cat("小样本合格率比较（alpha=", alpha, "）：\n", sep = "")
if (length(detected_sparse) > 0) {
  cat("检出差异的组对：", paste(detected_sparse, collapse = ", "), "\n")
} else {
  cat("无检出差异的组对。\n")
}
sparse_pairs <- sparse_pairs[c("pair", "p_raw", "p_holm", "min_expected",
                               "OR", "OR_low", "OR_high")]
package_sparse <- rstatix::pairwise_fisher_test(
  sparse_tab, p.adjust.method = "holm", alternative = "two.sided"
)
stopifnot(max(abs(package_sparse$p.adj - sparse_pairs$p_holm)) < 1e-10)
print(sparse_pairs)
# p_raw来自精确检验，p_holm是三对合格率比较族的调整P；min_expected只作稀疏性检查。
# 实跑：最小E为1.5、4、4.5；Holm P为1、0.05929983、0.1395570。

# 6. 两类问题的调整P并列 ----
# 问题：同一组对在完整构成与指定类别比例问题下，结论是否相同？
comparison <- data.frame(
  pair = full_pairs$pair,
  full_p_holm = full_pairs$p_holm,
  category_p_holm = category_pairs$p_holm
)
print(comparison)
cat("两列Holm P分别属于完整构成三对族与指定类别比例三对族。\n")

# 7. 独立检验练习：Hochberg ----
# 先用两项独立检验的模型演算Hochberg的FWER：各项H0为真时P~Uniform(0,1)。
alpha_proof <- 0.05
half_alpha <- alpha_proof / 2
prob_any_small <- 1 - (1 - half_alpha)^2
prob_both_middle <- (alpha_proof - half_alpha)^2
fwer_two_null <- prob_any_small + prob_both_middle
event_table <- data.frame(
  event = c("至少一个P≤alpha/2", "两P都在(alpha/2, alpha]", "两项H0真时FWER"),
  probability = c(prob_any_small, prob_both_middle, fwer_two_null),
  row.names = NULL
)
print(event_table)
cat("这是模型的理论误报概率，不是模拟结果或某一次观测的P值。\n")
cat("若只有1项H0为真，误拒绝该项必须满足P≤alpha，因此上界为alpha。\n")

# 这是单独给定的3项相互独立检验，不来自前面F-G、F-H、G-H共享数据的原例。
# 原始P为test1/test2/test3=.020/.023/.040，按FWER=.05作判断。
p_demo <- c(test1 = 0.020, test2 = 0.023, test3 = 0.040)
practice_table <- data.frame(
  test = unname(names(p_demo)),
  raw = unname(p_demo),
  holm = stats::p.adjust(p_demo, method = "holm"),
  hochberg = stats::p.adjust(p_demo, method = "hochberg"),
  row.names = seq_along(p_demo)
)
print(practice_table)
practice_alpha <- 0.05
practice_reject_holm <- practice_table$test[practice_table$holm <= practice_alpha]
practice_reject_hochberg <- practice_table$test[practice_table$hochberg <= practice_alpha]
cat("test1/test2/test3分别对应检验甲/乙/丙。\n")
cat("Hochberg目标是控制FWER；适用于相互独立或符合特定正依赖条件的检验。\n")
cat("本练习采用相互独立前提，原三组对分析继续使用Holm。\n")
cat("Holm拒绝项数：", length(practice_reject_holm), "；Hochberg拒绝项数：", length(practice_reject_hochberg), "。\n", sep = "")

# 对照：最大P首关通过时全部拒绝，但各项调整P仍可不同。
p_compare <- c(test1 = 0.001, test2 = 0.020, test3 = 0.040)
compare_table <- data.frame(
  test = unname(names(p_compare)),
  raw = unname(p_compare),
  hochberg = stats::p.adjust(p_compare, method = "hochberg"),
  row.names = seq_along(p_compare)
)
compare_reject_hochberg <- compare_table$test[compare_table$hochberg <= practice_alpha]
print(compare_table)
if (max(p_compare) <= practice_alpha) {
  cat("最大P=", max(p_compare), "≤", practice_alpha,
      "，Hochberg首关通过，拒绝全部三项。\n", sep = "")
} else {
  cat("最大P=", max(p_compare), ">", practice_alpha,
      "，首关未通过，继续按Hochberg规则判断。\n", sep = "")
}
cat("各项调整P由整组P从后往前取累计最小计算，不一定相同。\n")
cat("对照例Hochberg拒绝项数：", length(compare_reject_hochberg), "。\n", sep = "")

# 8. G-H效应量与区间报告 ----
# 只报告一个示例组对，区间为单项区间，P值仍按三对家族调整。
gh <- category_pairs[category_pairs$pair == "G-H", ]
cat(sprintf("G-H电极脱落比例差=%.6f（前组减后组），单项%.0f%% CI [%.6f, %.6f]，Holm P=%.6f。\n",
            gh$diff, conf_level * 100, gh$diff_low, gh$diff_high, gh$p_holm))
cat("区间没有作三项同时覆盖调整；Holm只用于同一问题的三对P值。\n")

# 9. 可选保存 ----
# 默认不写文件。需保存时先运行 Sys.setenv(ARTICLE37_SAVE_RESULTS = "1")，再运行本节。
save_results <- identical(Sys.getenv("ARTICLE37_SAVE_RESULTS"), "1")
if (save_results) {
  output_dir <- file.path("results", "R")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list(main_overall = main_overall, full_pairs = full_pairs,
                  category_pairs = category_pairs, sparse_pairs = sparse_pairs)
  for (name in names(results)) {
    write.csv(results[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  }
  cat("已保存四张表到", output_dir, "\n")
}
