# 统计课堂41：双向有序表的条件模拟。UTF-8。
# 目标：用固定行列边际的完整枚举和Monte Carlo理解双向有序表。
# 共享数据为医学统计学第5版p116–117例7-9（主例）及另设8例假设教学资料；输入为data/ordinal_tables.csv。
# 实际依赖：DescTools、coin；结果默认在Console/Plots显示，ARTICLE41_SAVE_RESULTS=1时写入results/和figures/。
# RStudio打开本项目Rproj，打开此脚本，从01节依次运行；也可Source。
# 输入是共享计数CSV；输出先显示在Console和Plots，保存默认关闭。

# 01 数据、参数与函数 ----
stopifnot(file.exists("data/ordinal_tables.csv"))
source("scripts/helpers/ordinal_helpers.R", encoding = "UTF-8")
raw <- read.csv("data/ordinal_tables.csv", stringsAsFactors = FALSE)
cat("输入字段示例：\n")
print(utils::head(raw, 3))
if (identical(Sys.getenv("ARTICLE41_SHOW_DETAILS"), "1")) print(raw)
B <- 19999L
SEED <- 20261122L
row_scores <- 1:4
col_scores <- 1:4
SAVE_RESULTS <- identical(Sys.getenv("ARTICLE41_SAVE_RESULTS"), "1")
SHOW_DETAILS <- identical(Sys.getenv("ARTICLE41_SHOW_DETAILS"), "1")
cat("参数：B=", B, ", SEED=", SEED, ", SAVE_RESULTS=", as.integer(SAVE_RESULTS), ", SHOW_DETAILS=", as.integer(SHOW_DETAILS), "\n", sep = "")
get_table <- function(name) {
  part <- raw[raw$example == name, ]
  part <- part[order(part$row_score, part$col_score), ]
  matrix(part$count, nrow = 4, byrow = TRUE)
}
main <- get_table("main")
mini <- get_table("mini")
if (SHOW_DETAILS) { print(main); print(mini) }

# 02 两张观察表：T、r、条件方差与渐近P ----
observed <- rbind(main = ordinal_stats(main, row_scores, col_scores),
                  mini = ordinal_stats(mini, row_scores, col_scores))
cat("观察摘要（chi2_p为χ²1渐近P；不是完整条件尾部P）\n")
print(observed[, c("N", "T", "r", "Q", "chi2_p")], digits = 11)
cat("T是中心化评分积和；Lxx/Lyy是两项评分的离差平方和，VarT是固定边际下T的方差。\n")
cat("r是标准化评分关联，Q=(N-1)r²用于χ²1近似；chi2_p不是完整条件尾部P。\n")

formula_steps <- function(tab, ans, name, show_details = FALSE) {
  n <- sum(tab)
  row_total <- sum(tab * row_scores[row(tab)])
  col_total <- sum(tab * col_scores[col(tab)])
  row_mean <- row_total / n
  col_mean <- col_total / n
  row1_lxx <- sum(tab[1, ] * (row_scores[1] - row_mean)^2)
  cells <- data.frame(
    example = name,
    quantity = c("row_score_total", "col_score_total", "row_mean", "col_mean",
                 "sqrt(Lxx*Lyy)", "VarT", "SDT", "Z", "row_1_Lxx_contribution",
                 "(1,1)", "(1,4)"),
    value = c(row_total, col_total, row_mean, col_mean, sqrt(ans["Lxx"] * ans["Lyy"]),
              ans["VarT"], sqrt(ans["VarT"]), ans["T"] / sqrt(ans["VarT"]), row1_lxx, NA, NA),
    count = c(rep(NA, 9), tab[1, 1], tab[1, 4]),
    centered_product = c(rep(NA, 9), (row_scores[1] - row_mean) * (col_scores[1] - col_mean),
                         (row_scores[1] - row_mean) * (col_scores[4] - col_mean)),
    T_contribution = c(rep(NA, 9), tab[1, 1] * (row_scores[1] - row_mean) * (col_scores[1] - col_mean),
                       tab[1, 4] * (row_scores[1] - row_mean) * (col_scores[4] - col_mean))
  )
  if (!show_details) return(cells)
  cat("\n[", name, "] 演算中间量\n", sep = "")
  cat("两评分总和=", row_total, ", ", col_total, "；均值=", format(row_mean, digits = 12),
      ", ", format(col_mean, digits = 12), "\n", sep = "")
  cat("第一行（评分", row_scores[1], "）对Lxx贡献=", format(row1_lxx, digits = 12), "\n", sep = "")
  cat("格(1,1)：一人的中心化积=", format(cells$centered_product[10], digits = 12),
      "；乘人数后的T贡献=", format(cells$T_contribution[10], digits = 12), "\n", sep = "")
  cat("格(1,4)：一人的中心化积=", format(cells$centered_product[11], digits = 12),
      "；乘人数后的T贡献=", format(cells$T_contribution[11], digits = 12), "\n", sep = "")
  cat("sqrt(Lxx*Lyy)=", format(sqrt(ans["Lxx"] * ans["Lyy"]), digits = 12),
      "；VarT=", format(ans["VarT"], digits = 12), "；SDT=", format(sqrt(ans["VarT"]), digits = 12),
      "；Z=", format(ans["T"] / sqrt(ans["VarT"]), digits = 12), "\n", sep = "")
  cells
}
formula_rows <- rbind(formula_steps(main, observed["main", ], "main", SHOW_DETAILS),
                      formula_steps(mini, observed["mini", ], "mini", SHOW_DETAILS))

# 03 小表的完整枚举；每张表按条件概率加权 ----
tabs <- enumerate_tables(rowSums(mini), colSums(mini))
cols <- colSums(mini)
Tvals <- vapply(tabs, function(z) ordinal_stats(z, row_scores, col_scores)["T"], numeric(1))
probabilities <- vapply(tabs, table_probability, numeric(1))
enum_df <- data.frame(T = Tvals, probability = probabilities)
enum_df$extreme <- abs(Tvals) >= abs(observed["mini", "T"]) - 1e-10
exact_p <- sum(probabilities[enum_df$extreme])
cat("完整枚举摘要：可能表数：", length(tabs), "，概率和：", sum(probabilities),
    "，极端表数：", sum(enum_df$extreme), "，加权精确P：", exact_p, "\n")
# 预期数只用来独立核对，结果由输入和枚举计算。
stopifnot(length(tabs) == 282, abs(sum(probabilities) - 1) < 1e-12)
stopifnot(abs(exact_p - 4 / 105) < 1e-12)
if (SHOW_DETAILS) cat("精确P回答预定评分的双侧条件尾部，不能用极端表张数/总张数代替。\n")
tail_source <- enum_df[enum_df$extreme, ]
tail_summary <- aggregate(list(probability_sum = tail_source$probability), list(T = tail_source$T), sum)
counts <- aggregate(list(table_count = tail_source$T), list(T = tail_source$T), length)
names(tail_summary)[2] <- "probability_sum"
names(counts)[2] <- "table_count"
tail_summary <- merge(counts, tail_summary, by = "T")
tail_summary <- tail_summary[order(tail_summary$T), ]
total_label_sequences <- factorial(sum(cols)) / prod(factorial(cols))
tail_summary$label_sequence_count <- tail_summary$probability_sum * total_label_sequences
if (SHOW_DETAILS) { cat("\n[mini] 极端T分组（由实际enum_df归并）\n"); print(tail_summary, row.names = FALSE, digits = 12) }
if (SHOW_DETAILS) cat("极端标签序列=", round(sum(tail_summary$label_sequence_count)),
    "/总标签序列=", total_label_sequences, "（核对常数2520）；理论尾部期望超越次数(B×精确p)=", format(B * exact_p, digits = 12),
    "，不等于某次实际k。\n", sep = "")

# 04 固定边际Monte Carlo及计算精度 ----
mc_mini <- conditional_mc(mini, B, SEED, row_scores, col_scores)
mc_main <- conditional_mc(main, B, SEED, row_scores, col_scores)
mc_summary <- rbind(mini = mc_mini$summary, main = mc_main$summary)
print(mc_summary, digits = 11)
cat("k_over_B估计固定观察表的条件尾概率；加一P把观察结果计入，用作本文Monte Carlo检验P值。\n")
cat("MCSE和CP95是尾部概率的计算精度，不是r的置信区间；零超越不等于真实P为0。\n")
mini_first <- mc_mini$first_table
cat("小例第一轮模拟表（同一批B轮中的第一张，Sum为行列合计）：\n")
print(addmargins(mini_first))
first_T <- unname(ordinal_stats(mini_first, row_scores, col_scores)["T"])
stopifnot(isTRUE(all.equal(first_T, unname(mc_mini$T_sim[1]))))
cat("第一张表T=", format(first_T, digits = 12), "，与mini_sim_T首行相同。\n", sep = "")
main_report <- paste0(as.integer(observed["main", "N"]), "例资料，按当前预设评分，r=",
  format(observed["main", "r"], digits = 6), "，Q=", format(observed["main", "Q"], digits = 6),
  "；固定边际双侧条件模拟B=", B, "，极端", mc_main$summary["k"], "次，加一P=",
  format(mc_main$summary["plus1_over_plus1"], digits = 8), "；")
if (mc_main$summary["k"] == 0) main_report <- paste0(main_report,
  "条件尾概率的单侧95%二项上限=", format(mc_main$summary["zero_exceedance_one_sided_95_upper"], digits = 8), "。")
cat("主例报告：", main_report, "\n", sep = "")

# 05 完整检验函数：输入观察表，返回检验结果 ----
# 若包缺失会提示安装，不把缺包伪装成已执行。
if (!requireNamespace("DescTools", quietly = TRUE)) stop("请安装DescTools包。")
if (!requireNamespace("coin", quietly = TRUE)) stop("请安装coin包。")
print(DescTools::MHChisqTest(main, srow = row_scores, scol = col_scores))
print(DescTools::MHChisqTest(mini, srow = row_scores, scol = col_scores))
# 本篇默认评分1:4。表对象的维度名称标明两边，避免依赖字母序。
dimnames(mini) <- list(X = as.character(1:4), Y = as.character(1:4))
set.seed(SEED)
fit <- coin::lbl_test(as.table(mini), alternative = "two.sided",
  scores = list(X = row_scores, Y = col_scores),
  distribution = coin::approximate(nresample = B))
cat("coin标准化Z：", as.numeric(coin::statistic(fit, type = "standardized")),
    "；P：", as.numeric(coin::pvalue(fit)), "\n")
if (SHOW_DETAILS) {
  exact_message <- tryCatch({
    coin::lbl_test(as.table(mini), distribution = "exact")
    "exact调用成功"
  }, error = function(e) conditionMessage(e))
  cat("当前4×4的exact接口：", exact_message, "\n")
}

# 06 先查看分布图；需要时保存CSV和图 ----
if (SAVE_RESULTS) {
  dir.create("figures", showWarnings = FALSE)
  png("figures/ordinal_conditional_mc_R.png", width = 1200, height = 1000)
}
# Rscript没有交互图窗时用空设备计算图；设保存开关才写PNG。
if (!SAVE_RESULTS && !interactive()) grDevices::pdf(file = NULL)
plot_distribution(enum_df, mc_mini$T_sim, observed["mini", "T"])
if (!SAVE_RESULTS && !interactive()) dev.off()
if (SAVE_RESULTS) {
  dev.off()
  dir.create("results", showWarnings = FALSE)
  write.csv(observed, "results/R_observed_stats.csv")
  write.csv(enum_df, "results/R_mini_exact_distribution.csv", row.names = FALSE)
  write.csv(formula_rows, "results/R_formula_steps.csv", row.names = FALSE)
  write.csv(tail_summary, "results/R_mini_tail_summary.csv", row.names = FALSE)
  write.csv(mc_summary, "results/R_mc_results.csv")
  write.csv(data.frame(T = mc_mini$T_sim), "results/R_mini_sim_T.csv", row.names = FALSE)
  write.table(mc_mini$first_table, "results/R_mini_first_table.csv", sep = ",", row.names = FALSE, col.names = FALSE)
}
# 练习：将B改为4999/79999，从01节重跑；比较MCSE，保持观察资料不变。
