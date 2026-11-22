# 统计课堂39：有序多分类的加权κ。UTF-8；数据为假设的100名独立患者。
# 指标为过去24小时疼痛对日常活动的影响，A/B各给一次等级，盲评同次就诊资料。
# 无=没有疼痛；轻=有疼痛但不影响活动；中=影响活动但仍能完成；重=部分活动无法完成。
# 打开统计课堂39.Rproj，在项目根目录按节运行；Source可运行全文件。
# 依赖DescTools、vcd。默认在Console显示；可保存到运行结果/R/。
# ARTICLE39_SAVE_RESULTS=1用于批量保存；文章配图可单独重建。

# ---- 1. 读入完整共享配对表 ----
root <- getwd()
if (basename(root) == "代码") root <- dirname(root)
if (!file.exists(file.path(root, "数据", "paired_tables.csv"))) {
  stop("请先打开本篇Rproj，使工作目录为项目根。")
}
source(file.path(root, "代码", "kappa_parts.R"), encoding="UTF-8")
raw <- read.csv(file.path(root, "数据", "paired_tables.csv"))
print(head(raw)) # A/B编码0无、1轻、2中、3重；count为该联合格人数。
case_names <- c("main", "control")
expected_n <- 100 # 本题每张表100人；计算分母由实际计数求和。
stopifnot(setequal(unique(raw$case), case_names))
tables <- list()
for (case in case_names) {
  data <- raw[raw$case == case, ] # $取case列；[]按条件筛选当前主例或对照。
  stopifnot(nrow(data) == 16, all(data$A %in% 0:3), all(data$B %in% 0:3))
  stopifnot(all(data$count >= 0), all(data$count == floor(data$count)))
  stopifnot(!anyDuplicated(data[c("A", "B")]))
  tab <- matrix(0, 4, 4, dimnames=list(A=c("无","轻","中","重"), B=c("无","轻","中","重")))
  # CSV编码从0开始，R位置从1开始；按联合格而非CSV行号放入表。
  tab[cbind(data$A + 1, data$B + 1)] <- data$count # +1将编码变成1—4位置，cbind将行、列位置配成一对。
  if (sum(tab) != expected_n) stop(paste(case, "总人数应为", expected_n, "，实际为", sum(tab)))
  tables[[case]] <- tab
  print(case)
  print(tab)
}
# 日常练习先只看Console；在此处把这个布尔值改为TRUE即可保存结果。
save_results <- FALSE

# ---- 2. 构造一致权重 ----
# outer逐一作两边等级的差；abs取绝对值，得到16格的等级步数。
distance <- abs(outer(0:3, 0:3, "-"))
weights <- list(unweighted=diag(4), linear=1-distance/3, quadratic=1-(distance/3)^2)
print(weights$linear)
print(weights$quadratic)

# ---- 3.1 主例线性κ：把正文手算逐步写成代码 ----
main_manual <- tables$main
manual_n <- sum(main_manual)
manual_exact <- sum(diag(main_manual)) / manual_n # diag取主对角线，即A/B完全同级格。
manual_p <- main_manual / manual_n
manual_row <- rowSums(manual_p) # rowSums横向求和，得到A行边际。
manual_col <- colSums(manual_p) # colSums纵向求和，得到B列边际。
manual_chance <- outer(manual_row, manual_col) # 外积构造A、B独立时的16格参照比例。
manual_weight <- weights$linear
manual_po <- sum(manual_weight * manual_p) # 逐格相乘一致权重与真实比例，再求和。
manual_pe <- sum(manual_weight * manual_chance) # 同样逐格加权独立参照，再求和。
manual_kappa <- (manual_po - manual_pe) / (1 - manual_pe)
distance_steps <- data.frame()
for (d in 0:3) {
  mask <- distance == d # mask筛选同距离格；同距离格共享同一权重。
  patients <- sum(main_manual[mask]) # patients是该距离的患者人数。
  chance_cells <- sum(mask) # chance_cells是该距离的格子数。
  chance_probability <- sum(manual_chance[mask]) # 边际外积参照比例之和。
  linear_weight <- manual_weight[mask][1] # 同距离权重相同，取首格即可。
  Po_term <- linear_weight * patients / manual_n
  Pe_term <- linear_weight * chance_probability
  distance_steps <- rbind(distance_steps, data.frame(distance=d, patients=patients,
    chance_cells=chance_cells, chance_probability=chance_probability,
    linear_weight=linear_weight, Po_term=Po_term, Pe_term=Pe_term))
}
stopifnot(abs(sum(distance_steps$Po_term)-manual_po) < 1e-12, abs(sum(distance_steps$Pe_term)-manual_pe) < 1e-12)
linear_steps <- data.frame(n=manual_n, exact_agreement=manual_exact,
                           Po=manual_po, Pe=manual_pe, kappa=manual_kappa)
cat(sprintf("3.1主例：n=%d；完全同级比例=%.2f%%（回答有多少人两次评级相同）。\n", manual_n, 100 * manual_exact))
print("每格联合比例p：")
print(manual_p)
print("行边际r与列边际c：")
print(manual_row)
print(manual_col)
print("独立机会参照r×c：")
print(manual_chance)
cat(sprintf("线性Po=%.7f（真实配对的加权一致量）；Pe=%.7f（边际独立时的参照量）；κ=%.7f（相对机会基准的加权一致性）。\n", manual_po, manual_pe, manual_kappa))
print("linear_steps（3.1逐步结果）：")
print(linear_steps)
print(distance_steps)

# ---- 3.2 全表比较、一般SE/Wald及原生核对 ----
# DescTools默认未加权；这里每次显式传规则，conf.level=.95要求区间。
control <- tables$control
distance_comparison <- data.frame()
for (d in 0:3) {
  mask <- distance == d
  distance_comparison <- rbind(distance_comparison, data.frame(distance=d,
    main_count=sum(main_manual[mask]), control_count=sum(control[mask])))
}
control_checks <- data.frame(
  same_diagonal=all(diag(main_manual) == diag(control)),
  same_row_totals=all(rowSums(main_manual) == rowSums(control)),
  same_column_totals=all(colSums(main_manual) == colSums(control)))
print(distance_comparison)
print(control_checks)
if (!all(control_checks[1,])) {
  stop("control应保持主例的对角和行列合计，请检查分歧移动。")
}

package_weights <- c(unweighted="Unweighted", linear="Equal-Spacing", quadratic="Fleiss-Cohen")
results <- data.frame()
for (case in names(tables)) {
  tab <- tables[[case]]
  exact_agreement <- sum(diag(tab)) / sum(tab) # diag取完全同级格；除以总人数得到比例。
  cat(sprintf("%s完全同级比例：%.2f%%（对角人数/总人数）。\n", case, 100*exact_agreement))
  for (method in names(weights)) {
    # tables[[case]]取当前表；weights[[method]]取当前规则，循环覆盖2表×3权重。
    parts <- kappa_parts(tab, weights[[method]]) # [[按名称取当前权重矩阵。
    observed_loss <- unname(1 - parts["Po"])
    expected_loss <- unname(1 - parts["Pe"])
    loss_ratio <- if (expected_loss > 1e-14) observed_loss / expected_loss else NA_real_
    native <- c(kappa=NA_real_, lwr.ci=NA_real_, upr.ci=NA_real_)
    if (is.finite(parts["kappa"])) {
      native <- DescTools::CohenKappa(tab, weights=package_weights[[method]], conf.level=.95)
      stopifnot(max(abs(parts[c("kappa", "CI_low", "CI_high")] - native)) < 1e-10)
    }
    # parts位置/字段依次是Po、Pe、kappa、SE、Wald下限、上限。
    row <- data.frame(case=case, method=method, t(parts), exact_agreement=exact_agreement,
                      Lo=observed_loss, Le=expected_loss, loss_ratio=loss_ratio,
                      DescTools_kappa=unname(native["kappa"]))
    results <- rbind(results, row)
    print(row, row.names=FALSE)
    cat(sprintf("实际平均损失Lo=%.6f，机会平均损失Le=%.6f，Lo/Le=%.6f。\n", observed_loss, expected_loss, loss_ratio))
    if (method != "unweighted" && is.finite(parts["kappa"])) {
      fit <- vcd::Kappa(tab, weights=package_weights[[method]])
      stopifnot(abs(fit$Weighted["value"] - parts["kappa"]) < 1e-10)
      stopifnot(abs(fit$Weighted["ASE"] - parts["SE"]) < 1e-10)
      print(confint(fit, level=.95)) # Weighted行给本节加权区间。
    }
  }
}
cat("完全同级比例另由对角人数/总人数计算；κ描述指定规则下相对机会基准的一致性。\n")
# &要求两个条件同时成立，选出主例采用线性权重的那一行。
main_linear <- results[results$case == "main" & results$method == "linear", ][1, ]
main_quadratic <- results[results$case == "main" & results$method == "quadratic", ][1, ]
stopifnot(max(abs(c(linear_steps$Po, linear_steps$Pe, linear_steps$kappa) -
                 c(main_linear$Po, main_linear$Pe, main_linear$kappa))) < 1e-12)
function_comparison <- data.frame(case=results$case, method=results$method,
                                   kappa_formula=results$kappa,
                                   kappa_native=results$DescTools_kappa)
print("function_comparison（自编公式与DescTools原生点估计）：")
print(function_comparison)
# report_text为由当前结果格式化生成的报告句，不调用新的report软件。
report_text <- sprintf("对%d名患者按无、轻、中、重顺序分析；预设线性一致权重为1−|i−j|/3。完全同级比例为%.1f%%，线性Po=%.4f、Pe=%.4f、κ=%.3f，95%%Wald置信区间为[%.3f, %.3f]；二次敏感性分析κ=%.3f，95%%Wald置信区间为[%.3f, %.3f]。",
                       manual_n, 100*main_linear$exact_agreement, main_linear$Po, main_linear$Pe, main_linear$kappa,
                       main_linear$CI_low, main_linear$CI_high, main_quadratic$kappa, main_quadratic$CI_low, main_quadratic$CI_high)
cat("report_text：", report_text, "\n")

# ---- 4. 主例成对Bootstrap，每轮重新算边际 ----
tab <- tables$main
# R按列展开；每个编号对应一个联合格，次数是该格患者数。
pair_cells <- rep(seq_len(16), as.vector(tab))
repetitions <- 10000
set.seed(20261122)
boot_values <- matrix(NA_real_, repetitions, 2)
for (b in seq_len(repetitions)) {
  # 对100个患者编号有放回抽100次，同人的A/B由格编号一起保留。
  sampled_ids <- sample(pair_cells, size=length(pair_cells), replace=TRUE) # replace=TRUE表示患者有放回重抽。
  sampled <- matrix(tabulate(sampled_ids, nbins=16), 4, 4) # tabulate计数后按同一列展开顺序重构4×4表。
  boot_values[b, 1] <- kappa_parts(sampled, weights$linear)["kappa"]
  boot_values[b, 2] <- kappa_parts(sampled, weights$quadratic)["kappa"]
}
bootstrap_results <- data.frame()
for (j in 1:2) {
  valid <- is.finite(boot_values[, j])
  if (!any(valid)) stop("全部重抽κ无定义，请检查数据。")
  # .025、.975是双侧95%区间的端点概率；type=7在位置1+(B-1)q处线性插值。
  # 这一分位数规则与Python/MATLAB配套一致。
  interval <- quantile(boot_values[valid, j], c(.025, .975), type=7, names=FALSE)
  row <- data.frame(method=c("linear", "quadratic")[j], B=repetitions,
                    valid=sum(valid), undefined=sum(!valid), CI_low=interval[1], CI_high=interval[2])
  bootstrap_results <- rbind(bootstrap_results, row)
}
print(bootstrap_results, row.names=FALSE)

# ---- 5. 可选保存；不改变默认首次练习的Console查看方式 ----
# 日常只改第1节的save_results；环境变量1用于批量复现时覆盖为TRUE。
if (Sys.getenv("ARTICLE39_SAVE_RESULTS") == "1") save_results <- TRUE
if (save_results) {
  out <- file.path(root, "运行结果", "R")
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  write.csv(results, file.path(out, "summary.csv"), row.names=FALSE)
  write.csv(bootstrap_results, file.path(out, "bootstrap.csv"), row.names=FALSE)
}
# 练习：先预测control表的κ，再运行核对；权重按研究损失选，不能按κ高低选。
