# 两个独立样本率之差的精确检验与置信区间
# 精确方法也适用于大样本，无需先满足正态近似门槛。
# 本脚本用小样本演示函数接口；扩大计算规模时关注运行时间与数值精度。
# 首次使用时，在R控制台安装下面两个包。
# install.packages(c("exact2x2", "ExactCIdiff"))
# 打开本篇R项目后运行，或将工作目录设为本篇“配套项目”。

required_packages <- c("exact2x2", "ExactCIdiff")
packages_available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(packages_available)) {
  stop("请先安装exact2x2和ExactCIdiff，再运行本脚本。")
}

here <- getwd()
if (basename(here) != "代码") here <- file.path(here, "代码")
project_root <- normalizePath(file.path(here, ".."))
output_directory <- file.path(project_root, "运行结果", "R")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# 先算观察到的样本率差。两组对象独立，组内采用独立同概率的二项模型。
positive_first <- 2
sample_size_first <- 10
positive_second <- 0
sample_size_second <- 10
rate_first <- positive_first / sample_size_first
rate_second <- positive_second / sample_size_second
observed_difference <- rate_first - rate_second
cat("第一组样本率、第二组样本率及两者之差\n")
print(c(first = rate_first, second = rate_second, difference = observed_difference))

# 检验两总体率相等时，按共同率计算零假设下每组的近似检查值。
pooled_rate <- (positive_first + positive_second) / (sample_size_first + sample_size_second)
null_check <- c(sample_size_first, sample_size_second) * pooled_rate * (1 - pooled_rate)
cat("零假设下两组计数方差的估计，均为0.9，低于Rosner经验门槛5\n")
print(null_check)

# exact2x2接口计算第2组减第1组，交换输入后对应本例第一组减第二组。
# score用于给可能的结果排序，P值由无条件精确方法的数值搜索得到。
# central指定双侧方式，midp=FALSE保留标准精确方法设置。
exact_result <- exact2x2::uncondExact2x2(
  x1 = positive_second, n1 = sample_size_second,
  x2 = positive_first, n2 = sample_size_first,
  parmtype = "difference", nullparm = 0,
  alternative = "two.sided", method = "score",
  tsmethod = "central", midp = FALSE,
  conf.int = TRUE, conf.level = 0.95,
  control = exact2x2::ucControl(nPgrid = 1000, nCIgrid = 1000)
)
print(exact_result)

# ExactCIdiff接口直接计算x/n1-y/n2，按本例组别顺序输入。
# 该函数采用另一种区间构造，端点可以与上一个函数不同。
ci_result <- ExactCIdiff::BinomCI(
  n1 = sample_size_first, n2 = sample_size_second,
  x = positive_first, y = positive_second,
  conf.level = 0.95, CItype = "Two.sided",
  grid.one = 60, grid.two = 40
)
print(ci_result)

# 两种区间都包含0，显示这批小样本仍有较大的率差估计不确定性。
# BinomCI输出区间，其P值栏保留NA。
results <- data.frame(
  method = c("exact2x2_score_central", "ExactCIdiff_inductive"),
  estimate = c(as.numeric(exact_result$estimate), ci_result$estimate),
  lower = c(exact_result$conf.int[1], ci_result$ExactCI[1]),
  upper = c(exact_result$conf.int[2], ci_result$ExactCI[2]),
  p_value = c(exact_result$p.value, NA_real_)
)
print(results)
write.csv(results, file.path(output_directory, "率差精确方法比较.csv"), row.names = FALSE)
cat("包版本\n")
print(c(exact2x2 = as.character(packageVersion("exact2x2")), ExactCIdiff = as.character(packageVersion("ExactCIdiff"))))
