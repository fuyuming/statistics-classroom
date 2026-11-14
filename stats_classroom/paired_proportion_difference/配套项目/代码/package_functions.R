# 统计课堂27：多项分布函数核对（R）
# 问题：q+=(阴转阳)=0.05、q0(不变)=0.70、q−(阳转阴)=0.25时，
# 给定三类人数安排的概率是多少？这里比较n=2的两组和n=20主例。

# 问题/条件/数据身份：定位项目并查看共享数据 ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else ""
if (nzchar(script_path) && !file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
script_dir <- if (nzchar(script_path) && file.exists(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
candidates <- unique(c(getwd(), script_dir, dirname(script_dir)))
valid <- vapply(candidates, function(path) file.exists(file.path(path, "数据", "paired_20.csv")), logical(1))
if (!any(valid)) stop("找不到项目根目录；请从项目根或代码目录运行。")
project <- normalizePath(candidates[which(valid)[1]], mustWork = TRUE)
out_dir <- file.path(project, "运行结果", "R")
people <- read.csv(file.path(project, "数据", "paired_20.csv"), check.names = FALSE)
names(people)
head(people, 5)

# 条件/参数：人数列顺序为(+1,0,-1)，且三类人数总和等于n。
q_plus <- 0.05
q_zero <- 0.70
q_minus <- 0.25
cases <- data.frame(model = c("n2_020", "n2_101", "n20_1145"),
                    n = c(2L, 2L, 20L), n_plus = c(0L, 1L, 1L),
                    n_zero = c(2L, 0L, 14L), n_minus = c(0L, 1L, 5L))
stopifnot(all(cases$n_plus + cases$n_zero + cases$n_minus == cases$n))

# 手动计算：三组完整人数概率 ----
# 组合数给出类别安排数；一个具体安排的概率是三类概率幂的乘积。
cases$arrangements <- choose(cases$n, cases$n_plus) *
  choose(cases$n - cases$n_plus, cases$n_zero)
cases$single_arrangement_probability <- q_plus^cases$n_plus * q_zero^cases$n_zero * q_minus^cases$n_minus
cases$manual_probability <- cases$arrangements * cases$single_arrangement_probability
cases[, c("model", "n", "n_plus", "n_zero", "n_minus", "manual_probability")]

# 调用函数：stats::dmultinom ----
# counts与prob的顺序必须相同；size是总人数，而非重复次数。
cases$probability <- vapply(seq_len(nrow(cases)), function(i) {
  counts <- c(cases$n_plus[i], cases$n_zero[i], cases$n_minus[i])
  dmultinom(counts, size = cases$n[i], prob = c(q_plus, q_zero, q_minus))
}, numeric(1))
cases[, c("model", "n", "n_plus", "n_zero", "n_minus", "probability")]

# 结果核对：手算与函数逐组对应 ----
stopifnot(max(abs(cases$manual_probability - cases$probability)) < 1e-12)
result <- cases[, c("model", "n", "n_plus", "n_zero", "n_minus", "probability")]
result
# 可选保存：默认不写文件；改TRUE后更新一份概率CSV ----
save_results <- FALSE
if (save_results) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(out_dir, "package_functions_probabilities.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# 思考/可选练习：改一组人数，先检查三数之和，再比较手算和函数结果。
