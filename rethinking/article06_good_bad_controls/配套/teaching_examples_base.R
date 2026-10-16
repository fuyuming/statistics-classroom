#!/usr/bin/env Rscript
# 四个完全假设数据的确定性教学算例（R 基础版）。
# 只使用枚举、加权平均和 OLS；不调用 quap，也不做重复模拟。

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "teaching_examples_base.R"
root <- normalizePath(dirname(script_path), mustWork = TRUE)
input_path <- file.path(root, "input_data.csv")
output_path <- file.path(root, "results", "r_base_results.csv")

dat <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE,
                na.strings = c("", "NA"))
weighted_mean <- function(d, value = "mean_value", weight = "count") {
  sum(d[[value]] * d[[weight]], na.rm = TRUE) / sum(d[[weight]], na.rm = TRUE)
}
add_result <- function(example, metric, value) {
  data.frame(example = example, metric = metric, value = as.numeric(value),
             stringsAsFactors = FALSE)
}

out <- list()

# 算例 1：原始分层比例与固定为各 50% 的标准化均值。
ex1 <- dat[dat$example == 1, ]
raw_trained <- weighted_mean(ex1[ex1$trained == 1, ])
raw_untrained <- weighted_mean(ex1[ex1$trained == 0, ])
std_trained <- mean(ex1$mean_value[ex1$trained == 1])
std_untrained <- mean(ex1$mean_value[ex1$trained == 0])
out <- c(out, list(
  add_result(1, "raw_trained_mean", raw_trained),
  add_result(1, "raw_untrained_mean", raw_untrained),
  add_result(1, "raw_difference", raw_trained - raw_untrained),
  add_result(1, "standardized_trained_mean", std_trained),
  add_result(1, "standardized_untrained_mean", std_untrained),
  add_result(1, "standardized_difference", std_trained - std_untrained)
))

# 算例 2：条件概率 P(乙运动 | 甲运动) 与 P(乙运动 | 甲不运动)。
ex2 <- dat[dat$example == 2, ]
for (phase in c("before", "after")) {
  d <- ex2[ex2$phase == phase, ]
  cell <- setNames(d$count, d$friend_cell)
  p_b_given_a <- cell[["both"]] / (cell[["both"]] + cell[["a_only"]])
  p_b_given_not_a <- cell[["b_only"]] / (cell[["b_only"]] + cell[["neither"]])
  out <- c(out, list(
    add_result(2, paste0(phase, "_P_B_exercises_given_A_exercises"), p_b_given_a),
    add_result(2, paste0(phase, "_P_B_exercises_given_A_not_exercises"), p_b_given_not_a)
  ))
}

# 算例 3：Y~X 与 Y~X+Z 的 OLS 系数；不输出显著性或 SE。
ex3 <- dat[dat$example == 3, ]
fit_unadjusted <- lm(y ~ x, data = ex3)
fit_adjusted <- lm(y ~ x + z, data = ex3)
out <- c(out, list(
  add_result(3, "ols_y_on_x_x_coefficient", coef(fit_unadjusted)[["x"]]),
  add_result(3, "ols_y_on_x_plus_z_x_coefficient", coef(fit_adjusted)[["x"]])
))

# 算例 4：按培训率加权得到两年的总均值，并固定培训状态比较经验差。
ex4 <- dat[dat$example == 4, ]
year1 <- ex4[ex4$experience_years == 1, ]
year2 <- ex4[ex4$experience_years == 2, ]
rate1 <- unique(year1$training_rate)
rate2 <- unique(year2$training_rate)
mean1 <- sum(year1$mean_value * ifelse(year1$trained == 1, rate1, 1 - rate1))
mean2 <- sum(year2$mean_value * ifelse(year2$trained == 1, rate2, 1 - rate2))
untrained <- setNames(ex4$mean_value[ex4$trained == 0], ex4$experience_years[ex4$trained == 0])
trained <- setNames(ex4$mean_value[ex4$trained == 1], ex4$experience_years[ex4$trained == 1])
out <- c(out, list(
  add_result(4, "year1_overall_mean", mean1),
  add_result(4, "year2_overall_mean", mean2),
  add_result(4, "overall_difference", mean2 - mean1),
  add_result(4, "experience_difference_fixed_untrained", untrained[["2"]] - untrained[["1"]]),
  add_result(4, "experience_difference_fixed_trained", trained[["2"]] - trained[["1"]]),
  add_result(4, "training_difference_year1", trained[["1"]] - untrained[["1"]]),
  add_result(4, "training_difference_year2", trained[["2"]] - untrained[["2"]])
))

result <- do.call(rbind, out)
expected <- c(raw_difference = 28, standardized_difference = 10,
  before_P_B_exercises_given_A_exercises = 0.5,
  before_P_B_exercises_given_A_not_exercises = 0.5,
  after_P_B_exercises_given_A_exercises = 0.8,
  after_P_B_exercises_given_A_not_exercises = 0.2,
  ols_y_on_x_x_coefficient = 1 / 130,
  ols_y_on_x_plus_z_x_coefficient = 0.2,
  overall_difference = 6, experience_difference_fixed_untrained = 2,
  experience_difference_fixed_trained = 2)
for (name in names(expected)) {
  observed <- result$value[result$metric == name]
  stopifnot(length(observed) == 1, isTRUE(all.equal(as.numeric(observed), expected[[name]], tolerance = 1e-10)))
}

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(result, output_path, row.names = FALSE, fileEncoding = "UTF-8")
cat("R 基础版：四个算例的预期数字核对通过。\n")
cat("结果已写入：", output_path, "\n", sep = "")
