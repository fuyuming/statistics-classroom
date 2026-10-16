#!/usr/bin/env Rscript
# 四个完全假设数据的确定性教学算例（R tidyverse 版）。
# 只使用 tidyverse 的数据整理语法和 stats::lm；不调用 quap，也不做重复模拟。

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  stop("本脚本需要 tidyverse；未安装时请按 README 做静态检查，勿在脚本中自动安装依赖。")
}
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "teaching_examples_tidyverse.R"
root <- normalizePath(dirname(script_path), mustWork = TRUE)
dat <- readr::read_csv(file.path(root, "input_data.csv"), show_col_types = FALSE)
output_path <- file.path(root, "results", "r_tidyverse_results.csv")

add_result <- function(example, metric, value) {
  tibble(example = example, metric = metric, value = as.numeric(value))
}
out <- list()

# 算例 1：先按真实人数加权，再把两档各设为 50%。
ex1 <- dat %>% filter(example == 1)
raw <- ex1 %>% group_by(trained) %>% summarise(mean_value = weighted.mean(mean_value, count), .groups = "drop")
std <- ex1 %>% group_by(trained) %>% summarise(mean_value = mean(mean_value), .groups = "drop")
raw_trained <- raw %>% filter(trained == 1) %>% pull(mean_value)
raw_untrained <- raw %>% filter(trained == 0) %>% pull(mean_value)
std_trained <- std %>% filter(trained == 1) %>% pull(mean_value)
std_untrained <- std %>% filter(trained == 0) %>% pull(mean_value)
out <- c(out, list(
  add_result(1, "raw_trained_mean", raw_trained),
  add_result(1, "raw_untrained_mean", raw_untrained),
  add_result(1, "raw_difference", raw_trained - raw_untrained),
  add_result(1, "standardized_trained_mean", std_trained),
  add_result(1, "standardized_untrained_mean", std_untrained),
  add_result(1, "standardized_difference", std_trained - std_untrained)
))

# 算例 2：条件概率 P(乙运动 | 甲运动) 与 P(乙运动 | 甲不运动)。
ex2 <- dat %>% filter(example == 2)
for (phase in c("before", "after")) {
  # .env$phase 明确引用循环变量，避免与数据列 phase 同名。
  d <- ex2 %>% filter(.data$phase == .env$phase) %>% select(friend_cell, count)
  both <- d %>% filter(friend_cell == "both") %>% pull(count)
  a_only <- d %>% filter(friend_cell == "a_only") %>% pull(count)
  b_only <- d %>% filter(friend_cell == "b_only") %>% pull(count)
  neither <- d %>% filter(friend_cell == "neither") %>% pull(count)
  out <- c(out, list(
    add_result(2, paste0(phase, "_P_B_exercises_given_A_exercises"), both / (both + a_only)),
    add_result(2, paste0(phase, "_P_B_exercises_given_A_not_exercises"), b_only / (b_only + neither))
  ))
}

# 算例 3：有限四行 OLS；真实效应设为零，Y 的差异由 U 造成。
ex3 <- dat %>% filter(example == 3)
fit_unadjusted <- stats::lm(y ~ x, data = ex3)
fit_adjusted <- stats::lm(y ~ x + z, data = ex3)
out <- c(out, list(
  add_result(3, "ols_y_on_x_x_coefficient", coef(fit_unadjusted)[["x"]]),
  add_result(3, "ols_y_on_x_plus_z_x_coefficient", coef(fit_adjusted)[["x"]])
))

# 算例 4：每一年的总体均值由该年的培训率混合得到。
ex4 <- dat %>% filter(example == 4)
overall <- ex4 %>%
  group_by(experience_years) %>%
  summarise(overall_mean = sum(mean_value * if_else(trained == 1, training_rate, 1 - training_rate)), .groups = "drop")
fixed <- ex4 %>% select(experience_years, trained, mean_value) %>% pivot_wider(names_from = experience_years, values_from = mean_value)
mean1 <- overall %>% filter(experience_years == 1) %>% pull(overall_mean)
mean2 <- overall %>% filter(experience_years == 2) %>% pull(overall_mean)
untrained <- fixed %>% filter(trained == 0)
trained <- fixed %>% filter(trained == 1)
out <- c(out, list(
  add_result(4, "year1_overall_mean", mean1),
  add_result(4, "year2_overall_mean", mean2),
  add_result(4, "overall_difference", mean2 - mean1),
  add_result(4, "experience_difference_fixed_untrained", untrained$`2` - untrained$`1`),
  add_result(4, "experience_difference_fixed_trained", trained$`2` - trained$`1`),
  add_result(4, "training_difference_year1", trained$`1` - untrained$`1`),
  add_result(4, "training_difference_year2", trained$`2` - untrained$`2`)
))

result <- bind_rows(out)
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
  observed <- result %>% filter(metric == name) %>% pull(value)
  stopifnot(length(observed) == 1, isTRUE(all.equal(observed, expected[[name]], tolerance = 1e-10)))
}

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(result, output_path)
cat("R tidyverse 版：四个算例的预期数字核对通过。\n")
cat("结果已写入：", output_path, "\n", sep = "")
