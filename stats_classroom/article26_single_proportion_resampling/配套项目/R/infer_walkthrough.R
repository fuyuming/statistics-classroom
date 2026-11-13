# 统计课堂26：用 infer 逐步理解单样本率的 Bootstrap 与原假设模拟。
# 本节按“指定变量→生成样本→计算每轮统计量→得到分布与推断”的顺序编写，适合初学者。

# 1. 定位项目根目录、加载依赖与准备数据 ----
project_root <- getwd()
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  project_root <- normalizePath(file.path(dirname(script_path), ".."))
}
if (!file.exists(file.path(project_root, "data", "cases.csv"))) {
  stop("请从本项目运行；未找到 data/cases.csv。")
}
if (!requireNamespace("infer", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("本教程需要 infer 与 ggplot2。")
}
out_dir <- file.path(project_root, "results")
data_dir <- file.path(project_root, "data")
article_dir <- dirname(project_root)
dir.create(file.path(article_dir, "文章配图"), recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)
options(digits = 8)
font_family <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  available_families <- systemfonts::system_fonts()$family
  font_family <- if ("Heiti SC" %in% available_families) {
    "Heiti SC"
  } else if ("Noto Sans CJK SC" %in% available_families) {
    "Noto Sans CJK SC"
  } else if ("PingFang SC" %in% available_families) {
    "PingFang SC"
  } else {
    "sans"
  }
}
cat("Plot font family:", font_family, "\n")

sink(file.path(out_dir, "infer_v4_console.txt"), split = TRUE)
cat("R:", R.version.string, "\n")
cat("infer:", as.character(packageVersion("infer")), "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("Seed used for each independent block: 20261110\n")

# yes/no 用因子保存，并明确 no/yes 的水平顺序；success='yes' 即事件定义。
dat <- data.frame(
  outcome = factor(c(rep("yes", 117), rep("no", 63)), levels = c("no", "yes"))
)
cat("\nData counts:\n")
print(table(dat$outcome))

# 2. 先指定统计问题：观察到的达标率 ----
obs_stat <- dat |>
  infer::specify(response = outcome, success = "yes") |>
  infer::calculate(stat = "prop")
spec <- dat |>
  infer::specify(response = outcome, success = "yes")
cat("\nSpecification object (spec):\n")
print(spec)
cat("Observed statistic (obs_stat):\n")
print(as.data.frame(obs_stat), row.names = FALSE)

# 3. 用三轮 Bootstrap 看清 generate 与 calculate ----
set.seed(20261110)
demo_samples <- spec |>
  infer::generate(reps = 3, type = "bootstrap")
cat("\nThree-round bootstrap sample rows:\n")
cat("nrow(demo_samples):", nrow(demo_samples), "\n")
print(table(demo_samples$replicate))
demo_stats <- demo_samples |>
  infer::calculate(stat = "prop")
print(as.data.frame(demo_stats), row.names = FALSE)
write.csv(as.data.frame(demo_stats), file.path(data_dir, "infer_v4_demo_stats.csv"), row.names = FALSE)

# 4. 生成10000轮 Bootstrap 分布并计算百分位区间 ----
set.seed(20261110)
boot_dist <- dat |>
  infer::specify(response = outcome, success = "yes") |>
  infer::generate(reps = 10000, type = "bootstrap") |>
  infer::calculate(stat = "prop")
ci <- infer::get_confidence_interval(boot_dist, level = .95, type = "percentile")
cat("\nBootstrap confidence interval (percentile, 95%):\n")
print(as.data.frame(ci), row.names = FALSE)
write.csv(as.data.frame(boot_dist), file.path(data_dir, "infer_v4_bootstrap_reference.csv"), row.names = FALSE)
write.csv(as.data.frame(ci), file.path(data_dir, "infer_v4_bootstrap_ci.csv"), row.names = FALSE)

# native infer 图层先生成；额外的中文轴标题、观察率线和范围属于本脚本修饰。
native_boot_plot <- infer::visualize(boot_dist, bins = 25) + infer::shade_confidence_interval(ci)
boot_plot <- native_boot_plot +
  ggplot2::labs(title = "infer：从当前样本重抽 Bootstrap", x = "每轮达标率", y = "轮数") +
  ggplot2::geom_vline(xintercept = .65, linetype = "dotted", linewidth = .8, color = "#2345AB") +
  ggplot2::annotate("text", x = .65, y = Inf, label = "观察率 0.65", vjust = 1.5, color = "#2345AB", family = font_family) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_cartesian(xlim = c(.25, .80)) +
  ggplot2::theme_minimal(base_size = 16, base_family = font_family)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggplot2::ggsave(file.path(article_dir, "文章配图", "01_infer区间_v4.png"), boot_plot, width = 7, height = 4.8, dpi = 170, device = ragg::agg_png)
} else {
  ggplot2::ggsave(file.path(article_dir, "文章配图", "01_infer区间_v4.png"), boot_plot, width = 7, height = 4.8, dpi = 170)
}

# 5. H0:p=.45 下生成原假设分布，计算右尾P值 ----
set.seed(20261110)
null_dist <- dat |>
  infer::specify(response = outcome, success = "yes") |>
  infer::hypothesize(null = "point", p = .45) |>
  infer::generate(reps = 10000, type = "draw") |>
  infer::calculate(stat = "prop")
p_warning <- NULL
p_value <- withCallingHandlers(
  infer::get_p_value(null_dist, obs_stat = obs_stat, direction = "greater"),
  warning = function(w) {
    p_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
null_hits <- sum(null_dist$stat >= 117 / 180)
null_hit_rate <- null_hits / nrow(null_dist)
cat("\nNull distribution p-value (greater):\n")
print(as.data.frame(p_value), row.names = FALSE)
cat("Actual null hits (stat >= 117/180):", null_hits, "of", nrow(null_dist), "\n")
cat("Actual null hit proportion:", null_hit_rate, "\n")
if (!is.null(p_warning)) cat("Captured get_p_value warning:", p_warning, "\n") else cat("Captured get_p_value warning: none\n")
write.csv(as.data.frame(null_dist), file.path(data_dir, "infer_v4_null_reference.csv"), row.names = FALSE)
write.csv(as.data.frame(p_value), file.path(data_dir, "infer_v4_null_p_value.csv"), row.names = FALSE)

native_null_plot <- infer::visualize(null_dist, bins = 25) +
  infer::shade_p_value(obs_stat = obs_stat, direction = "greater")
null_plot <- native_null_plot +
  ggplot2::labs(title = "infer：H0:p=0.45 的原假设分布", x = "每轮达标率", y = "轮数") +
  ggplot2::annotate("text", x = .65, y = Inf, label = "观察率 0.65", vjust = 1.5, color = "red2", family = font_family) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::coord_cartesian(xlim = c(.25, .70)) +
  ggplot2::theme_minimal(base_size = 16, base_family = font_family)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggplot2::ggsave(file.path(article_dir, "文章配图", "02_infer原假设_v4.png"), null_plot, width = 7, height = 4.8, dpi = 170, device = ragg::agg_png)
} else {
  ggplot2::ggsave(file.path(article_dir, "文章配图", "02_infer原假设_v4.png"), null_plot, width = 7, height = 4.8, dpi = 170)
}

# 6. 单组记录换序：每轮仍保留117个yes与63个no ----
# 与正文一致，用sample不放回换序；两组置换检验见permutation_walkthrough.R。
set.seed(20261110)
shuffle_rates <- replicate(
  3, mean(sample(dat$outcome, replace = FALSE) == "yes")
)
cat("\nSingle-group shuffle rates:\n")
print(shuffle_rates)
stopifnot(all(shuffle_rates == obs_stat$stat))
write.csv(
  data.frame(replicate = 1:3, stat = shuffle_rates),
  file.path(data_dir, "infer_v4_shuffle_stats.csv"), row.names = FALSE
)
cat("\nEvidence paths are project-relative: data/infer_v4_* and results/infer_v4_console.txt\n")
sink()
