# 统计课堂42：三种支架设计的残差分布图（UTF-8）
# 一、从项目根目录读数据；也可在RStudio打开本课Rproj后逐节运行。 ----
options(stringsAsFactors = FALSE, digits = 17)
project_root <- getwd()
if (basename(project_root) == "代码") project_root <- dirname(project_root)
source(file.path(project_root, "代码", "residplot.R"))
data <- read.csv(file.path(project_root, "数据", "anova_example.csv"))
print(head(data))
save_outputs <- identical(Sys.getenv("RESIDPLOT_SAVE_OUTPUTS"), "1")
output_dir <- file.path(project_root, "运行结果", "R")
if (save_outputs) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_case <- function(case_data, case_name) {
  # 分类组别拟合的是每个组的均数。
  fit <- lm(value ~ group, data = case_data)
  if (save_outputs) {
    png(file.path(output_dir, paste0(case_name, "_residplot.png")), width = 1200, height = 800)
  }
  result <- residplot(fit)
  if (save_outputs) dev.off()
  z <- result$studentized_residuals
  summary <- data.frame(n = length(z), df_resid = fit$df.residual,
                        bandwidth = result$bandwidth, mean_studentized = mean(z),
                        sd_studentized = sd(z), max_abs_studentized = max(abs(z)))
  cat(case_name, "\n")
  print(summary, row.names = FALSE)
  if (save_outputs) {
    residual_table <- data.frame(case_data, fitted = fitted(fit), residual = residuals(fit),
                                 leverage = hatvalues(fit), studentized = z)
    tables <- list(residuals = residual_table, histogram = result$histogram,
                   curves = result$curves, summary = summary)
    for (name in names(tables)) {
      write.csv(tables[[name]], file.path(output_dir, paste0(case_name, "_", name, ".csv")),
                row.names = FALSE, quote = FALSE)
    }
  }
  invisible(list(fit = fit, plot = result))
}
# 二、原例：直方图、拟合正态线、KDE和底部短线分别展示什么？ ----
base <- run_case(data, "base")
# 三、改动练习：只把A3第5个值改成44，另拟合一个模型。 ----
outlier_data <- data
selected <- outlier_data$group == "A3" & outlier_data$replicate == 5
outlier_data$value[selected] <- 44
outlier <- run_case(outlier_data, "outlier")
cat("比较两张图：哪一处残差位置最突出？最大绝对值由数据重新计算。\n")
