# 统计课堂26：两组率差的 infer 置换检验（新手可运行版）
# 目标：在 H0“分组与结局独立”下，观察 A-B 达标率差是否足够大。

local({
  # 保存运行日志，并在脚本退出时关闭输出重定向。
  project_root <- getwd()
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
    project_root <- normalizePath(file.path(dirname(script_path), ".."))
  }
  data_dir <- file.path(project_root, "data")
  out_dir <- file.path(project_root, "results")
  article_dir <- dirname(project_root)
  if (!file.exists(file.path(data_dir, "cases.csv"))) {
    stop("请从本项目运行；未找到 data/cases.csv。")
  }
  if (!requireNamespace("infer", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("本教程需要 infer 与 ggplot2。")
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(article_dir, "文章配图"), showWarnings = FALSE, recursive = TRUE)

  console_path <- file.path(out_dir, "permutation_v10_console.txt")
  sink(console_path, split = TRUE)
  on.exit({
    while (sink.number() > 0) sink()
  }, add = TRUE)

  options(digits = 8)
  cat("R:", R.version.string, "\n")
  cat("infer:", as.character(packageVersion("infer")), "\n")
  cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
  if (requireNamespace("ragg", quietly = TRUE)) {
    cat("ragg:", as.character(packageVersion("ragg")), "\n")
  } else {
    cat("ragg: not installed; using default ggplot2 device\n")
  }
  cat("Seed used for each independent block: 20261110\n")

  # 8 名独立受试者：A 组 3/4 达标，B 组 1/4 达标。
  study <- data.frame(
    group = factor(c(rep("A", 4), rep("B", 4)), levels = c("A", "B")),
    outcome = factor(c("yes", "yes", "yes", "no", "yes", "no", "no", "no"),
                     levels = c("no", "yes"))
  )
  cat("\nStudy (8 rows):\n")
  print(study, row.names = FALSE)
  cat("\nStudy table: group x outcome\n")
  print(with(study, table(group, outcome)))
  write.csv(study, file.path(data_dir, "permutation_v10_raw8.csv"), row.names = FALSE)

  # 观察到的 A-B 差：3/4 - 1/4 = 0.5。
  obs_diff <- study |>
    infer::specify(outcome ~ group, success = "yes") |>
    infer::calculate(stat = "diff in props", order = c("A", "B"))
  cat("\nObserved difference (A-B):\n")
  print(as.data.frame(obs_diff), row.names = FALSE)

  # H0：分组与结局独立；置换 outcome 在固定 group 名额内的标签。
  set.seed(20261110)
  perm_spec <- study |>
    infer::specify(outcome ~ group, success = "yes") |>
    infer::hypothesize(null = "independence")
  demo_samples <- perm_spec |>
    infer::generate(reps = 5, type = "permute")
  demo_stats <- demo_samples |>
    infer::calculate(stat = "diff in props", order = c("A", "B"))
  cat("\nFive-round permute sample details:\n")
  cat("nrow(demo_samples):", nrow(demo_samples), "\n")
  cat("replicate x group counts (each cell should be 4):\n")
  print(with(demo_samples, table(replicate, group)))
  cat("Five-round calculated differences:\n")
  print(as.data.frame(demo_stats), row.names = FALSE)
  write.csv(as.data.frame(demo_samples), file.path(data_dir, "permutation_v10_demo_samples.csv"), row.names = FALSE)
  write.csv(as.data.frame(demo_stats), file.path(data_dir, "permutation_v10_demo_stats.csv"), row.names = FALSE)

  # 10000 轮 Monte Carlo 置换分布；每个独立随机块都重新固定种子。
  set.seed(20261110)
  perm_dist <- perm_spec |>
    infer::generate(reps = 10000, type = "permute") |>
    infer::calculate(stat = "diff in props", order = c("A", "B"))
  cat("\n10000-round permutation distribution:\n")
  cat("nrow(perm_dist):", nrow(perm_dist), "\n")
  print(summary(perm_dist$stat))
  perm_hits <- sum(perm_dist$stat >= as.numeric(obs_diff$stat))
  perm_p <- infer::get_p_value(perm_dist, obs_stat = obs_diff, direction = "greater")
  cat("Actual hits (stat >= 0.5):", perm_hits, "of", nrow(perm_dist), "\n")
  cat("Actual hit proportion:", perm_hits / nrow(perm_dist), "\n")
  cat("get_p_value(direction='greater'):\n")
  print(as.data.frame(perm_p), row.names = FALSE)
  write.csv(as.data.frame(perm_dist), file.path(data_dir, "permutation_v10_perm_dist.csv"), row.names = FALSE)
  write.csv(as.data.frame(perm_p), file.path(data_dir, "permutation_v10_p_value.csv"), row.names = FALSE)

  # 原生 infer 图层保留；最终图用百分点横轴与离散质量整柱标记尾部。
  available_families <- character()
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    available_families <- systemfonts::system_fonts()$family
  }
  font_family <- if ("Heiti SC" %in% available_families) {
    "Heiti SC"
  } else if ("Noto Sans CJK SC" %in% available_families) {
    "Noto Sans CJK SC"
  } else if ("PingFang SC" %in% available_families) {
    "PingFang SC"
  } else {
    "sans"
  }
  cat("Plot font family:", font_family, "\n")
  base_plot <- infer::visualize(perm_dist, bins = 5)
  # 保留原生 shade 图层对象，供课堂代码对照；最终图改用离散质量整柱标记，
  # 避免连续面积在 .5 处只染半根柱而误读尾部计数。
  native_plot <- base_plot +
    infer::shade_p_value(obs_stat = obs_diff, direction = "greater")
  stat_freq <- as.data.frame(table(perm_dist$stat), stringsAsFactors = FALSE)
  names(stat_freq) <- c("stat", "n")
  stat_freq$stat <- as.numeric(stat_freq$stat)
  stat_freq$tail <- ifelse(stat_freq$stat >= as.numeric(obs_diff$stat), "≥50百分点", "其余")
  perm_plot <- base_plot +
    ggplot2::geom_col(
      data = stat_freq,
      mapping = ggplot2::aes(x = stat, y = n, fill = tail),
      width = .45, alpha = 1, inherit.aes = FALSE
    ) +
    ggplot2::geom_vline(xintercept = as.numeric(obs_diff$stat),
                        color = "red2", linewidth = .9) +
    ggplot2::labs(title = "infer：两组率差的置换分布",
                  x = "每轮率差（A−B，百分点）", y = "轮数", fill = "率差") +
    ggplot2::scale_x_continuous(breaks = seq(-1, 1, .5), labels = function(x) x * 100) +
    ggplot2::scale_fill_manual(values = c("其余" = "grey65", "≥50百分点" = "red2")) +
    ggplot2::theme_minimal(base_size = 16, base_family = font_family)
  plot_path <- file.path(article_dir, "文章配图", "03_置换分布_v10.png")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(plot_path, perm_plot, width = 7, height = 4.8, dpi = 170, device = ragg::agg_png)
  } else {
    ggplot2::ggsave(plot_path, perm_plot, width = 7, height = 4.8, dpi = 170)
  }
  cat("Plot saved:", plot_path, "\n")

  # 70 个等可能 A 组组合的精确置换分布，作为 Monte Carlo 的可核对基准。
  y <- as.integer(study$outcome == "yes")
  exact_idx <- combn(seq_len(nrow(study)), 4)
  exact_stat <- apply(exact_idx, 2, function(idx) mean(y[idx]) - mean(y[-idx]))
  exact <- data.frame(permutation = seq_along(exact_stat), stat = exact_stat)
  cat("\nExact enumeration:\n")
  cat("nrow(exact):", nrow(exact), "\n")
  print(table(exact$stat))
  exact_hits <- sum(exact$stat >= as.numeric(obs_diff$stat))
  cat("Exact hits (stat >= 0.5):", exact_hits, "of", nrow(exact), "\n")
  cat("Exact P:", exact_hits / nrow(exact), "\n")
  write.csv(exact, file.path(data_dir, "permutation_v10_exact70.csv"), row.names = FALSE)

  # 最终自动检查，便于父代理和读者复核。
  stopifnot(nrow(demo_samples) == 40,
            all(as.vector(with(demo_samples, table(replicate, group))) == 4),
            nrow(perm_dist) == 10000,
            as.numeric(perm_p$p_value) == perm_hits / nrow(perm_dist),
            exact_hits == 17,
            identical(as.integer(table(exact$stat)), c(1L, 16L, 36L, 16L, 1L)))
  cat("\nCHECK: PASS (demo cells=4; MC P=hits/10000; exact hits=17/70; exact frequencies=1,16,36,16,1)\n")
})
