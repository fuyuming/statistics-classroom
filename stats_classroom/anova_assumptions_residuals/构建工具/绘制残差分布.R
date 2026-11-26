# 从正文共享数据重建residplot辅助图；在任意目录用Rscript运行本文件。
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]))
article_root <- dirname(dirname(script_file))
project_root <- file.path(article_root, "配套项目")
source(file.path(project_root, "代码", "residplot.R"))
dat <- read.csv(file.path(project_root, "数据", "residual_cases.csv"))
fit_a <- lm(value ~ group, data = subset(dat, scenario == "A"))
png(file.path(article_root, "文章配图", "04_残差分布.png"),
    width = 1800, height = 1200, res = 220, pointsize = 20, type = "cairo-png")
# 缩小图例，保留轴与标题字号，使图例不遮挡密度峰或柱形。
par(mar = c(5.8, 5.8, 3.0, 1.2), mgp = c(3.2, 1.0, 0), cex = 0.72,
    cex.axis = 1.1 / 0.72, cex.lab = 1.25 / 0.72, cex.main = 1.2 / 0.72)
residplot(fit_a, nbreaks = 5L, draw = TRUE)
dev.off()
