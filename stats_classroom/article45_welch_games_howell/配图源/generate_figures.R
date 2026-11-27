# 公众号正文两幅数据图统一由R绘制；公式PNG由generate_formulas.py排版。
# 命令行：Rscript 配图源/generate_figures.R；或在本交付根目录Source本文件。
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg)) {
  script_file <- sub("^--file=", "", script_arg[1])
  root <- dirname(dirname(normalizePath(script_file)))
} else {
  root <- getOption("statistics45.root", getwd())
}
source(file.path(root, "练习项目", "R", "plot_article.R"), encoding = "UTF-8")
draw_article(data_file = file.path(root, "练习项目", "data", "hetero_normal.csv"),
             out_dir = file.path(root, "文章配图"),
             coordinates_file = file.path(root, "配图源", "区间坐标.csv"))
