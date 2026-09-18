# 2026-09-18｜第二次课·统计描述：分类无序变量的特征与可视化（频数、构成比、众数）
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-4-分类数据集中趋势与可视化.R
# 运行：在RStudio打开本文件后Source，或命令行Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("showtext", "tidyverse", "DescTools", "table1", "gt", "readxl", "eoffice", "flextable", "rvg", "patchwork"), repos="https://cloud.r-project.org")
# options(digits=3)设置有效数字显示，不改变原数据；统计解释要保留单位和分母。
# 结果默认保存到本文件夹；只想到控制台看时先运行 options(medstats.export=FALSE)。
local({
  # 优先定位被source的文件，其次Rscript参数和RStudio当前编辑文件。
  paths <- Filter(function(x) is.character(x) && length(x)==1 && nzchar(x),
                  lapply(sys.frames(), function(x) x$ofile))
  args <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  path <- if(length(paths)) tail(paths,1)[[1]] else if(length(args)) sub("^--file=", "", args[1]) else ""
  if (!nzchar(path) && requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable())
    path <- rstudioapi::getActiveDocumentContext()$path
  if(nzchar(path)) {
    root <- dirname(normalizePath(path, mustWork=TRUE))
    if(basename(root)=="R" && file.exists(file.path(dirname(root),"医学统计学.Rproj"))) root <- dirname(root)
    setwd(root)
  }
  packages <- c("showtext", "tidyverse", "DescTools", "table1", "gt", "readxl", "eoffice", "flextable", "rvg", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) {
    message("首次运行，正在安装缺失的包：", paste(missing,collapse=", "))
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  failed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(failed)) stop("以下包尚未安装成功，请检查网络及安装日志：",paste(failed,collapse=", "))
})

###### 环境设置#####
options(digits = 3)
# 设置数字的显示精度为 3 位。
Sys.setenv(LANGUAGE ="en")
# 设置系统环境变量 LANGUAGE 为 "en"，目的是让 R 在报错信息等情况下显示英文内容，方便理解和与其他英文资料对照。
options(scipen = 200)
# 设置科学计数法的显示规则。scipen是科学计数法的显示惩罚参数；值大时更偏向普通表示，并非200位限制

# 防止中文出现乱码
#install.packages("showtext")
# 如果未安装 showtext 包，运行这行代码安装。该包用于在 R 图形中正确显示中文等字体，避免出现乱码问题。
library(showtext)
# 加载 showtext 包，以便在后续代码中使用其功能。

showtext_auto()
# 开启自动使用 showtext 功能，确保在图形中可以正确显示字体。

# 清理工作环境
# 原脚本清空环境语句已停用：请用新的R会话运行，避免删除其他分析对象。


###### 加载必需包 #####
# 包的用途说明：
# tidyverse: 数据处理和可视化工具集（包含dplyr, ggplot2等）
# DescTools: 描述性统计分析工具
# table1: 制作学术论文标准表格
# gt: 制作美观的统计表格
# xlsx: 读取Excel文件
# eoffice: 导出图表到Office文档

required_packages <- c("tidyverse", "DescTools", "table1", "gt", "readxl", "eoffice")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}


#关于众数
#
color_1 <- c("赤","橙","黄","绿","青","蓝","紫")
DescTools::Mode(color_1)

color_2 <- rep(x=c("赤","橙","黄","绿","青","蓝","紫"),
               times=c(2,6,1,10,3,10,4))
color_2
DescTools::Mode(color_2)


###### 数据读取与预处理 #####
# 读取数据
# 公开版：性别记录以 CSV 形式随仓库提供（与课堂用 01-.xlsx 数值一致）
data1 <- data.frame(X1 = read.csv("./data/sex_2025.csv")$sex)

# 初步查看数据结构
glimpse(data1) # 发现数据变量名称是X1，数据类型是chr，不符合分类特性，需要改变量名字和改变数据类型为因子型

# 预处理：重命名并设置因子
data1 <- dplyr::rename(data1, sex = X1)
data1$sex <- factor(data1$sex, levels = c("男","女"))

# 再次查看数据结构
glimpse(data1)


###### 1. 基础频数分析 #####
# 频数分布
tab <- table(data1$sex)
tab
# 百分比分布
pct <- prop.table(tab) * 100
round(pct, 1)
# 累计频数
cumsum(tab)
# 累计百分比
round(cumsum(pct), 1)


###### 2. DescTools 快速摘要 #####
# 生成描述性统计摘要和图表
Desc(data1$sex)
# 导出图表到Word文档
# === 结果保存开关：默认打开（结果直接存到本文件夹）===
# 只想到控制台看、不生成文件时，先运行 options(medstats.export = FALSE) 再跑本脚本。
if (is.null(getOption("medstats.export"))) options(medstats.export = TRUE)

# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)


if (isTRUE(getOption("medstats.export", FALSE))) try(topptx(ggplot(data1, aes(x = sex)) + geom_bar(), filename = "./output/R/性别分布图.pptx", width = 8, height = 6))

#发现集中趋势的众数
DescTools::Mode(data1$sex)

###### 3. 学术化表格制作 #####
# 方法A: table1包（简洁版）
table1_result <- table1(~ sex, data = data1)
# 使用 table1 包中的 table1 函数生成三线表。
# 这里的参数 ~sex 表示描述sex这一变量；竖线后才指定分层变量，如~age | sex。
# data = data1 指定了用于制表的数据框为 data1。
# table1 函数通常会生成一个简洁的三线表，包含标题行、分隔线和数据行。
table1_result
# 导出图表到Word文档
library(flextable)
# 转换并导出
table1_flex <- t1flex(table1_result)
if (isTRUE(getOption("medstats.export", FALSE))) try(save_as_docx(table1_flex, path = "./output/R/性别分布表.docx"))


# 方法B: gt包（自定义版）
# 准备数据
freq_summary <- data.frame(
  Category = names(tab),
  Frequency = as.numeric(tab),
  Percentage = round(as.numeric(pct), 1),
  Cumulative = round(cumsum(as.numeric(pct)), 1)
)

# 创建gt表格
publication_table <- freq_summary %>%
  gt() %>%
  tab_header(
    title = "Table 1. Distribution of Gender",
    subtitle = "Descriptive Statistics (N = 12)"
  ) %>%
  cols_label(
    Category = "Gender",
    Frequency = "Count",
    Percentage = "Percent (%)",
    Cumulative = "Cumulative (%)"
  ) %>%
  # 添加总计行
  summary_rows(
    groups = NULL,
    columns = Frequency,
    fns = list(Total = ~ sum(.)),
    fmt = ~ fmt_number(., decimals = 0)
  ) %>%
  # 格式化数字
  fmt_number(columns = c(Frequency), decimals = 0) %>%
  fmt_number(columns = c(Percentage, Cumulative), decimals = 1) %>%
  # 样式设置
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "#f0f0f0")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) %>%
  # 添加脚注
  tab_footnote(
    footnote = "Data source: data/sex_2025.csv",
    locations = cells_title("subtitle")
  )

#在Rstudio中显示
print(publication_table)
# 保存为 Word 文档
if (isTRUE(getOption("medstats.export", FALSE))) try(gtsave(publication_table, "./output/R/性别分布表_gt版.docx"))



###### 4. 数据可视化 #####

# 直条图(频数图)
p1 <- ggplot(data1, aes(x = sex)) +
  geom_bar(fill = "steelblue", alpha = 0.7, color = "black") +
  geom_text(stat = 'count', aes(label = after_stat(count)), 
            vjust = -0.5, size = 4) +
  labs(
    title = "Gender Distribution - Frequency",
    x = "Gender", 
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

p1

# 百分比条形图(频率图)
p2 <- ggplot(data1, aes(x = sex)) +
  geom_bar(aes(y = after_stat(count/sum(count))),
           fill = "lightcoral", alpha = 0.7, color = "black") +
  geom_text(aes(y = after_stat(count/sum(count)), 
                label = paste0(round(after_stat(count/sum(count))*100, 1), "%")),
            stat = 'count', vjust = -0.5, size = 4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Gender Distribution - Percentage",
    x = "Gender", 
    y = "Percentage"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )

p2

# === 导出到PowerPoint和Adobe Illustrator ===
# 方法1: 导出到PowerPoint
library(eoffice)    #直接使用 eoffice 包导出图形至于PPTX内部
library(rvg)       # 用于矢量图形

# 直接导出单个图形到PowerPoint
if (isTRUE(getOption("medstats.export", FALSE))) eoffice::topptx(p1, file = "./output/R/分类直条图.pptx", 
                width = 10, height = 6)

# 方法2: 导出为SVG格式（可在Adobe Illustrator中编辑）
# SVG是矢量格式，在AI中可以完全编辑
if (isTRUE(getOption("medstats.export", FALSE))) ggsave("./output/R/分类直条图.svg", plot = p1, 
       width = 10, height = 6, units = "in", dpi = 300)

###### 使用 patchwork 拼接图形 #####
# 安装和加载 patchwork
# install.packages("patchwork")
library(patchwork)

# 方法1A: 水平拼接（并排显示）
combined_plot1 <- p1 + p2
combined_plot1

# 方法1B: 垂直拼接（上下显示）
combined_plot2 <- p1 / p2
combined_plot2

# 方法1C: 添加整体标题和标签
combined_plot3 <- p1 + p2 + 
  plot_annotation(
    title = "Figure 1. Gender Distribution Analysis",
    subtitle = "Frequency and Percentage Distribution (N = 12)",
    caption = "Data source: 01-.xlsx file",
    tag_levels = 'A'  # 添加子图标签 A, B
  ) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

combined_plot3



