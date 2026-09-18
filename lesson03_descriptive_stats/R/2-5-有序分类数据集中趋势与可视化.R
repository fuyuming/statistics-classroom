# 2026-09-18｜第二次课·统计描述：有序分类变量的特征与可视化（构成比、累计比例、中位数）
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-5-有序分类数据集中趋势与可视化.R
# 运行：在RStudio打开本文件后Source，或命令行Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("showtext", "haven", "tidyverse", "table1", "DescTools", "scales", "flextable", "gt", "eoffice", "rvg", "patchwork"), repos="https://cloud.r-project.org")
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
  packages <- c("showtext", "haven", "tidyverse", "table1", "DescTools", "scales", "flextable", "gt", "eoffice", "rvg", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) {
    message("首次运行，正在安装缺失的包：", paste(missing,collapse=", "))
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  failed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(failed)) stop("以下包尚未安装成功，请检查网络及安装日志：",paste(failed,collapse=", "))
})

#####============================================================#####
#####  环境设置与基础配置
#####============================================================#####

options(digits = 3)
# 设置数字的显示精度为 3 位，避免过多小数位影响阅读
# 这在处理百分比和比例数据时特别有用

Sys.setenv(LANGUAGE = "en")
# 设置系统环境变量 LANGUAGE 为 "en"，目的是让 R 在报错信息等情况下显示英文内容
# 便于理解错误信息，也方便与其他英文资料对照学习

options(scipen = 200)
# 设置科学计数法的显示规则。scipen是科学计数法的显示惩罚参数；值大时更偏向普通表示，并非200位限制
# 例如：显示 0.001 而不是 1e-03，提高数据的可读性

# 防止绘图中文出现乱码的解决方案
# install.packages("showtext")  # 修正：packages拼写，如果未安装则取消注释运行
# showtext包专门用于在 R 图形中正确显示中文等非ASCII字体，避免出现乱码问题

library(showtext)
# 加载 showtext 包，以便在后续代码中使用其字体渲染功能

showtext_auto()
# 开启自动使用 showtext 功能，确保在所有图形输出中都能正确显示字体
# 这个函数会自动检测并应用合适的字体渲染方式

#####============================================================#####
#####  必需包的检查与加载
#####============================================================#####

# 清理工作环境，确保分析的纯净性
# 原脚本清空环境语句已停用：请用新的R会话运行，避免删除其他分析对象。

# 定义本次有序分类数据分析需要的包列表
required_packages <- c(
  "haven",         # 读取SPSS、Stata、SAS等统计软件的数据文件
  "tidyverse",     # 数据科学工具集，包含dplyr、ggplot2、tidyr等核心包
  "table1",        # 制作学术级描述性统计表格（三线表）
  "DescTools",     # 描述性统计和数据分析的扩展工具
  "scales"         # 图形标度和标签的格式化工具
)

# 智能包管理：检查并安装缺失的包
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    message(paste("正在安装包:", pkg))
    install.packages(pkg)
    library(pkg, character.only = TRUE)
    message(paste("包", pkg, "安装并加载完成"))
  }
}

#####============================================================#####
#####  数据读取与预处理
#####============================================================#####

# === 读取SPSS数据文件 ===
# 使用haven包读取.sav格式的SPSS数据文件
# 注意：确保工作目录下存在 data 文件夹与 education_2025.csv
# 公开版：教育等级编码以 CSV 形式随仓库提供（与课堂用 02.sav 的取值一致，1=小学 … 5=研究生）
data2 <- data.frame(学历 = as.character(read.csv("./data/education_2025.csv")$education_code))

# 初步查看数据结构,
glimpse(data2)
# 这里显示的是 <chr+lbl>，表示：
# chr = character (字符串类型)
# lbl = labeled (带标签的字符串)

# 1. 查看原始值（字符串）
print(data2$学历)
# 2. 查看标签信息（如果有的话）
attributes(data2$学历)


# === 创建有序因子变量 ===
data2$education_ordered <- ordered(data2$学历,
                                   levels = c("1", "2", "3", "4", "5"),
                                   labels = c("小学", "初中", "高中", "本科", "研究生"))
#再次查看数据结构
glimpse(data2)
# 因子转换后的数据结构 
str(data2$education_ordered)    # 查看转换后的因子结构

# 有序因子的优势：
# 1. 保留了类别间的顺序关系
# 2. 在统计分析中会考虑顺序信息
# 3. 可以计算中位数等位置统计量

#集中趋势
#众数
DescTools::Mode(data2$education_ordered)

#中位数
DescTools::Median(data2$education_ordered)

#####============================================================#####
#####  统计表与统计图
#####============================================================#####
# === 详细频数统计分析 ===
# 使用tidyverse方法进行更详细的频数分析
freq_table <- data2 %>%                                  # 使用管道符开始数据处理
  count(education_ordered) %>%                                         # 计算每个学历水平的频数
  arrange(education_ordered) %>%                                   # 按教育等级顺序排列，累计值才有顺序含义
  mutate(
    freq = n / sum(n),                                   # 计算相对频数（比例）
    percent = round(freq * 100, 2),                      # 计算百分比并保留2位小数
    cum_freq = cumsum(n),                                # 计算累计频数
    cum_percent = round(cumsum(n) / sum(n) * 100, 2)              # 计算累计百分比
  )

#详细频数统计表
freq_table

# === 学术级三线表制作 ===
# # 方法A: 使用table1包快速生成符合学术发表标准的描述性统计表
academic_table <- table1::table1(~education_ordered, data = data2)
academic_table
# 导出图表到Word文档
library(flextable)
# 转换并导出
table1_flex <- t1flex(academic_table)
# === 结果保存开关：默认打开（结果直接存到本文件夹）===
# 只想到控制台看、不生成文件时，先运行 options(medstats.export = FALSE) 再跑本脚本。
if (is.null(getOption("medstats.export"))) options(medstats.export = TRUE)

# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)


if (isTRUE(getOption("medstats.export", FALSE))) try(save_as_docx(table1_flex, path = "./output/R/学历统计表.docx"))
# table1包的优势：
# 1. 自动计算频数、百分比
# 2. 格式符合学术期刊要求
# 3. 可处理连续变量和分类变量
# 4. 支持分组分析


# 方法B: gt包（自定义版）
# === 方法B: gt包（自定义版）===
library(gt)
# 首先计算基础统计
tab <- table(data2$education_ordered)
pct <- prop.table(tab) * 100

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
    title = "Table 1. Distribution of Education Level",
    subtitle = paste("Descriptive Statistics (N =", nrow(data2), ")")
  ) %>%
  cols_label(
    Category = "Education Level",
    Frequency = "Count",
    Percentage = "Percent (%)",
    Cumulative = "Cumulative (%)"
  ) %>%
  # 使用新版本的grand_summary_rows()替代summary_rows()
  grand_summary_rows(
    columns = Frequency,
    fns = list(Total = ~ sum(.)),
    fmt = ~ fmt_number(., decimals = 0)
  ) %>%
  # 格式化数字
  fmt_number(columns = c(Frequency), decimals = 0) %>%
  fmt_number(columns = c(Percentage, Cumulative), decimals = 1) %>%
  # 样式设置 - 表头样式
  tab_style(
    style = list(
      cell_text(weight = "bold"),
      cell_fill(color = "#f0f0f0")
    ),
    locations = cells_column_labels()
  ) %>%
  # 表体居中对齐
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body()
  ) %>%
  # 汇总行样式（修正这里的错误）
  tab_style(
    style = list(
      cell_text(weight = "bold", align = "center"),
      cell_fill(color = "#e6f3ff")
    ),
    locations = cells_grand_summary()  # 使用正确的位置函数
  ) %>%
  # 添加脚注
  tab_footnote(
    footnote = "Data source: Survey data analysis",
    locations = cells_title("subtitle")
  ) %>%
  # 设置表格选项
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12,
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    heading.border.bottom.style = "solid",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid"
  )

# 在RStudio中显示
publication_table

# 保存为 Word 文档
if (isTRUE(getOption("medstats.export", FALSE))) try(gtsave(publication_table, "./output/R/教育水平分布表_gt版.docx"))
# 保存为 PNG 图片
if (isTRUE(getOption("medstats.export", FALSE))) try(gtsave(publication_table, "./output/R/教育水平分布表_gt版.png"))



#####============================================================#####
#####  数据可视化：频数与频率图
#####============================================================#####

# === 频数直条图 ===
# 创建专业的频数分布图
freq_plot <- ggplot(data = freq_table, aes(x = education_ordered, y = n)) +
  geom_col(
    fill = "steelblue",                                  # 专业的蓝色
    alpha = 0.8,                                         # 适度透明
    width = 0.7                                          # 柱子宽度
  ) +
  geom_text(
    aes(label = n),                                      # 在柱子上显示频数
    vjust = -0.5,                                        # 文本位置调整
    size = 3.5,                                          # 文本大小
    color = "black"
  ) +
  labs(
    title = "研究对象学历分布",
    subtitle = paste0("样本量 N = ", sum(freq_table$n)),
    x = "学历水平",
    y = "频数",
    caption = "数据来源：研究样本"
  ) +
  scale_y_continuous(
    limits = c(0, max(freq_table$n) * 1.2),             # 动态设置y轴上限
    breaks = seq(0, max(freq_table$n), by = 1),         # 动态设置刻度
    expand = expansion(mult = c(0, 0.1))                 # 调整y轴扩展
  ) +
  theme_minimal() +                                      # 使用简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0),               # x轴标签不旋转（学历标签较短）
    panel.grid.minor = element_blank(),                  # 移除次要网格线
    panel.grid.major.x = element_blank()                 # 移除x轴主网格线
  )

freq_plot

# === 频率直条图（百分比图）===
# 创建频率分布图
freq_rate_plot <- ggplot(data = freq_table, aes(x = education_ordered, y = freq)) +
  geom_col(
    fill = "lightcoral", 
    alpha = 0.8,
    width = 0.7
  ) +
  geom_text(
    aes(label = paste0(percent, "%")),                   # 显示百分比标签
    vjust = -0.5,
    size = 3.5,
    color = "black"
  ) +
  labs(
    title = "研究对象学历分布（频率）",
    subtitle = paste0("样本量 N = ", sum(freq_table$n)),
    x = "学历水平",
    y = "频率",
    caption = "数据来源：研究样本"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::percent_format(accuracy = 1),       # 使用scales包格式化百分比
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

freq_rate_plot

# === 导出到PowerPoint和Adobe Illustrator ===
# 方法1: 导出到PowerPoint
library(eoffice)    #直接使用 eoffice 包导出图形至于PPTX内部
library(rvg)       # 用于矢量图形

# 直接导出单个图形到PowerPoint
if (isTRUE(getOption("medstats.export", FALSE))) eoffice::topptx(freq_plot, file = "./output/R/频数分布图.pptx", 
          width = 10, height = 6)

# 方法2: 导出为SVG格式（可在Adobe Illustrator中编辑）
# SVG是矢量格式，在AI中可以完全编辑
if (isTRUE(getOption("medstats.export", FALSE))) ggsave("./output/R/频数分布图.svg", plot = freq_plot, 
       width = 10, height = 6, units = "in", dpi = 300)


# === 使用 patchwork 拼接图形 ===
# 安装和加载 patchwork
# install.packages("patchwork")
library(patchwork)

# 方法1A: 水平拼接（并排显示）
combined_plot1 <- freq_plot + freq_rate_plot
combined_plot1

# 方法1B: 垂直拼接（上下显示）
combined_plot2 <- freq_plot / freq_rate_plot
combined_plot2

# 方法1C: 添加整体标题和标签（修正版）
combined_plot3 <- freq_plot + freq_rate_plot + 
  plot_annotation(
    title = "Figure 1. 学历分布分析",
    subtitle = paste0("频数和百分比分布 (N = ", sum(freq_table$n), ")"),
    caption = "数据来源：研究样本",
    tag_levels = 'A'  # 添加子图标签 A, B
  ) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

combined_plot3


# === 导出拼接图 ===
# 导出为矢量格式
if (isTRUE(getOption("medstats.export", FALSE))) ggsave("./output/R/学历分布综合图.svg", plot = combined_plot3, 
       width = 12, height = 6, units = "in", dpi = 300)
# 导出为高分辨率PNG）
if (isTRUE(getOption("medstats.export", FALSE))) ggsave("./output/R/学历分布综合图.png", plot = combined_plot3, 
       width = 12, height = 6, units = "in", dpi = 300)

