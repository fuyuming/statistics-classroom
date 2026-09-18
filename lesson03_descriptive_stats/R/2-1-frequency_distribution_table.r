# 2026-09-18｜第二次课·统计描述：分类与连续数据的频数表（例2-1 诊室问卷；138 例红细胞）
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-1-frequency_distribution_table.r
# 运行：在RStudio打开本文件后Source，或命令行Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("showtext", "sjPlot", "gt", "tidyverse", "DescTools", "vcd", "readxl"), repos="https://cloud.r-project.org")
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
  packages <- c("showtext", "sjPlot", "gt", "tidyverse", "DescTools", "vcd", "readxl")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) {
    message("首次运行，正在安装缺失的包：", paste(missing,collapse=", "))
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  failed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(failed)) stop("以下包尚未安装成功，请检查网络及安装日志：",paste(failed,collapse=", "))
})

### 环境设置#####
options(digits = 7)
# 设置数字的显示精度为 7 位
Sys.setenv(LANGUAGE ="en")
# 设置系统环境变量 LANGUAGE 为 "en"，让 R 在报错信息等情况下显示英文内容
options(scipen = 200)
# 设置科学计数法的显示规则。scipen是科学计数法的显示惩罚参数；值大时更偏向普通表示，并非200位限制

# 防止中文出现乱码
#install.packages("showtext")  # 修正：packages拼写错误
# 如果未安装 showtext 包，运行这行代码安装。该包用于在 R 图形中正确显示中文等字体
library(showtext)
# 加载 showtext 包
showtext_auto()
# 开启自动使用 showtext 功能，确保在图形中可以正确显示字体

# ============================================================================
# 频数表的制作示例
# ============================================================================

#####============================================================#####
#####  第一部分：数据准备与包加载
#####============================================================#####

# 清理工作环境
# 原脚本清空环境语句已停用：请用新的R会话运行，避免删除其他分析对象。

# 加载必需的包
required_packages <- c(
  "sjPlot",        # 制作发表级统计表格和图形
  "gt",            # 制作精美的表格
  "tidyverse",     # 数据处理和可视化工具集（包含dplyr, ggplot2等）
  "DescTools",     # 描述性统计和数据分析工具
  "vcd",           # 可视化分类数据（Visualizing Categorical Data）
  "readxl"         # 读取Excel文件
)

# 检查并安装缺失的包
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}



#####============================================================#####
#####  第二部分：类别数据频数表制作
#####============================================================#####
# 读取分类数据
example2_2 <- read.csv("./example2_2.csv", stringsAsFactors = TRUE)

# 查看数据结构
str(example2_2)
head(example2_2)


# === 2.1 一维频数表 ===

# 制作简单一维频数表
tab1 <- table(example2_2$诊室)
tab1  # 诊室频数分布

# 转换为百分比
tab1_percent <- prop.table(tab1) * 100
round(tab1_percent, 2)  # 诊室百分比分布

# === 2.2 二维列联表 ===

# 制作态度和诊室的二维列联表
attach(example2_2)  # 附加数据框，方便直接使用变量名

tab2 <- table(态度, 诊室)
tab2  # 态度 × 诊室 列联表

# 添加边际合计
addmargins(tab2)  # 带边际合计的列联表

# 转换为百分比
round(addmargins(prop.table(tab2) * 100), 2)  # 百分比列联表

detach(example2_2)  # 取消附加

# === 2.3 多维列联表 ===

# 使用ftable函数制作三维列联表
tab3 <- ftable(example2_2, 
               row.vars = c("性别", "态度"), 
               col.vars = "诊室")
tab3  # 三维列联表（性别+态度 × 诊室）

# 带边际合计的三维列联表
ftable(addmargins(table(example2_2$性别, 
                        example2_2$态度, 
                        example2_2$诊室)))

# 另一种排列方式：行变量为"诊室"，列变量为"性别"、"态度"
ftable(example2_2, 
       row.vars = c("诊室"), 
       col.vars = c("性别", "态度"))

# 使用vcd包的structable函数（公式语法）
structable(性别 + 态度 ~ 诊室, data = example2_2)

# === 2.5 使用sjPlot制作学术级二维列联表 ===

# 制作发表级二维列联表
contingency_table_sjplot <- sjPlot::tab_xtab(
  var.row = example2_2$性别,
  var.col = example2_2$态度,
  title = "Table 1. Cross-tabulation of Gender and Attitude",
  var.labels = c("Gender", "Attitude"),
  show.row.prc = TRUE,    # 显示行百分比
  show.col.prc = TRUE,    # 显示列百分比
  show.summary = TRUE,    # 显示卡方检验结果
  encoding = "UTF-8"
)

contingency_table_sjplot  # 显示表格（在RStudio中会在Viewer面板显示）


#####============================================================#####
#####  gt包处理高维度列联表的完整方案
#####============================================================#####

# 加载gt表格制作相关的包
library(gt)          # 制作精美的发表级表格
library(dplyr)       # 数据操作和管道操作符
library(tidyr)       # 数据整理，包括pivot_wider/longer等函数
library(purrr)       # 函数式编程工具，用于map系列函数
library(stringr)     # 字符串处理工具

# === 方法1：gt的分组功能处理三维列联表 ===

# 1.1 准备三维数据
three_way_data <- example2_2 %>%
  # 计算每个组合的频数
  count(诊室, 性别, 态度) %>%
  # 确保所有组合都存在，缺失的组合用0填充
  complete(诊室, 性别, 态度, fill = list(n = 0)) %>%
  # 按诊室和性别分组，计算百分比
  group_by(诊室, 性别) %>%
  mutate(
    total_by_gender = sum(n),                              # 计算每个性别在该诊室的总数
    row_pct = if_else(total_by_gender > 0,                 # 计算行百分比
                      round(n / total_by_gender * 100, 1), 
                      0),
    display = paste0(n, " (", row_pct, "%)")               # 创建显示格式：频数(百分比)
  ) %>%
  ungroup() %>%
  # 转换为宽格式，态度作为列
  dplyr::select(诊室, 性别, 态度, display) %>%
  pivot_wider(names_from = 态度, values_from = display, values_fill = "0 (0.0%)") %>%
  # 添加总计列
  left_join(
    example2_2 %>%
      count(诊室, 性别) %>%                                # 计算每个诊室-性别组合的总数
      mutate(Total = paste0(n, " (100.0%)")),             # 总计始终是100%
    by = c("诊室", "性别")
  ) %>%
  # 创建分组标识，用于gt表格的分组显示
  mutate(
    诊室_label = paste("诊室", 诊室, "- 性别与态度分布"),
    .before = 1                                            # 将新列放在第一列
  ) %>%
  dplyr::select(-诊室, -n)                                        # 移除不需要的列

# 1.2 创建三维gt表格
three_way_gt <- three_way_data %>%
  gt(
    groupname_col = "诊室_label",                          # 指定分组列
    rowname_col = "性别"                                   # 指定行名列
  ) %>%
  # 添加表格标题和副标题
  tab_header(
    title = "Table 1. Three-way Cross-tabulation Analysis",
    subtitle = "Gender × Attitude × Department Distribution"
  ) %>%
  # 设置行名标题
  tab_stubhead(label = "Gender") %>%
  # 设置列标题（将中文转为英文）
  cols_label(
    赞成 = "Approve",
    反对 = "Oppose", 
    Total = "Total"
  ) %>%
  # 设置分组标题的样式：深蓝色背景，白色粗体文字
  tab_style(
    style = list(
      cell_fill(color = "#2F5597"),                        # 背景色
      cell_text(color = "white", weight = "bold", size = px(12))  # 文字样式
    ),
    locations = cells_row_groups()                         # 应用到分组行
  ) %>%
  # 设置"赞成"列的样式：浅绿色背景，深绿色粗体文字
  tab_style(
    style = list(
      cell_fill(color = "#E8F5E8"),
      cell_text(weight = "bold", color = "#006100")
    ),
    locations = cells_body(columns = 赞成)                 # 应用到赞成列的所有单元格
  ) %>%
  # 设置"反对"列的样式：浅红色背景，深红色粗体文字
  tab_style(
    style = list(
      cell_fill(color = "#FFF2F2"),
      cell_text(weight = "bold", color = "#9C0006")
    ),
    locations = cells_body(columns = 反对)                 # 应用到反对列的所有单元格
  ) %>%
  # 设置所有列居中对齐
  cols_align(align = "center", columns = everything()) %>%
  # 添加脚注说明
  tab_footnote(
    footnote = "Values shown as count (row percentage within department and gender)",
    locations = cells_title(groups = "subtitle")           # 脚注位置在副标题处
  )

three_way_gt

# === 方法2：嵌套表格结构 ===

# 2.1 创建带小计的嵌套结构
nested_three_way <- example2_2 %>%
  # 计算基础频数
  count(诊室, 性别, 态度) %>%
  complete(诊室, 性别, 态度, fill = list(n = 0)) %>%
  # 计算性别内部的百分比
  group_by(诊室, 性别) %>%
  mutate(
    total_gender = sum(n),                                 # 每个性别在该诊室的总数
    row_pct = round(n / total_gender * 100, 1),           # 性别内部百分比
    display = paste0(n, " (", row_pct, "%)")              # 显示格式
  ) %>%
  ungroup() %>%
  # 转换为宽格式
  dplyr::select(诊室, 性别, 态度, display) %>%
  pivot_wider(names_from = 态度, values_from = display, values_fill = "0 (0.0%)") %>%
  # 添加性别小计
  left_join(
    example2_2 %>%
      count(诊室, 性别) %>%
      mutate(Gender_Total = paste0(n, " (100.0%)")),      # 性别小计
    by = c("诊室", "性别")
  ) %>%
  # 为每个诊室添加总计行
  group_by(诊室) %>%
  group_modify(~ {
    # 计算诊室总计
    dept_total <- example2_2 %>%
      filter(诊室 == .y$诊室) %>%                         # 筛选当前诊室的数据
      count(态度) %>%                                      # 按态度计数
      pivot_wider(names_from = 态度, values_from = n, values_fill = 0) %>%
      mutate(
        total_dept = 赞成 + 反对,                          # 诊室总人数
        赞成 = paste0(赞成, " (", round(赞成/total_dept*100, 1), "%)"),  # 诊室内赞成百分比
        反对 = paste0(反对, " (", round(反对/total_dept*100, 1), "%)"),  # 诊室内反对百分比
        Gender_Total = paste0(total_dept, " (100.0%)"),   # 诊室总计
        性别 = "部门总计"                                  # 标识这是总计行
      ) %>%
      dplyr::select(性别, 赞成, 反对, Gender_Total)
    
    # 将原数据和总计行合并
    bind_rows(.x, dept_total)
  }) %>%
  ungroup() %>%
  # 创建诊室标签
  mutate(
    诊室_label = paste("诊室", 诊室),
    .before = 1
  ) %>%
  dplyr::select(-诊室)

# 2.2 创建嵌套gt表格
nested_gt <- nested_three_way %>%
  gt(
    groupname_col = "诊室_label",                          # 按诊室分组
    rowname_col = "性别"                                   # 性别作为行名
  ) %>%
  # 表格标题
  tab_header(
    title = "Table 2. Nested Three-way Analysis with Department Totals",
    subtitle = "Hierarchical Display of Gender × Attitude × Department"
  ) %>%
  tab_stubhead(label = "Gender") %>%
  # 列标题
  cols_label(
    赞成 = "Approve",
    反对 = "Oppose",
    Gender_Total = "Gender Total"
  ) %>%
  # 分组标题样式：蓝色背景
  tab_style(
    style = list(
      cell_fill(color = "#4472C4"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_row_groups()
  ) %>%
  # 部门总计行样式：浅蓝色背景，深蓝色斜体粗体文字
  tab_style(
    style = list(
      cell_fill(color = "#D9E2F3"),
      cell_text(weight = "bold", style = "italic", color = "#1F4E79")
    ),
    locations = cells_body(rows = 性别 == "部门总计")       # 只应用到总计行
  ) %>%
  # 普通行的"赞成"列样式
  tab_style(
    style = cell_fill(color = "#E8F5E8"),
    locations = cells_body(columns = 赞成, rows = 性别 != "部门总计")  # 排除总计行
  ) %>%
  # 普通行的"反对"列样式
  tab_style(
    style = cell_fill(color = "#FFF2F2"),
    locations = cells_body(columns = 反对, rows = 性别 != "部门总计")  # 排除总计行
  ) %>%
  cols_align(align = "center", columns = everything())

nested_gt



#####============================================================#####
#####  第三部分：连续变量频数表制作
#####============================================================#####

# === 3.1 连续变量分组频数表 ===

# 读取连续变量数据（红细胞数据）
# 公开版：红细胞原始数据以 CSV 形式随仓库提供（与课堂用 xlsx 数值一致）
rbc_df <- data.frame(x = read.csv("./data/rbc.csv")$rbc)

# 查看数据结构和基本信息
str(rbc_df)                                                # 查看数据类型和结构
head(rbc_df)                                               # 查看前6行数据

# 使用DescTools包的Freq函数 - 默认分组
rbc_freq_default <- Freq(rbc_df$x)
rbc_freq_default                                           # 显示红细胞数默认分组频数表

# === 3.2 自定义分组间距（基于统计学原理）===

# 为什么按照0.3作为间距分组？基于Sturges规则计算：

# 第一步：计算基本统计量
n     <- length(rbc_df$x)  ; n                             # 样本量：观测值个数
min_x <- min(rbc_df$x)   ;  min_x                          # 最小值：数据的下界
max_x <- max(rbc_df$x)   ;  max_x                          # 最大值：数据的上界
R     <- max_x - min_x   ; R                               # 数据范围宽度：极差

# 第二步：应用Sturges规则计算理论分组数
# Sturges规则：k = 1 + log2(n)，适用于正态分布或近正态分布的数据
k_sturges_raw <- 1 + log2(n)                              # 理论分组数（可能是小数）
k_sturges     <- round(k_sturges_raw) ; k_sturges         # 四舍五入取整，得到实际分组数k=8

# 第三步：根据分组数计算组距
h_sturges     <- R / k_sturges ;h_sturges                 # 组距 = 极差/分组数，约0.29875

# 第四步：选择合适的组距
# 0.29875非常接近0.3，为了便于理解和计算，选择0.3作为最终组距
h_final <- 0.3

# 第五步：确定分组边界
# 下限：向下取整到0.3的整数倍，确保包含最小值
low_edge  <- floor(min_x / h_final) * h_final ;  low_edge    # 分组下界
# 上限：向上取整到0.3的整数倍，确保包含最大值  
high_edge <- ceiling(max_x / h_final) * h_final ; high_edge  # 分组上界

# 第六步：生成分组断点
# 从下界到上界，以0.3为步长生成断点序列
breaks <- seq(low_edge, high_edge, by = h_final) ; breaks

# 第七步：生成分组标签
# 创建左闭右开区间标签 [a, b)
labels <- paste0("[", head(breaks, -1), "-", tail(breaks, -1), ")") ; labels
# 将最后一个区间改为闭区间 [a, b]，确保包含最大值
labels[length(labels)] <- sub("\\)$", "]", labels[length(labels)])

# 第八步：使用自定义分组创建频数表
rbc_freq_custom <- Freq(
  rbc_df$x,                                                # 要分析的变量
  breaks = breaks,                                         # 自定义分组断点
  right  = FALSE,                                          # 左闭右开区间 [a,b)
  include.lowest = TRUE,                                   # 包含最小端点
  labels = labels                                          # 自定义区间标签
)

rbc_freq_custom                                            # 显示自定义分组频数表

# === 3.3 制作发表级红细胞数频数分布表 ===

# 第一步：数据预处理和分组
rbc_summary <- rbc_df %>%
  mutate(
    # 使用cut函数创建红细胞数分组
    rbc_group = cut(x,                                     # 要分组的连续变量
                    breaks = seq(3.0, 5.7, by = 0.3),     # 分组断点：从3.0到5.7，步长0.3
                    right = FALSE,                         # 左闭右开区间
                    include.lowest = TRUE,                 # 包含最低值
                    labels = c("[3.0-3.3)", "[3.3-3.6)", "[3.6-3.9)", 
                               "[3.9-4.2)", "[4.2-4.5)", "[4.5-4.8)", 
                               "[4.8-5.1)", "[5.1-5.4)", "[5.4-5.7)"))  # 自定义标签
  ) %>%
  filter(!is.na(rbc_group)) %>%                            # 移除缺失值（如果有的话）
  count(rbc_group, name = "Frequency") %>%                 # 计算每组的频数
  mutate(
    # 计算各种统计指标
    Percentage = round(Frequency / sum(Frequency) * 100, 1),        # 百分比：该组频数/总频数×100%
    `Cumulative Frequency` = cumsum(Frequency),                     # 累积频数：当前组及之前所有组的频数之和
    `Cumulative Percentage` = round(cumsum(Frequency) / sum(Frequency) * 100, 1)  # 累积百分比
  ) %>%
  dplyr::rename(`RBC Count (×10¹²/L)` = rbc_group)                # 重命名分组列为更专业的名称

# 第二步：使用gt包制作发表级三线表
publication_rbc_table <- rbc_summary %>%
  gt() %>%                                                 # 创建gt表格对象
  # 添加表格标题和副标题
  tab_header(
    title = "Table 4. Frequency Distribution of Red Blood Cell Count",
    subtitle = "Distribution among 138 normal adult women"
  ) %>%
  # 设置列标题样式：上边框加粗，文字加粗（三线表的第一条线）
  tab_style(
    style = list(
      cell_borders(sides = "top", weight = px(2)),         # 顶部粗边框
      cell_text(weight = "bold")                           # 列标题加粗
    ),
    locations = cells_column_labels()                      # 应用到所有列标题
  ) %>%
  # 设置表格底部边框（三线表的第三条线）
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(2)),
    locations = cells_body(rows = nrow(rbc_summary))       # 应用到最后一行
  ) %>%
  # 设置列标题下边框（三线表的第二条线）
  tab_style(
    style = cell_borders(sides = "bottom", weight = px(1)),
    locations = cells_column_labels()
  ) %>%
  # 设置列对齐：第一列左对齐，其他列居中对齐
  cols_align(align = "center", columns = -1) %>%          # -1表示除第一列外的所有列
  # 添加脚注说明
  tab_footnote(
    footnote = "RBC: Red Blood Cell; Data collected by random sampling",
    locations = cells_title(groups = "subtitle")           # 脚注位置在副标题处
  ) %>%
  # 设置表格字体为学术期刊常用的Times New Roman
  opt_table_font(font = "Times New Roman") %>%
  # 设置表格选项，形成标准三线表效果
  tab_options(
    table.border.top.style = "hidden",                    # 隐藏表格最外层上边框
    table.border.bottom.style = "hidden",                 # 隐藏表格最外层下边框
    column_labels.border.bottom.width = px(1),            # 列标题下边框宽度
    table_body.border.bottom.width = px(2)                # 表格主体底部边框宽度
  )

publication_rbc_table                                      # 显示发表级红细胞数频数分布表

# === 结果导出（可选）：需要时先运行 options(medstats.export=TRUE)，再重跑本段 ===
# 默认只在控制台显示、不落盘；打开开关后生成下面三个文件（落在本文件夹根目录）。
# === 结果保存开关：默认打开（结果直接存到本文件夹）===
# 只想到控制台看、不生成文件时，先运行 options(medstats.export = FALSE) 再跑本脚本。
if (is.null(getOption("medstats.export"))) options(medstats.export = TRUE)

# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)


if (isTRUE(getOption("medstats.export", FALSE))) {
  try(openxlsx::write.xlsx(as.data.frame(rbc_summary), "./output/R/红细胞频数分布表.xlsx"))   # 频数表数据，便于核对
  try(gtsave(publication_rbc_table, "./output/R/红细胞频数分布表_gt版.docx"))                 # 发表级三线表（Word）
  try(gtsave(publication_rbc_table, "./output/R/红细胞频数分布表_gt版.png"))                  # 发表级三线表（图片）
}


# === 3.4 补充：其他分组方法示例 ===

# 方法1：等宽分组（已在上面演示）
# 优点：组距相等，便于比较；缺点：可能导致某些组频数过少

# 方法2：等频分组（分位数分组）
# 将数据分为频数大致相等的几组
rbc_quantile_groups <- rbc_df %>%
  mutate(
    # 使用分位数创建4个等频组
    rbc_quartile = cut(x,
                       breaks = quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1)),  # 四分位数
                       include.lowest = TRUE,
                       labels = c("Q1 (25%)", "Q2 (25%)", "Q3 (25%)", "Q4 (25%)"))
  ) %>%
  count(rbc_quartile, name = "Frequency") %>%
  mutate(Percentage = round(Frequency / sum(Frequency) * 100, 1))

# 使用gt制作等频分组表
quartile_table <- rbc_quantile_groups %>%
  gt() %>%
  tab_header(
    title = "Table 5. Quartile Distribution of Red Blood Cell Count",
    subtitle = "Equal frequency grouping (approximately 25% each)"
  ) %>%
  cols_label(
    rbc_quartile = "Quartile Group",
    Frequency = "Frequency",
    Percentage = "Percentage (%)"
  ) %>%
  # 为每个四分位数组设置不同的背景色
  tab_style(
    style = cell_fill(color = "#E8F5E8"),                 # 浅绿色
    locations = cells_body(rows = rbc_quartile == "Q1 (25%)")
  ) %>%
  tab_style(
    style = cell_fill(color = "#E8F0FF"),                 # 浅蓝色
    locations = cells_body(rows = rbc_quartile == "Q2 (25%)")
  ) %>%
  tab_style(
    style = cell_fill(color = "#FFF8E8"),                 # 浅黄色
    locations = cells_body(rows = rbc_quartile == "Q3 (25%)")
  ) %>%
  tab_style(
    style = cell_fill(color = "#FFE8E8"),                 # 浅红色
    locations = cells_body(rows = rbc_quartile == "Q4 (25%)")
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_footnote(
    footnote = "Each quartile contains approximately 25% of the observations",
    locations = cells_title(groups = "subtitle")
  )

quartile_table

# === 3.5 频数表的统计描述 ===

# 基于频数表计算描述性统计量
rbc_descriptive <- rbc_df %>%
  summarise(
    `Sample Size` = n(),                                   # 样本量
    `Mean` = round(mean(x), 2),                           # 均值
    `Median` = round(median(x), 2),                       # 中位数
    `Std Dev` = round(sd(x), 2),                          # 标准差
    `Min` = round(min(x), 2),                             # 最小值
    `Max` = round(max(x), 2),                             # 最大值
    `Range` = round(max(x) - min(x), 2),                  # 极差
    `Q1` = round(quantile(x, 0.25), 2),                   # 第一四分位数
    `Q3` = round(quantile(x, 0.75), 2),                   # 第三四分位数
    `IQR` = round(IQR(x), 2)                              # 四分位距
  ) %>%
  # 转换为长格式以便制作表格
  pivot_longer(everything(), names_to = "Statistic", values_to = "Value")

# 制作描述性统计表
descriptive_table <- rbc_descriptive %>%
  gt() %>%
  tab_header(
    title = "Table 6. Descriptive Statistics of Red Blood Cell Count",
    subtitle = "Summary statistics for 138 normal adult women"
  ) %>%
  cols_label(
    Statistic = "Statistical Measure",
    Value = "Value (×10¹²/L)"
  ) %>%
  # 为不同类型的统计量设置不同样式
  tab_style(
    style = list(
      cell_fill(color = "#F0F8FF"),                        # 浅蓝色背景
      cell_text(weight = "bold")
    ),
    locations = cells_body(rows = Statistic %in% c("Sample Size", "Mean", "Median"))
  ) %>%
  tab_style(
    style = cell_fill(color = "#F5F5F5"),                 # 浅灰色背景
    locations = cells_body(rows = Statistic %in% c("Std Dev", "Range", "IQR"))
  ) %>%
  cols_align(align = "center", columns = Value) %>%
  tab_footnote(
    footnote = "All values except Sample Size are in units of ×10¹²/L",
    locations = cells_column_labels(columns = Value)
  )

descriptive_table
