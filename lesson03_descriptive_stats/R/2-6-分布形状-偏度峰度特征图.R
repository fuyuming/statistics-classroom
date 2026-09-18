# 2026-09-18｜第二次课·统计描述：分布形态与偏度（正态／右偏／左偏，均数 vs 中位数）
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-6-分布形状-偏度峰度特征图.R
# 运行：在RStudio打开本文件后Source，或命令行Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("ggplot2", "gridExtra", "dplyr"), repos="https://cloud.r-project.org")
# options(digits=3)设置有效数字显示，不改变原数据；统计解释要保留单位和分母。
# 结果保存目录：output/R（与 output/Python、output/MATLAB 并列）
dir.create("output/R", showWarnings = FALSE, recursive = TRUE)

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
  packages <- c("ggplot2", "gridExtra", "dplyr")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) {
    message("首次运行，正在安装缺失的包：", paste(missing,collapse=", "))
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  failed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(failed)) stop("以下包尚未安装成功，请检查网络及安装日志：",paste(failed,collapse=", "))
})

# ===== 偏度示意图 (Skewness Demonstration) =====

# 加载必要的包
library(ggplot2)   # 用于绘图
library(gridExtra) # 用于拼接多个图形

# 定义绘图函数：展示不同偏度的分布
plot_skew <- function(skew_type = c("right", "normal", "left")){
  # match.arg() 确保输入参数只能是指定的三个值之一
  skew_type <- match.arg(skew_type)
  
  # 创建x轴数据点，从-5到5，共1000个点
  x <- seq(-5, 5, length.out = 1000)
  
  # 根据偏度类型生成不同的概率密度
  if(skew_type == "right"){
    # 右偏分布：使用对数正态分布
    # dlnorm() 生成对数正态分布的密度函数
    # 这里X=Y-5（Y为对数正态变量），因此均数、中位数也必须减5
    y <- dlnorm(x + 5, meanlog = 0, sdlog = 0.5)   
    title <- "右偏"
    mean_pos <- exp(0.5^2/2) - 5    # 均值位置（右偏时均值>中位数）
    median_pos <- exp(0) - 5  # 中位数位置
  } else if(skew_type == "left"){
    # 左偏分布：使用翻转的对数正态分布
    # -x+5 实现水平翻转
    y <- dlnorm(-x + 5, meanlog = 0, sdlog = 0.5)  
    title <- "左偏"
    mean_pos <- 5 - exp(0.5^2/2)   # 均值位置（左偏时均值<中位数）
    median_pos <- 5 - exp(0) # 中位数位置
  } else {
    # 正态分布：完全对称，无偏度
    # dnorm() 生成正态分布的密度函数
    y <- dnorm(x, mean = 0, sd = 1)                
    title <- "不偏"
    mean_pos <- 0      # 均值位置（对称时均值=中位数）
    median_pos <- 0    # 中位数位置
  }
  
  # 将x和y组合成数据框，便于ggplot使用
  df <- data.frame(x = x, y = y)
  
  # 使用ggplot2绘制图形
  ggplot(df, aes(x, y)) +
    # 绘制密度曲线
    geom_line(color = "red", linewidth = 1) +
    # 添加均值的垂直线（紫色虚线）
    geom_vline(xintercept = mean_pos, color = "purple", 
               linetype = "dotted", linewidth = 1) +
    # 添加中位数的垂直线（蓝色虚线）
    geom_vline(xintercept = median_pos, color = "blue", 
               linetype = "dotted", linewidth = 1) +
    # 添加均值标签
    annotate("text", x = mean_pos, y = max(y)*0.8, 
             label = "均值", color = "purple", 
             angle = 90, vjust = -0.5) +
    # 添加中位数标签
    annotate("text", x = median_pos, y = max(y)*0.8, 
             label = "中位数", color = "blue", 
             angle = 90, vjust = -0.5) +
    # 设置图形标题
    ggtitle(title) +
    # 使用简洁主题，字体大小16
    theme_minimal(base_size = 16) +
    # 标题居中
    theme(plot.title = element_text(hjust = 0.5))
}

# 生成三种不同偏度的图形
p1 <- plot_skew("right")   # 右偏图
p2 <- plot_skew("normal")  # 正态图
p3 <- plot_skew("left")    # 左偏图

# 使用grid.arrange将三个图形水平排列
# ncol = 3 表示排成一行三列
grid.arrange(p1, p2, p3, ncol = 3)

# ===== 峰度示意图 (Kurtosis Demonstration) =====

# 加载必要的包
library(ggplot2) # 绘图包
library(dplyr)   # 数据处理包

# 自定义拉普拉斯分布的概率密度函数
# R语言没有内置的拉普拉斯分布函数，所以需要自己定义
dlaplace <- function(x, mu = 0, b = 1){
  # 拉普拉斯分布的PDF公式：f(x) = 1/(2b) * exp(-|x-μ|/b)
  # mu: 位置参数（均值）
  # b: 尺度参数（控制分布的宽度）
  return(1/(2*b) * exp(-abs(x - mu)/b))
}

# 定义x轴范围：从-6到6，共1000个数据点
x <- seq(-6, 6, length.out = 1000)

# 计算三种不同峰度分布的密度值
y_laplace <- dlaplace(x, mu = 0, b = 1)   # 尖峭峰：峰度 > 3
y_norm    <- dnorm(x, mean = 0, sd = 1)   # 正态峰：峰度 = 3
y_unif    <- dunif(x, min = -3, max = 3)  # 平阔峰：峰度 < 3

# 将数据整理成长格式，便于ggplot绘制多条线
df <- data.frame(
  # x值重复3次（对应3种分布）
  x = rep(x, 3),
  # y值：依次是拉普拉斯、正态、均匀分布的密度
  y = c(y_laplace, y_norm, y_unif),
  # 分布类型标签（包含峰度信息）
  type = factor(rep(c("尖峭峰 Laplace (Kurtosis = 6)",      # 高峰度
                      "正态峰 Normal (Kurtosis = 3)",        # 标准峰度
                      "平阔峰 Uniform (Kurtosis = 1.8)"),    # 低峰度
                    each = length(x)))  # 每个标签重复1000次
)

# 绘制峰度对比图
ggplot(df, aes(x = x, y = y, color = type, linetype = type)) +
  # 绘制线条，线宽为1
  geom_line(linewidth = 1) +
  # 设置图形标题和坐标轴标签
  labs(title = "峰度示意图 (Kurtosis Comparison)",
       x = "数值", 
       y = "密度") +
  # 手动设置线条颜色
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  # 使用简洁主题，基础字体大小24
  theme_minimal(base_size = 24) +
  theme(
    # 移除图例标题
    legend.title = element_blank(),
    # 主标题居中
    plot.title = element_text(hjust = 0.5),
    # 图例位置：右上角（0.95, 0.95表示相对位置）
    legend.position = c(0.95, 0.95),              
    # 图例对齐方式：右上对齐
    legend.justification = c("right", "top")       
  )

# ===== 知识点总结 =====
# 
# 偏度 (Skewness)：
# - 右偏（正偏）：均值 > 中位数，长尾向右
# - 左偏（负偏）：均值 < 中位数，长尾向左  
# - 无偏：均值 = 中位数，分布对称
#
# 峰度 (Kurtosis)：
# - 尖峭峰 (Leptokurtic)：峰度 > 3，分布较尖锐
# - 正态峰 (Mesokurtic)：峰度 = 3，标准正态分布
# - 平阔峰 (Platykurtic)：峰度 < 3，分布较平坦
#
# 函数说明：
# - seq()：生成等差数列
# - rep()：重复元素
# - factor()：创建因子变量
# - match.arg()：参数匹配检查
# - annotate()：在图上添加文本注释
# - geom_vline()：添加垂直参考线

# === 结果保存（默认开）：把本脚本最后绘制的一张图存到本文件夹 ===
# 只想到控制台看、不生成图片时，先运行 options(medstats.export = FALSE)。
if (isTRUE(getOption("medstats.export", TRUE))) {
  try(ggplot2::ggsave("./output/R/2-6-分布形状-偏度峰度特征图.png", width = 8, height = 5, dpi = 150))
}
