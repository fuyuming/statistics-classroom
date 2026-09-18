# 2026-09-18｜第二次课·统计描述：变异系数（CV）；本节因课堂时间关系未展开，脚本一并公开。
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 页码：99–100；脚本名：2-10-变异系数.R
# 运行：在RStudio打开本文件后Source，或命令行Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("showtext", "tidyverse", "openxlsx", "DescTools", "patchwork", "rcompanion", "psych"), repos="https://cloud.r-project.org")
# options(digits=3)设置有效数字显示，不改变原数据；统计解释要保留单位和分母。
# 图表导出为可选步骤：需要时先运行 options(medstats.export=TRUE)。
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
  packages <- c("showtext", "tidyverse", "openxlsx", "DescTools", "patchwork", "rcompanion", "psych")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) {
    message("首次运行，正在安装缺失的包：", paste(missing,collapse=", "))
    install.packages(missing, repos="https://cloud.r-project.org")
  }
  failed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if(length(failed)) stop("以下包尚未安装成功，请检查网络及安装日志：",paste(failed,collapse=", "))
})

###### 环境设置 #####
options(digits = 3)
Sys.setenv(LANGUAGE = "en")
options(scipen = 200)

# 防止中文出现乱码
library(showtext)
showtext_auto()

# ============================================================================
# 变异系数分析流程
# 重点：比较不同变量的相对离散程度
# ============================================================================

#####============================================================#####
#####  第一部分：环境准备与数据加载
#####============================================================#####

# 清理工作环境
# 原脚本清空环境语句已停用：请用新的R会话运行，避免删除其他分析对象。

# 加载必需的包
required_packages <- c(
  "tidyverse",     # 数据处理和可视化
  "openxlsx",      # Excel文件读取
  "DescTools",     # 描述性统计（主要使用这个包）
  "patchwork",     # 图形拼接
  "rcompanion",    # 正态性可视化
  "psych"          # 心理统计
)

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#####============================================================#####
#####  第二部分：数据读取与准备
#####============================================================#####

# 读取新生儿数据
# 公开版：新生儿身高体重以 CSV 形式随仓库提供（与课堂用 06.xlsx 数值一致）
data6 <- read.csv("./data/newborn_2025.csv")

# 查看数据基本信息
glimpse(data6)

# 只保留身高（height_cm）与体重（weight_kg）两列
data6 <- data6 %>% 
  dplyr::select(height_cm, weight_kg)

# 提取向量用于后续分析
height_vec <- data6$height_cm    # 提取身高向量
weight_vec <- data6$weight_kg    # 提取体重向量
n_sample <- nrow(data6)          # 样本数量

#####============================================================#####
#####  第三部分：基本统计量分析
#####============================================================#####

# === 3.1 身高的基本统计量 ===
(height_mean <- Mean(height_vec))         # 身高算术均值
(height_median <- Median(height_vec))     # 身高中位数
(height_sd <- sd(height_vec))             # 身高标准差
(height_var <- Var(height_vec))           # 身高方差

# === 3.2 体重的基本统计量 ===
(weight_mean <- Mean(weight_vec))         # 体重算术均值
(weight_median <- Median(weight_vec))     # 体重中位数
(weight_sd <- sd(weight_vec))             # 体重标准差
(weight_var <- Var(weight_vec))           # 体重方差

# === 3.3 分布形态指标 ===
height_skewness <- Skew(height_vec)       # 身高偏度
height_kurtosis <- Kurt(height_vec)       # 身高峰度
weight_skewness <- Skew(weight_vec)       # 体重偏度
weight_kurtosis <- Kurt(weight_vec)       # 体重峰度

# 偏度判断
height_skew_interpretation <- ifelse(abs(height_skewness) < 0.5, "近似对称分布",
                                     ifelse(abs(height_skewness) < 1, "中度偏态", "高度偏态"))
height_skew_interpretation
weight_skew_interpretation <- ifelse(abs(weight_skewness) < 0.5, "近似对称分布",
                                     ifelse(abs(weight_skewness) < 1, "中度偏态", "高度偏态"))
weight_skew_interpretation
#####============================================================#####
#####  第四部分：数据分布可视化
#####============================================================#####

# === 4.1 核密度图对比 ===
# p1: 绘制身高（height_cm）的核密度图，包含数据点的rug plot
p1 <- ggplot(data = data6, mapping = aes(x = height_cm)) +
  geom_rug(alpha = 0.6) +  # 显示数据点的位置
  geom_density(aes(y = after_stat(density)), col = "blue", lwd = 0.8) +  # 绘制核密度曲线
  geom_vline(xintercept = height_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = height_median, color = "green", linetype = "dashed", linewidth = 1) +
  annotate("text", x = height_mean, y = max(density(height_vec)$y) * 0.8, 
           label = paste0("均值 = ", round(height_mean, 2)), 
           color = "red", hjust = -0.1, size = 3) +
  annotate("text", x = height_median, y = max(density(height_vec)$y) * 0.6, 
           label = paste0("中位数 = ", round(height_median, 2)), 
           color = "green", hjust = -0.1, size = 3) +
  labs(title = "身高分布",
       subtitle = paste0("偏度 = ", round(height_skewness, 3)),
       x = "身高 (cm)", y = "密度") +
  theme_minimal()

# p2: 绘制体重（weight_kg）的核密度图
p2 <- ggplot(data = data6, mapping = aes(x = weight_kg)) +
  geom_rug(alpha = 0.6) +  # 显示数据点的位置
  geom_density(aes(y = after_stat(density)), col = "blue", lwd = 0.8) +  # 绘制核密度曲线
  geom_vline(xintercept = weight_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = weight_median, color = "green", linetype = "dashed", linewidth = 1) +
  annotate("text", x = weight_mean, y = max(density(weight_vec)$y) * 0.8, 
           label = paste0("均值 = ", round(weight_mean, 2)), 
           color = "red", hjust = -0.1, size = 3) +
  annotate("text", x = weight_median, y = max(density(weight_vec)$y) * 0.6, 
           label = paste0("中位数 = ", round(weight_median, 2)), 
           color = "green", hjust = -0.1, size = 3) +
  labs(title = "体重分布",
       subtitle = paste0("偏度 = ", round(weight_skewness, 3)),
       x = "体重 (kg)", y = "密度") +
  theme_minimal()

# 合并两个图像，左右展示
p1 | p2

# === 4.2 正态性检验可视化 ===
# 将两个正态性检验图放在一张图里
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# 绘制height_cm的直方图并拟合正态分布
plotNormalHistogram(data6$height_cm, 
                    main = "身高正态性检验", 
                    col = "lightblue")

# 绘制weight_kg的直方图并拟合正态分布
plotNormalHistogram(data6$weight_kg, 
                    main = "体重正态性检验", 
                    col = "lightcoral")

# 重置绘图参数
par(mfrow = c(1, 1))


#####============================================================#####
#####  第五部分：正态性检验
#####============================================================#####

# === 5.1 Shapiro-Wilk正态性检验 ===
height_shapiro_result <- shapiro.test(data6$height_cm)  # 身高正态性检验
weight_shapiro_result <- shapiro.test(data6$weight_kg)  # 体重正态性检验

# 正态性检验结果汇总
normality_results <- data.frame(
  Variable = c("身高", "体重"),
  Statistic = c(
    round(height_shapiro_result$statistic, 4),
    round(weight_shapiro_result$statistic, 4)
  ),
  P_value = c(
    round(height_shapiro_result$p.value, 4),
    round(weight_shapiro_result$p.value, 4)
  ),
  Interpretation = c(
    ifelse(height_shapiro_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据"),
    ifelse(weight_shapiro_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据")
  )
)

normality_results

#####============================================================#####
#####  第六部分：详细统计描述
#####============================================================#####

# === 6.1 使用psych包进行详细统计描述 ===
# 使用psych包中的describe函数生成每个变量的统计汇总
# 包括均值、中位数、标准差等信息
psych_results <- psych::describe(data6, type = 3)
psych_results

#####============================================================#####
#####  第七部分：变异系数计算与比较
#####============================================================#####

# === 7.1 变异系数手工计算 ===
# 变异系数定义: CV = SD/Mean × 100%
height_cv_manual <- height_sd / height_mean
weight_cv_manual <- weight_sd / weight_mean

# === 7.2 使用DescTools包计算变异系数 ===
height_cv_desc <- CoefVar(data6$height_cm)
weight_cv_desc <- CoefVar(data6$weight_kg)

# === 7.3 使用apply函数批量计算 ===
cv_batch <- apply(data6, 2, function(x) { sd(x) / mean(x) })

# 验证计算结果的一致性
cv_calculations <- data.frame(
  变量 = c("身高", "体重"),
  手工计算CV = c(height_cv_manual, weight_cv_manual),
  DescTools_CV = c(height_cv_desc, weight_cv_desc),
  Apply函数CV = cv_batch,
  CV百分比 = c(height_cv_desc * 100, weight_cv_desc * 100)
)

cv_calculations

#####============================================================#####
#####  第八部分：统计量综合比较
#####============================================================#####

# === 8.1 基本统计量比较表 ===
basic_stats_comparison <- data.frame(
  统计量 = c("样本数", "均值", "中位数", "标准差", "方差", "偏度", "峰度"),
  身高 = c(n_sample, 
         round(height_mean, 3), 
         round(height_median, 3), 
         round(height_sd, 3), 
         round(height_var, 3),
         round(height_skewness, 3),
         round(height_kurtosis, 3)),
  体重 = c(n_sample, 
         round(weight_mean, 3), 
         round(weight_median, 3), 
         round(weight_sd, 3), 
         round(weight_var, 3),
         round(weight_skewness, 3),
         round(weight_kurtosis, 3))
)

basic_stats_comparison

# === 8.2 离散程度比较 ===
dispersion_comparison <- data.frame(
  变量 = c("身高", "体重"),
  单位 = c("cm", "kg"),
  均值 = c(round(height_mean, 3), round(weight_mean, 3)),
  标准差 = c(round(height_sd, 3), round(weight_sd, 3)),
  变异系数 = c(round(height_cv_desc, 4), round(weight_cv_desc, 4)),
  变异系数百分比 = c(paste0(round(height_cv_desc * 100, 2), "%"), 
              paste0(round(weight_cv_desc * 100, 2), "%")),
  相对离散程度 = c(
    ifelse(height_cv_desc > weight_cv_desc, "较大", "较小"),
    ifelse(weight_cv_desc > height_cv_desc, "较大", "较小")
  )
)

dispersion_comparison

#####============================================================#####
#####  第九部分：结论与解释
#####============================================================#####

# === 9.1 变异系数比较结论 ===
cv_difference <- abs(height_cv_desc - weight_cv_desc)
higher_cv_var <- ifelse(height_cv_desc > weight_cv_desc, "身高", "体重")
higher_cv_value <- max(height_cv_desc, weight_cv_desc)

conclusion_summary <- data.frame(
  项目 = c("变异系数较大的变量", "最大变异系数", "变异系数差值", "结论"),
  结果 = c(higher_cv_var, 
         round(higher_cv_value, 4),
         round(cv_difference, 4),
         paste0(higher_cv_var, "的相对离散程度更大"))
)

conclusion_summary

# === 9.2 统计学意义解释 ===
cat("=== 变异系数分析结论 ===\n")
cat("身高变异系数:", round(height_cv_desc, 4), "(", round(height_cv_desc * 100, 2), "%)\n")
cat("体重变异系数:", round(weight_cv_desc, 4), "(", round(weight_cv_desc * 100, 2), "%)\n")
cat("结论:", higher_cv_var, "的相对离散程度更大\n")
cat("\n=== 统计学解释 ===\n")
cat("变异系数消除了量纲和数量级的影响，\n")
cat("可以直接比较不同变量的相对离散程度。\n")
cat("变异系数越大，表示数据的相对离散程度越大。\n")

# === 9.3 为什么不能直接比较标准差 ===
cat("\n=== 为什么不能直接比较标准差 ===\n")
cat("身高标准差:", round(height_sd, 3), "cm\n")
cat("体重标准差:", round(weight_sd, 3), "kg\n")
cat("由于单位不同，无法直接比较标准差的大小。\n")
cat("变异系数通过除以均值，消除了单位的影响，\n")
cat("使得不同变量间的离散程度具有可比性。\n")

# CV适用于具有实质零点、均数为正且远离0的变量；本例为身高和体重。
# 改变乘法单位（cm→m）不改变CV，但加常数（摄氏→开尔文）会改变CV。
# 正态性不是CV的计算前提；偏态时需结合分布图解释代表性。
