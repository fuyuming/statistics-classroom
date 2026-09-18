# 2026-09-18｜第二次课·统计描述：几何均数（GM 与 GSD）；本节因课堂时间关系未展开，脚本一并公开。
# 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
# 脚本名：2-9-几何均数.R
# 运行：在RStudio打开本文件后Source，或命令行Rscript 本文件。
# 首次使用需要联网安装缺包；只安装缺失项，不反复升级已有包。
# 可手动执行 install.packages(c("showtext", "tidyverse", "readr", "car", "DescTools", "nortest", "bestNormalize", "rcompanion", "moments", "EnvStats"), repos="https://cloud.r-project.org")
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
  packages <- c("showtext", "tidyverse", "readr", "car", "DescTools", "nortest", "bestNormalize", "rcompanion", "moments", "EnvStats")
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
# 偏态分布数据的描述性统计分析流程
# 重点：对数变换 → 几何均数与几何标准差的计算与应用
# ============================================================================

#####============================================================#####
#####  第一部分：环境准备与数据加载
#####============================================================#####

# 清理工作环境
# 原脚本清空环境语句已停用：请用新的R会话运行，避免删除其他分析对象。

# 加载必需的包
required_packages <- c(
  "tidyverse",     # 数据处理和可视化
  "readr",         # 数据读取
  "car",           # QQ图和统计检验
  "DescTools",     # 描述性统计（主要使用这个包）
  "nortest"        # 正态性检验
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

# 读取浓度数据
data5 <- read.table("./data/05.txt", header = TRUE)

# 查看数据基本信息 - 使用glimpse更简洁
glimpse(data5)

# 提取浓度向量用于后续分析
concentration_vec <- data5$Concentration #提取数据框中的一列为向量(变量)
stopifnot(is.numeric(concentration_vec), all(is.finite(concentration_vec)), all(concentration_vec > 0))
# 对数和几何均数要求正值；零/检出限不能随意加常数处理。
n_sample <- length(concentration_vec) #样本数量

#####============================================================#####
#####  第三部分：原始数据的分布形态分析
#####============================================================#####

# === 3.1 基本统计量 ===
(mean_val <- Mean(concentration_vec))       # 算术均值
median_val <- Median(concentration_vec)    # 中位数
median_val
mode_val <- Mode(concentration_vec)        # 众数
mode_val
# === 3.2 分布形态指标 ===
skewness_val <- Skew(concentration_vec)    # 偏度：|偏度|<0.5对称，<1中度偏态，>=1高度偏态
skewness_val
kurtosis_val <- Kurt(concentration_vec)    # 峰度：接近0正常峰态，>0尖峰态，<0平峰态
kurtosis_val
# 偏度判断
skew_interpretation <- ifelse(abs(skewness_val) < 0.5, "近似对称分布",
                             ifelse(abs(skewness_val) < 1, "中度偏态", "高度偏态"))
skew_interpretation
skew_direction <- ifelse(skewness_val > 0, "右偏", "左偏")
skew_direction

# === 3.3 原始数据可视化 ===

# 柱状图显示算术均值的问题
bar_plot <- ggplot(data = data5, aes(x = Concentration)) +
  geom_bar(fill = "red", alpha = 0.7) +
  geom_vline(xintercept = mean_val, color = "purple", linewidth = 1.2) +
  annotate("text", x = mean_val + 0.15, y = 3,
           label = paste0("算术平均数 ", round(mean_val, 3), "\n能否表示集中趋势？"),
           size = 8, color = "blue") +
  annotate("segment", x = mean_val, y = 2.5, xend = mean_val + 0.05, yend = 3,
           arrow = arrow(type = "closed", length = unit(0.1, "inches")),
           color = "blue") +
  labs(title = "原始数据分布 - 算术均值的代表性问题",
       subtitle = paste0("偏度 = ", round(skewness_val, 3), " (", skew_direction, "偏分布)"),
       x = "浓度", y = "频数") +
  theme_minimal()

bar_plot

# 直方图 + 核密度图
hist_density_plot <- ggplot(data = data5, mapping = aes(x = Concentration)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, 
                 fill = "pink", color = "black", alpha = 0.7) +
  geom_rug(alpha = 0.6) +
  geom_density(color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = mean_val, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median_val, color = "green", linetype = "dashed", linewidth = 1) +
  annotate("text", x = mean_val, y = max(density(concentration_vec)$y) * 0.8, 
           label = paste0("均值 = ", round(mean_val, 3)), 
           color = "red", hjust = -0.1, size = 8) +
  annotate("text", x = median_val, y = max(density(concentration_vec)$y) * 0.6, 
           label = paste0("中位数 = ", round(median_val, 3)), 
           color = "green", hjust = -0.1, size = 8) +
  labs(title = "原始数据分布形态分析",
       subtitle = "直方图 + 核密度曲线",
       x = "浓度", y = "密度") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"))

hist_density_plot

# QQ图检查正态性
car::qqPlot(concentration_vec,
       main = "Q-Q Plot for Normality Assessment", 
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles",
       col = "blue", pch = 16, cex = 1.2, id = TRUE)

qqnorm(concentration_vec,
       main = "Q-Q Plot for Normality Assessment",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles",
       col = "blue", pch = 16, cex = 1.2)
qqline(concentration_vec, col = "red", lwd = 2)


# === 3.4 正态性检验 ===
# Shapiro-Wilk检验（样本量限制：n≤5000）
shapiro_result <- shapiro.test(concentration_vec)

# Lilliefors校正：均数和SD由本样本估计，不能套普通单样本KS检验P值。
# 本例含重复/舍入值，连续模型近似也有局限；检验留作后续学习。
ks_result <- nortest::lillie.test(concentration_vec)

# Anderson-Darling检验（更敏感，来自nortest包）
ad_result <- ad.test(concentration_vec)

# Jarque-Bera检验（基于偏度和峰度的综合检验）
jarque_result <- JarqueBeraTest(concentration_vec)

# 正态性检验结果汇总
normality_results <- data.frame(
  Test = c("Shapiro-Wilk", "Lilliefors (estimated parameters)", "Anderson-Darling", "Jarque-Bera"),
  Statistic = c(
    round(shapiro_result$statistic, 4),
    round(ks_result$statistic, 4),
    round(ad_result$statistic, 4),
    round(jarque_result$statistic, 4)
  ),
  P_value = c(
    round(shapiro_result$p.value, 4),
    round(ks_result$p.value, 4),
    round(ad_result$p.value, 4),
    round(jarque_result$p.value, 4)
  ),
  Interpretation = c(
    ifelse(shapiro_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据"),
    ifelse(ks_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据"),
    ifelse(ad_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据"),
    ifelse(jarque_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据")
  )
)

normality_results

# 分布形态总结
distribution_summary <- data.frame(
  指标 = c("偏度", "峰度", "偏度解释", "峰度解释"),
  数值 = c(round(skewness_val, 4), round(kurtosis_val, 4), 
         paste0(skew_interpretation, "(", skew_direction, ")"),
         ifelse(abs(kurtosis_val) < 0.5, "正常峰态",
                ifelse(kurtosis_val > 0, "尖峰态", "平峰态")))
)

distribution_summary



#####============================================================#####
#####  第四部分：对数变换与变换后的分析
#####============================================================#####

# === 4.1 对数变换 ===
# 对右偏数据进行对数变换，目的：矫正右偏分布，使数据趋于正态
data5 <- data5 %>% 
  mutate(log_concentration = log(Concentration))

log_concentration_vec <- data5$log_concentration

# === 4.2 变换后的分布形态分析 ===
(log_mean <- Mean(log_concentration_vec))      # 对数均值
(log_median <- Median(log_concentration_vec))  # 对数中位数
(log_skewness <- Skew(log_concentration_vec))  # 对数偏度
(log_kurtosis <- Kurt(log_concentration_vec)) # 对数峰度

# === 4.3 变换后的可视化 ===
log_hist_plot <- ggplot(data = data5, mapping = aes(x = log_concentration)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, 
                 fill = "lightblue", color = "black", alpha = 0.7) +
  geom_rug(alpha = 0.6) +
  geom_density(color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = log_mean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = log_median, color = "green", linetype = "dashed", linewidth = 1) +
  annotate("text", x = log_mean, y = max(density(log_concentration_vec)$y) * 0.8, 
           label = paste0("对数均值 = ", round(log_mean, 3)), 
           color = "red", hjust = -0.1) +
  annotate("text", x = log_median, y = max(density(log_concentration_vec)$y) * 0.6, 
           label = paste0("对数中位数 = ", round(log_median, 3)), 
           color = "green", hjust = -0.1) +
  labs(title = "对数变换后的分布形态",
       subtitle = paste0("变换后偏度 = ", round(log_skewness, 3)),
       x = "ln(浓度)", y = "密度") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"))

log_hist_plot

# 变换后的QQ图
#car包
qqPlot(log_concentration_vec,
       main = "Q-Q Plot After Log Transformation", 
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles (log scale)",
       col = "darkgreen", pch = 16, cex = 1.2, id = TRUE)
#基础包
qqnorm(log_concentration_vec, main="对数变换 Q-Q图"); qqline(log_concentration_vec,col="red")


# === 4.4 变换后的正态性检验 ===
log_shapiro_result <- shapiro.test(log_concentration_vec)
log_shapiro_result
ShapiroFranciaTest(log_concentration_vec)

#####============================================================#####
#####  4.5 常见非正态数据数据变换方法总结
#####============================================================#####

library(DescTools)
library(car)
library(bestNormalize)
library(rcompanion)

# (1) 对数变换
data5$Concentration_log <- log(data5$Concentration)

# (2) 平方根变换
data5$Concentration_sqrt <- sqrt(data5$Concentration)

# (3) Box–Cox 最优变换
lambda_bc <- BoxCoxLambda(data5$Concentration) # 自动找最佳 Box–Cox λ
lambda_bc 
data5$Concentration_bc <- BoxCox(data5$Concentration, lambda = lambda_bc)

# (4) Yeo–Johnson 变换 - 使用car包统一实现
# 先找最优lambda
powerTransform_result <- powerTransform(data5$Concentration, family="yjPower")
lambda_yj <- powerTransform_result$lambda
lambda_yj  # 查看lambda值
data5$Concentration_yj <- yjPower(data5$Concentration, lambda = lambda_yj)

# 变换参数总结
# Box-Cox λ = lambda_bc
# Yeo-Johnson λ = lambda_yj

#####============================================================#####
#####  变换效果可视化对比
#####============================================================#####
#install.packages('rcompanion')
library(rcompanion)
# 使用 plotNormalHistogram 进行变换效果对比
par(mfrow=c(2,3), mar=c(4,4,3,1))

# 原始数据
plotNormalHistogram(data5$Concentration, 
                    main="原始数据", 
                    col="lightblue")
# 对数变换
plotNormalHistogram(data5$Concentration_log, 
                    main="对数变换", 
                    col="lightcoral")
# 平方根变换
plotNormalHistogram(data5$Concentration_sqrt, 
                    main="平方根变换", 
                    col="lightpink")
# Box-Cox变换
plotNormalHistogram(data5$Concentration_bc, 
                    main=paste0("Box-Cox (λ=", round(lambda_bc,2), ")"), 
                    col="lightgreen")
# Yeo-Johnson变换
plotNormalHistogram(data5$Concentration_yj, 
                    main=paste0("Yeo-Johnson (λ=", round(lambda_yj,2), ")"), 
                    col="orange")

# 重置绘图参数
par(mfrow=c(1,1))

# Q-Q图对比
par(mfrow=c(2,3), mar=c(4,4,3,1))
qqnorm(data5$Concentration, main="原始数据 Q-Q图"); qqline(data5$Concentration, col='red')
qqnorm(data5$Concentration_log, main="对数变换 Q-Q图"); qqline(data5$Concentration_log, col='red')
qqnorm(data5$Concentration_sqrt, main="平方根变换 Q-Q图"); qqline(data5$Concentration_sqrt, col='red')
qqnorm(data5$Concentration_bc, main="Box-Cox Q-Q图"); qqline(data5$Concentration_bc, col='red')
qqnorm(data5$Concentration_yj, main="Yeo-Johnson Q-Q图"); qqline(data5$Concentration_yj, col='red')

# 重置绘图参数
par(mfrow=c(1,1))

#####============================================================#####
#####  正态性检验对比
#####============================================================#####

# Shapiro-Wilk 正态性检验
transformations <- c("原始数据", "对数变换", "平方根变换", "Box-Cox", "Yeo-Johnson")
variables <- list(data5$Concentration, data5$Concentration_log, data5$Concentration_sqrt,
                  data5$Concentration_bc, data5$Concentration_yj)

# 创建结果表格
results <- data.frame(
  变换方法 = transformations,
  W统计量 = numeric(5),
  p值 = numeric(5),
  正态性 = character(5),
  stringsAsFactors = FALSE
)

for(i in 1:5) {
  test_result <- shapiro.test(variables[[i]])
  results$W统计量[i] <- round(test_result$statistic, 4)
  results$p值[i] <- round(test_result$p.value, 4)
  results$正态性[i] <- ifelse(test_result$p.value > 0.05, "未拒绝正态模型；不等于证明正态", "与正态模型不一致的证据")
}

print(results)



#####============================================================#####
#####  变换方法选择建议
#####============================================================#####

# 变换方法选择原则:
# 1. 对数变换: 适用于右偏数据，要求所有值 > 0
# 2. 平方根变换: 适用于轻度右偏数据，要求所有值 ≥ 0
# 3. Box-Cox变换: 自动寻找最优变换，要求所有值 > 0
# 4. Yeo-Johnson变换: 可处理负值，适用范围最广

# 后续方法备查：仅查看本组候选变换的检验输出，禁止按最大P值自动选择变换。
# best_method <- which.max(results$p值)  # 不执行按P值择优
# 基于Shapiro-Wilk检验的最佳变换方法及其p值
# results$变换方法[best_method]
# results$p值[best_method]




#####============================================================#####
#####  第五部分：几何均数与几何标准差计算
#####============================================================#####

# === 5.1 几何均数计算 ===
# 几何均数定义: G = (X₁ × X₂ × ... × Xₙ)^(1/n) = exp(Σln(Xᵢ)/n)
geometric_mean_manual <- exp(sum(log(concentration_vec)) / length(concentration_vec))
geometric_mean_desc <- Gmean(concentration_vec)  # 使用DescTools包

# === 5.2 几何标准差计算 ===
# 几何标准差定义: GSD = exp(SD(ln(X)))
geometric_sd_desc <- Gsd(concentration_vec)
log_sd <- sd(log_concentration_vec)
geometric_sd_manual <- exp(log_sd)

# 验证计算结果的一致性
geometric_calculations <- data.frame(
  方法 = c("手工计算几何均数", "DescTools几何均数", "手工计算几何标准差", "DescTools几何标准差"),
  结果 = c(geometric_mean_manual, geometric_mean_desc, geometric_sd_manual, geometric_sd_desc)
)

geometric_calculations

#####============================================================#####
#####  第六部分：统计量比较与选择建议
#####============================================================#####

# === 6.1 不同统计量的比较 ===
statistics_comparison <- data.frame(
  统计量 = c("算术均数", "中位数", "几何均数"),
  数值 = c(round(mean_val, 3), round(median_val, 3), round(geometric_mean_desc, 3)),
  适用条件 = c("加法尺度的平均水平", "有序或定量资料的位置", "正值、倍数尺度的中心"),
  本例适用性 = c("可计算；受高值影响较大", "描述典型位置", "描述倍数中心"),
  解释 = c("受极值影响大", "稳健统计量", "反映倍数关系")
)

statistics_comparison



# 保留检验输出供后续课使用，本次可跳讲。
# 样本偏度方向与P值大小是不同问题，不能用P值代替分布图。
moments::agostino.test(concentration_vec)
# 几何标准差为无量纲的倍数因子。
geometric.sd <- function(x) exp(sd(log(x)))
stopifnot(isTRUE(all.equal(geometric.sd(concentration_vec), EnvStats::geoSD(concentration_vec))))
print(c(GM = exp(mean(log(concentration_vec))), GSD = geometric.sd(concentration_vec)))
# GM/GSD到GM*GSD在近似对数正态模型下概括约一个log标准差范围，不是均值置信区间。
