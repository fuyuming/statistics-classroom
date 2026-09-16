# 统计课堂：大数定律的抛硬币演示（由首课脚本整理）
# 先双击根目录的医学统计学.Rproj，再从头运行本文件。
# 首次准备：install.packages(c("showtext", "ggplot2"))
# 模拟不涉及任何学生记录。固定种子只用于本环境复现。
if (!requireNamespace("showtext", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop('请先运行 install.packages(c("showtext", "ggplot2"))')
}
if (!file.exists("医学统计学.Rproj")) stop("请先打开配套R项目，再运行脚本。")
library(showtext)
library(ggplot2)
showtext_auto()
class_font <- "wqy-microhei"

# ============================================================================
# 大数定律可视化演示 - 抛硬币实验（简化版）
# ============================================================================
# 目的：通过抛硬币实验演示大数定律，展示样本均值如何收敛到理论概率
# 理论基础：独立同分布的重复试验下，频率向概率收敛（大数定律）。
# 概率是模型中的抽象量，频率是实际比例；有限次模拟不是定理的证明。
# ============================================================================



# ============================================================================
# 第一步：实验设计和参数设置
# ============================================================================

# 定义抛硬币实验的总次数
num_toss <- 1000  # 进行1000次抛硬币实验

# 设置随机种子以确保结果的可重复性
set.seed(2023)

print(paste("开始进行", num_toss, "次抛硬币实验"))
print("理论概率：正面朝上的概率 = 0.5")

# ============================================================================
# 第二步：生成随机实验数据
# ============================================================================

# 随机生成num_toss个抛硬币结果
# 0表示反面朝上（花朝上），1表示正面朝上（国徽朝上）
toss <- sample(0:1, size = num_toss, replace = TRUE)

# 创建实验次数的序列，用于x轴坐标
iteration <- 1:num_toss

print(paste("实际正面次数：", sum(toss)))
print(paste("实际正面频率：", round(mean(toss), 4)))

# ============================================================================
# 第三步：数据处理和累积统计
# ============================================================================

# 创建数据框，整合所有实验信息
df <- data.frame(
  iteration = iteration,  # 实验次数（1到1000）
  toss = toss            # 每次抛硬币的结果（0或1）
)

# 计算累积均值：展示大数定律的收敛过程
# cumsum(df$toss)：计算到第i次为止正面出现的总次数
# df$iteration：到第i次为止的总实验次数
# 累积均值 = 累积正面次数 / 累积实验次数
df$cum_mean <- cumsum(df$toss) / df$iteration

# ============================================================================
# 第四步：创建可视化 - 累积均值收敛图
# ============================================================================

# 绘制折线图，显示累积均值随实验次数的变化（大数定律演示）
convergence_plot <- ggplot(df, aes(x = iteration, y = cum_mean)) +

  # 绘制累积均值的折线
  geom_line(color = "darkblue",      # 深蓝色线条
            linewidth = 1.2) +            # 线条粗细

  # 添加理论概率的水平参考线
  geom_hline(yintercept = 0.5,       # y = 0.5的水平线
             color = "red",          # 红色
             linetype = "dashed",    # 虚线样式
             linewidth = 1.2) +           # 线条粗细

  # 添加理论概率线的文字标注
  annotate("text",
           x = num_toss * 0.75,      # 标注位置x坐标
           y = 0.62,                 # 标注位置y坐标
           label = "理论概率 = 0.5",
           color = "red",
           size = 4.5,
           fontface = "bold", family = class_font) +

  # 设置坐标轴范围和标签
  ylim(0, 1) +                      # y轴范围0到1
  ylab("累积均值（正面频率）") +
  xlab("实验次数") +

  # 添加标题和说明
  ggtitle("大数定律演示：累积均值收敛到理论概率",
          subtitle = paste("随着实验次数从1增加到", num_toss,
                           "；频率有波动，并非逐次接近0.5")) +

  # 应用简洁主题
  theme_minimal(base_family = class_font) +

  # 自定义主题元素
  theme(
    # 文字大小设置
    axis.text = element_text(size = 12),        # 坐标轴刻度标签
    axis.title = element_text(size = 14, face = "bold"),  # 坐标轴标题
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),    # 主标题
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"), # 副标题

    # 网格线设置
    panel.grid.minor = element_blank(),         # 移除次要网格线
    panel.grid.major = element_line(color = "gray70"),  # 淡化主网格线

    # 背景设置
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),

    # 边距设置
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(convergence_plot)

# ============================================================================
# 第五步：统计分析和结果总结
# ============================================================================

# 计算关键统计量
final_mean <- tail(df$cum_mean, 1)                    # 最终累积均值
theoretical_prob <- 0.5                               # 理论概率
deviation <- abs(final_mean - theoretical_prob)       # 与理论值的绝对偏差
relative_error <- deviation / theoretical_prob * 100  # 相对误差百分比
heads_count <- sum(df$toss)                           # 正面总次数
tails_count <- num_toss - heads_count                 # 反面总次数

print("\n=== 实验结果统计分析 ===")
print(paste("总实验次数：", num_toss))
print(paste("正面次数：", heads_count, "次"))
print(paste("反面次数：", tails_count, "次"))
print(paste("最终累积均值：", round(final_mean, 4)))
print(paste("理论概率：", theoretical_prob))
print(paste("绝对偏差：", round(deviation, 4)))
print(paste("相对误差：", round(relative_error, 2), "%"))

# ============================================================================
# 第六步：大数定律的核心概念解释
# ============================================================================

print("\n=== 大数定律核心概念 ===")
print("大数定律表明：当独立重复试验的次数足够大时，")
print("事件发生的频率将稳定在其概率附近。")
print("")
print("在本实验中：")
print("- 每次抛硬币都是独立的随机试验")
print("- 正面朝上的理论概率为0.5")
print(paste("- 经过", num_toss, "次试验后，实际频率为", round(final_mean, 4)))
print(paste("- 与理论概率的相对偏差为", round(relative_error, 2), "%"))
print("")
print("这条模拟路径演示大数定律；有限次模拟不是证明，误差也不会逐次单调减小。")


# 保存可核对的模拟记录与图片；不需要另外提供原始数据文件。
dir.create("output", showWarnings=FALSE)
write.csv(df, "output/02_coin_trials_R.csv", row.names=FALSE)
local({
  old_dpi <- showtext_opts()$dpi
  on.exit(showtext_opts(dpi=old_dpi))
  showtext_opts(dpi=150)
  ggsave("output/02_coin_R.png", convergence_plot, width=10, height=6, dpi=150)
})
# 练习：改变set.seed()或num_toss，比较路径。次数增加不保证误差逐步减小。
