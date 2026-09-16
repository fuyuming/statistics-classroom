# 统计课堂：本年度班级的专业构成（公开汇总版）
# 改编自首课 R/01_roster.R。公开文件每行是一类专业，不是一名学生。
# 原来的个体序号、专业代码、上课方式均不进入公开包。
# 先双击医学统计学.Rproj，再从头运行本文件。
# 首次准备：install.packages("showtext")
if (!requireNamespace("showtext", quietly = TRUE)) {
  stop('请先运行 install.packages("showtext")')
}
library(showtext)
showtext_auto()
class_font <- "wqy-microhei"
if (!file.exists("data/major_counts.csv")) stop("请先打开配套R项目，确认data文件夹存在。")

# 1. 导入CSV：major为文本，count为整数；不是个体名单。
d <- read.csv("data/major_counts.csv", fileEncoding = "UTF-8-BOM",
              colClasses = c("character", "integer"))
print(d)
str(d)
# 2. 检查缺失、重复类别、负数，再计算；不能把缺失人数自动补成零。
stopifnot(!anyNA(d), !anyDuplicated(d$major), all(d$count >= 0))
# 3. 分母必须是人数之和，nrow(d)只等于专业类别数。
N <- sum(d$count)
if (N == 0) stop("总人数为0，无法计算构成比。")
d$percent <- 100 * d$count / N
stopifnot(N == 18, nrow(d) == 3, d$count[d$major == "生物医学工程"] == 12)
print(d)
# 4. 保存结果；只在展示时四舍五入，计算保留原有精度。
dir.create("output", showWarnings = FALSE)
write.csv(d, "output/01_major_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")
# 5. 绘图函数：上方显示人数最多的专业，标签同时显示人数与构成比。
draw_major <- function(d) {
  d <- d[order(d$count), ]
  par(mar = c(4, 9, 4, 2), family = class_font)
  pos <- barplot(d$count, names.arg = d$major, horiz = TRUE, las = 1,
                 xlim = c(0, max(d$count) * 1.4), col = "#267D86", border = NA,
                 xlab = "人数", main = "本年度课堂的专业构成")
  text(d$count + 0.25, pos, labels = sprintf("%d人（%.2f%%）", d$count, d$percent),
       adj = 0, family = class_font)
  mtext(sprintf("本班合计%d人；公开数据仅保留专业人数汇总", sum(d$count)),
        side = 3, line = 0.4, cex = 0.85)
}
# 在RStudio显示图，再导出PNG；字体与分辨率统一。
draw_major(d)
local({
  old_dpi <- showtext_opts()$dpi
  showtext_opts(dpi = 160)
  png("output/01_major_composition.png", width = 1600, height = 900, res = 160)
  on.exit({dev.off(); showtext_opts(dpi = old_dpi)}, add = TRUE)
  draw_major(d)
})
# 解释：这是本班构成，不能自动推广到全校或其他年份。
