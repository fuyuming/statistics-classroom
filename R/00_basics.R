# R基础练习：每一段可以独立观察输出；数据是教学示例。
# 1. 类型：numeric、integer、character、logical。
weight <- 2.3
count <- 7L
label <- "control"
is_valid <- TRUE
print(typeof(weight)); print(typeof(count)); print(typeof(label)); print(typeof(is_valid))
# 2. 向量同类元素；list可混合；data.frame按列组织；factor表示类别。
x <- c(1.8, 2.0, NA, 2.4)
print(x[1]); print(x[x > 2 & !is.na(x)])
group <- factor(c("control", "inoculated"))
mat <- matrix(1:6, nrow = 2)
meta <- list(name = "demo", values = x)
tab <- data.frame(id = c("C01", "C02"), weight = c(1.8, 2.0))
print(group); print(mat); print(meta$values); print(tab$weight)
# 3. 条件判断：if需要一个非缺失的TRUE/FALSE。
valid_x <- x[!is.na(x)]
if (length(valid_x) >= 2) {
  print(sd(valid_x))
} else {
  message("有效观测不足两个")
}
# 4. 循环：依次取元素；for内部操作写在花括号中。
for (value in valid_x) {
  print(value * 1000) # 克转成毫克
}
# 5. 函数：输入x，返回去除缺失后的均值；全缺失时返回NA。
mean_valid <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}
print(mean_valid(c(1.8, 2.0, NA))) # 1.9
stopifnot(isTRUE(all.equal(mean_valid(c(1.8, 2.0, NA)), 1.9)))
