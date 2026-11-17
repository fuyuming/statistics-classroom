# 统计课堂32：RR、OR与log-Wald区间（R / RStudio，UTF-8）
# 数据是人为设定的教学数据，不是真实患者资料。按 0--5 节从上到下运行。
# 改为TRUE可保存；环境变量用于教师批量核对。
SAVE_RESULTS <- FALSE
if (identical(Sys.getenv("RR_OR_SAVE_RESULTS"), "1")) SAVE_RESULTS <- TRUE
project_dir <- getwd()
data_path <- file.path(project_dir, "数据", "rehab_rates.csv")
if (!file.exists(data_path)) stop("请把RStudio项目根目录设为当前工作目录。")

# 0 读取与检查
dat <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
print(dat)
if (!identical(names(dat), c("group", "success", "failure"))) stop("输入字段必须为group, success, failure")
if (nrow(dat) != 2 || !setequal(dat$group, c("new", "control")) || anyDuplicated(dat$group)) stop("group必须恰好有唯一的new和control两行")
counts <- dat[, c("success", "failure")]
if (anyNA(counts) || any(!is.finite(as.matrix(counts)))) stop("success/failure必须是有限且非缺失的数值")
if (any(as.matrix(counts) <= 0) || any(as.matrix(counts) %% 1 != 0)) stop("success/failure必须是非负整数；本练习同时计算RR和OR，请使用四格正整数；零计数需另选方法")
new <- dat[dat$group == "new", ]
control <- dat[dat$group == "control", ]
a <- new$success
b <- new$failure
c <- control$success
d <- control$failure
n1 <- a + b
n0 <- c + d
p1 <- a / n1
p0 <- c / n0

# 1 点估计：RD是概率差，RR是概率比，OR是odds之比
rd <- p1 - p0
rr <- p1 / p0
odds1 <- p1 / (1 - p1)
odds0 <- p0 / (1 - p0)
or_value <- odds1 / odds0
cat(sprintf("RD=%.6f（%.1f个百分点），RR=%.6f，OR=%.6f\n", rd, rd * 100, rr, or_value))
cat(sprintf("odds(new)=%.6f，odds(control)=%.6f\n", odds1, odds0))

# 2 SE：先算样本率方差、log/logit斜率和变换后方差
var_p1 <- p1 * (1 - p1) / n1
var_p0 <- p0 * (1 - p0) / n0
slope_log_p1 <- 1 / p1
slope_log_p0 <- 1 / p0
slope_logit_p1 <- 1 / (p1 * (1 - p1))
slope_logit_p0 <- 1 / (p0 * (1 - p0))
var_log_p1 <- slope_log_p1^2 * var_p1
var_log_p0 <- slope_log_p0^2 * var_p0
var_logit_p1 <- slope_logit_p1^2 * var_p1
var_logit_p0 <- slope_logit_p0^2 * var_p0
se_log_rr <- sqrt(var_log_p1 + var_log_p0)
se_log_or <- sqrt(var_logit_p1 + var_logit_p0)
se_log_rr_cells <- sqrt(1 / a - 1 / n1 + 1 / c - 1 / n0)
se_log_or_cells <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
cat(sprintf("var(p1)=%.8f, var(p0)=%.8f\n", var_p1, var_p0))
cat(sprintf("log斜率=%.6f, %.6f；logit斜率=%.6f, %.6f\n", slope_log_p1, slope_log_p0, slope_logit_p1, slope_logit_p0))
cat(sprintf("变换后方差：log p1=%.8f, log p0=%.8f；logit p1=%.8f, logit p0=%.8f\n", var_log_p1, var_log_p0, var_logit_p1, var_logit_p0))
cat(sprintf("SE_log_RR=%.15f（四格公式=%.15f）；SE_log_OR=%.15f（四格公式=%.15f）\n", se_log_rr, se_log_rr_cells, se_log_or, se_log_or_cells))

# 3 log尺度区间，再exp回原尺度；qnorm(.975)是双侧95%的正侧界值
z <- qnorm(.975)
log_rr <- log(rr)
log_or <- log(or_value)
log_rr_lower <- log_rr - z * se_log_rr
log_rr_upper <- log_rr + z * se_log_rr
log_or_lower <- log_or - z * se_log_or
log_or_upper <- log_or + z * se_log_or
# exp把log尺度端点还原为原比值尺度。
rr_lower <- exp(log_rr_lower)
rr_upper <- exp(log_rr_upper)
or_lower <- exp(log_or_lower)
or_upper <- exp(log_or_upper)
cat(sprintf("z=%.12f；RR原尺度区间=[%.8f, %.8f]；OR=[%.8f, %.8f]\n", z, rr_lower, rr_upper, or_lower, or_upper))

# 4 交换组别：实际重算SE与log区间，再核对倒数端点
swapped_rr <- p0 / p1
swapped_or <- odds0 / odds1
swapped_se_log_rr <- sqrt(1 / c - 1 / n0 + 1 / a - 1 / n1)
swapped_se_log_or <- sqrt(1 / c + 1 / d + 1 / a + 1 / b)
swapped_log_rr <- log(swapped_rr)
swapped_log_or <- log(swapped_or)
swapped_log_rr_lower <- swapped_log_rr - z * swapped_se_log_rr
swapped_log_rr_upper <- swapped_log_rr + z * swapped_se_log_rr
swapped_log_or_lower <- swapped_log_or - z * swapped_se_log_or
swapped_log_or_upper <- swapped_log_or + z * swapped_se_log_or
swapped_rr_lower <- exp(swapped_log_rr_lower)
swapped_rr_upper <- exp(swapped_log_rr_upper)
swapped_or_lower <- exp(swapped_log_or_lower)
swapped_or_upper <- exp(swapped_log_or_upper)
cat(sprintf("交换后RR=%.8f，OR=%.8f；重算SE=%.15f, %.15f\n", swapped_rr, swapped_or, swapped_se_log_rr, swapped_se_log_or))
cat(sprintf("RR区间倒数核对最大误差=%.3e；OR区间倒数核对最大误差=%.3e\n", max(abs(swapped_rr_lower - 1 / rr_upper), abs(swapped_rr_upper - 1 / rr_lower)), max(abs(swapped_or_lower - 1 / or_upper), abs(swapped_or_upper - 1 / or_lower))))

# 5 人数同时乘4：点估计不变，SE按1/sqrt(4)=1/2缩小
scale <- 4
se_rr4 <- sqrt(1 / (scale * a) - 1 / (scale * n1) + 1 / (scale * c) - 1 / (scale * n0))
se_or4 <- sqrt(1 / (scale * a) + 1 / (scale * b) + 1 / (scale * c) + 1 / (scale * d))
cat(sprintf("人数乘%g后SE_log_RR=%.15f，SE_log_OR=%.15f\n", scale, se_rr4, se_or4))

# R额外软件对照：这是软件接口比较，不将Fisher结果当作Wald验证
if (!requireNamespace("DescTools", quietly = TRUE)) stop("请先安装DescTools")
tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
rd_r <- DescTools::BinomDiffCI(a, n1, c, n0, conf.level = .95, method = "wald")
rr_r <- DescTools::RelRisk(tab, conf.level = .95, method = "wald", delta = 0)
or_r <- DescTools::OddsRatio(tab, conf.level = .95, method = "wald")
cat("DescTools delta=0结果：\n")
print(rbind(RD = rd_r, RR = rr_r, OR = or_r))
cat("RelRisk默认delta=.5（与本页delta=0不同）：\n")
print(DescTools::RelRisk(tab, conf.level = .95, method = "wald"))
cat("Fisher条件OR（不是Wald验证）：\n")
print(fisher.test(tab))

# 可选保存：统一metric,value schema
metrics <- c("p1", "p0", "RD", "RR", "OR", "z", "se_log_RR", "se_log_OR", "log_RR_lower", "log_RR_upper", "log_OR_lower", "log_OR_upper", "RR_lower", "RR_upper", "OR_lower", "OR_upper")
values <- c(p1, p0, rd, rr, or_value, z, se_log_rr, se_log_or, log_rr_lower, log_rr_upper, log_or_lower, log_or_upper, rr_lower, rr_upper, or_lower, or_upper)
results <- data.frame(metric = metrics, value = values)
if (SAVE_RESULTS) {
    dir.create(file.path(project_dir, "运行结果", "R"), recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(project_dir, "运行结果", "R", "results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}
