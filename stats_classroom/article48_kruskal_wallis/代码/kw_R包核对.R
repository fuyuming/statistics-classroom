# 统计课堂48：用现成R函数核对渐近、精确与Monte Carlo推断
# 从含本项目README和Rproj的根目录，按第1—4节顺序运行。
# 这是可选入口；原kw_practice.R继续用于学习H、并列校正与手写置换。

# 1 准备：包、数据与设置 ----
needed <- c("coin", "kSamples", "rstatix")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("请按README安装可选包后再运行：", paste(missing, collapse = ", "))
}
root <- getwd()
if (!file.exists(file.path(root, "数据", "rabbit.csv"))) {
  stop("请先打开根目录Rproj，或将工作目录设为含数据文件夹的项目根目录。")
}
save_results <- Sys.getenv("SAVE_RESULTS", "0") == "1"
mc_draws <- 199999L
seed <- 20261130L
rabbit <- read.csv(file.path(root, "数据", "rabbit.csv"))
rabbit$drug <- factor(rabbit$drug) # 分组变量明确为因子，避免同名函数混淆。
bottles <- read.csv(file.path(root, "数据", "six_bottles.csv"))
bottles$condition <- factor(bottles$condition)
cat("每行一只兔：四独立组，每组6只；每兔药物眼相对盐水眼改善指标，单位分。\n")
cat("六瓶教学资料：A={1,2}，B={3,4}，C={5,6} mg/L，每组2瓶，互相独立。\n")
cat("所有R包函数均写包名前缀；coin与rstatix有同名kruskal_test。\n")

# 2 基础R与rstatix：二者均采用卡方近似 ----
base_fit <- stats::kruskal.test(score ~ drug, data = rabbit)
# rstatix把基础检验结果整理成数据框；它不添加精确/MC切换。
tidy_fit <- rstatix::kruskal_test(rabbit, score ~ drug)
print(base_fit)
print(tidy_fit)
stopifnot(abs(unname(base_fit$statistic) - tidy_fit$statistic) < 1e-10,
          abs(base_fit$p.value - tidy_fit$p) < 1e-12)

# 3 kSamples：完整枚举六瓶资料的90种身份划分 ----
# test="KW"选择秩评分；Nsim在exact分支也要覆盖全部合法划分数。
# dist=TRUE保留完整分布，可与手写枚举的9个支持点及频数比较。
exact_fit <- kSamples::qn.test(
  concentration_mg_L ~ condition, data = bottles,
  test = "KW", method = "exact", Nsim = 90, dist = TRUE
)
print(exact_fit)
cat("实际使用方法：", exact_fit$method, "；枚举数：", exact_fit$Nsim, "\n")
cat("qn依次为：H、卡方近似P、精确P\n")
print(exact_fit$qn, digits = 12)
stopifnot(exact_fit$method == "exact", exact_fit$Nsim == 90,
          length(exact_fit$null.dist) == 90,
          abs(exact_fit$qn[1] - 32 / 7) < 1e-10,
          abs(exact_fit$qn[3] - 1 / 15) < 1e-12)

# 4 两个包的Monte Carlo：仍是同一兔例条件零分布 ----
# approximate明确模拟次数；statistic/pvalue分别读取H及模拟P。
set.seed(seed)
coin_fit <- coin::kruskal_test(
  score ~ drug, data = rabbit,
  distribution = coin::approximate(nresample = mc_draws)
)
print(coin_fit)
print(coin::pvalue(coin_fit)) # 本次附99%的二项Monte Carlo置信区间。
coin_H <- as.numeric(coin::statistic(coin_fit))
coin_P <- as.numeric(coin::pvalue(coin_fit)) # 去掉附带属性，再计算MCSE。

# simulated用均匀随机划分抽样；返回method/Nsim须读取确认。
set.seed(seed)
ks_fit <- kSamples::qn.test(score ~ drug, data = rabbit,
                           test = "KW", method = "simulated", Nsim = mc_draws)
print(ks_fit)
stopifnot(ks_fit$method == "simulated", ks_fit$Nsim == mc_draws)
ks_P <- as.numeric(ks_fit$qn[3])
stopifnot(abs(coin_H - unname(base_fit$statistic)) < 1e-9,
          abs(ks_fit$qn[1] - coin_H) < 1e-9)

# 当前记录版本中，这两个包均返回b/B，不自动使用(b+1)/(B+1)。
# 读取包的原始结果，同时列出同一批抽样按加1规则换算的数值。
P <- c(coin_P, ks_P)
b <- round(P * mc_draws)
stopifnot(all(abs(P * mc_draws - b) < 1e-7))
comparison <- data.frame(
  package = c("coin", "kSamples"), H = c(coin_H, as.numeric(ks_fit$qn[1])),
  simulations = mc_draws, exceedances = b, package_P = P,
  MCSE = sqrt(P * (1 - P) / mc_draws),
  P_plus1_same_draws = (b + 1) / (mc_draws + 1)
)
print(comparison, digits = 12)
cat("包版本：\n")
for (p in needed) cat(p, as.character(utils::packageVersion(p)), "\n")
cat("随机序列依具体函数与版本；比较模拟P需结合MC误差，不要求逐值相同。\n")
if (save_results) {
  dest <- file.path(root, "运行结果", "R包核对")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  write.csv(comparison, file.path(dest, "rabbit_package_MC.csv"), row.names = FALSE)
}
