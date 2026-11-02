# 统计课堂22：Tukey检验与学生化极差q
# 数据：第12课例9-1，三个小麦品种各10个独立地块；原题未给单位。
# RStudio打开“统计课堂22.Rproj”，再Source本文件。需要stats与car包。
# 结果写入运行结果/R；先运行完整文件，再修改后续练习。

# 1. 查看数据 ----
wide <- read.csv("数据/小麦产量.csv", fileEncoding = "UTF-8-BOM")
print(wide)
stopifnot(ncol(wide) == 3, nrow(wide) == 10, !anyNA(wide))
k <- ncol(wide)
dat <- stack(wide)
names(dat) <- c("yield", "variety")
print(head(dat))
# 同行记录不构成配对；三个品种是独立组。
out <- "运行结果/R"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) {
  write.csv(x, file.path(out, name), row.names = FALSE, fileEncoding = "UTF-8")
}

# 2. 先做条件诊断：正态性看组内形状，方差检验看组间离散程度。
# P>.05只是“未发现明显偏离”的线索，不能证明前提，也不机械决定是否做ANOVA。
diagnostics <- data.frame(group = 1:k, shapiro_P = sapply(wide, function(x) shapiro.test(x)$p.value))
save_table(diagnostics, "diagnostics.csv")
print(diagnostics)
png(file.path(out, "qq_diagnostics.png"), width = 1500, height = 550, res = 130)
par(mfrow = c(1, 3))
for (i in 1:ncol(wide)) {
  qqnorm(wide[[i]], main = paste("Variety", i))
  qqline(wide[[i]], col = "steelblue")
}
dev.off()
if (!requireNamespace("car", quietly = TRUE)) stop("需要安装 car 包：install.packages('car')")
# Levene比较观测值到组中心的绝对偏差；中位数中心较稳健。
levene_median <- car::leveneTest(yield ~ variety, data = dat, center = median)
levene_mean <- car::leveneTest(yield ~ variety, data = dat, center = mean)
# Bartlett使用样本方差，对偏离正态较敏感。
bartlett <- bartlett.test(yield ~ variety, data = dat)
variance_tests <- data.frame(
  method = c("Levene median", "Levene mean", "Bartlett"),
  statistic = c(levene_median[1, "F value"], levene_mean[1, "F value"], unname(bartlett$statistic)),
  P = c(levene_median[1, "Pr(>F)"], levene_mean[1, "Pr(>F)"], bartlett$p.value)
)
save_table(variance_tests, "variance_tests.csv")
print(variance_tests)

# Levene比较各组到中心的绝对偏差（这里分别用中位数、均值中心）；
# Bartlett检验方差相等但依赖正态性，三者都只是辅助诊断，不机械选检验。
# 3. ANOVA及Tukey主分析（诊断之后） ----
fit <- aov(yield ~ variety, data = dat)
print(summary(fit))
tukey <- TukeyHSD(fit, conf.level = 0.95)$variety
print(tukey)

# 4. MSE、自由度与均值标准误从哪里来？ ----
k <- ncol(wide)
n <- nrow(wide)
N <- n * k
nu <- N - k
means <- colMeans(wide)
# 每组减去自己的均值；每列平方和除以n-1得到该组样本方差。
within_ss <- numeric(k)
for (i in seq_len(k)) {
  residual_i <- wide[[i]] - means[i]
  within_ss[i] <- sum(residual_i^2)
}
sse <- sum(within_ss)
mse <- sse / nu
ssa <- n * sum((means - mean(as.matrix(wide)))^2)
f_value <- (ssa / (k - 1)) / mse
p_anova <- pf(f_value, k - 1, nu, lower.tail = FALSE)
scale_q <- sqrt(mse / n)
group_result <- data.frame(group = 1:k, n = n, mean = means,
                           variance = within_ss / (n - 1), within_ss = within_ss)
anova_result <- data.frame(SSA = ssa, SSE = sse, df1 = k - 1, df2 = nu,
                           MSE = mse, F = f_value, P = p_anova)
save_table(group_result, "groups.csv")
save_table(anova_result, "anova.csv")
stopifnot(abs(mse - deviance(fit) / df.residual(fit)) < 1e-10)


# 5. 三对差值：q、调整P、同时区间 ----
# 行顺序和TukeyHSD相同：2-1、3-1、3-2。
pair_index <- matrix(c(2, 1, 3, 1, 3, 2), ncol = 2, byrow = TRUE)
difference <- numeric(3)
for (r in 1:3) {
  difference[r] <- means[pair_index[r, 1]] - means[pair_index[r, 2]]
}
q_value <- abs(difference) / scale_q
# lower.tail=FALSE计算右尾；nmeans是整个比较族的组数。
p_adj <- ptukey(q_value, nmeans = k, df = nu, lower.tail = FALSE)
q_critical <- qtukey(0.95, nmeans = k, df = nu)
half_width <- q_critical * scale_q
pairs <- data.frame(pair = c("2-1", "3-1", "3-2"), difference = difference,
                    q = q_value, p_adj = p_adj,
                    lower = difference - half_width, upper = difference + half_width)
print(pairs)
save_table(pairs, "pairs.csv")
# 独立的手工分步计算，与标准函数的四列核对。
stopifnot(max(abs(pairs$difference - tukey[, "diff"])) < 1e-10)
stopifnot(max(abs(pairs$p_adj - tukey[, "p adj"])) < 1e-8)
stopifnot(max(abs(pairs$lower - tukey[, "lwr"])) < 1e-7)
stopifnot(max(abs(pairs$upper - tukey[, "upr"])) < 1e-7)
cat("品种2-1均值差-90，同时95%区间约[-137.59,-42.41]。\n")
cat("品种3-1的区间跨0；同时报告差值、区间和调整P。\n")

# 6. k=2时，q与双侧t的关系 ----
q_grid <- c(0, 1, 3, 6)
p_from_q <- ptukey(q_grid, nmeans = 2, df = nu, lower.tail = FALSE)
# q=sqrt(2)*abs(t)，所以t的右尾位置为q/sqrt(2)，双侧乘2。
p_from_t <- 2 * pt(q_grid / sqrt(2), df = nu, lower.tail = FALSE)
save_table(data.frame(q = q_grid, q_tail = p_from_q, t_two_sided = p_from_t), "k2_check.csv")
stopifnot(max(abs(p_from_q - p_from_t)) < 1e-7)

# 把单对门槛用于三对，会有多大整组误报概率？
# 最大绝对均值差就是极差；至少一对越界等价于最大差越界。
# FWER在总体均值全相等及模型成立时由q分布计算，每次抽样重新估计MSE。
q_limits <- c(sqrt(2) * qt(0.975, df = nu), q_critical)
family_control <- data.frame(
  method = c("Unadjusted pairs", "Tukey family"),
  qcrit = q_limits, diff_threshold = q_limits * scale_q,
  fwer_global_null = ptukey(q_limits, nmeans = k, df = nu, lower.tail = FALSE)
)
save_table(family_control, "family_control.csv")
print(family_control)
# 本例品种3与2相差72：用相同MSE比较单对P与Tukey调整P。
print(c(single_P = 2 * pt(-72 / (sqrt(2) * scale_q), df = nu),
        tukey_P = ptukey(72 / scale_q, nmeans = k, df = nu, lower.tail = FALSE)))

# 7. 模拟q分布：固定nu=27，比较组数2、3、5 ----
# 标准正态组均值误差Z独立于共同方差估计U/nu；U服从卡方(nu)。
# 正态模型下，Q=(max Z-min Z)/sqrt(U/nu)。
# 使用同一行前2、3、5个Z，以及同一个U，直观看到增加组数的影响。
set.seed(20261101)
B <- 100000
z <- matrix(rnorm(B * 5), nrow = B, ncol = 5)
u_chisq <- rchisq(B, df = nu)
denominator <- sqrt(u_chisq / nu)
k_values <- c(2, 3, 5)
simulation <- data.frame(k = k_values, df = nu, q95 = NA_real_,
                         simulated_q95 = NA_real_, exceedance = NA_real_)
q_samples <- matrix(NA_real_, nrow = B, ncol = 3)
for (j in seq_along(k_values)) {
  selected <- z[, seq_len(k_values[j]), drop = FALSE]
  sample_range <- apply(selected, 1, max) - apply(selected, 1, min)
  q_samples[, j] <- sample_range / denominator
  simulation$q95[j] <- qtukey(0.95, nmeans = k_values[j], df = nu)
  simulation$simulated_q95[j] <- quantile(q_samples[, j], 0.95)
  simulation$exceedance[j] <- mean(q_samples[, j] > simulation$q95[j])
}
colnames(q_samples) <- c("k2", "k3", "k5")
save_table(as.data.frame(q_samples), "q_samples.csv")
save_table(simulation, "simulation.csv")
print(simulation)
# 100000次独立模拟时，比例在.05附近的蒙特卡洛标准误约.000689。
mc_se <- sqrt(0.05 * 0.95 / B)
stopifnot(all(abs(simulation$exceedance - 0.05) < 5 * mc_se))

# 8. 改单位练习：所有产量乘10 ----
scaled_fit <- aov(I(yield * 10) ~ variety, data = dat)
scaled_tukey <- TukeyHSD(scaled_fit)$variety
unit_check <- data.frame(pair = pairs$pair, difference_x10 = scaled_tukey[, "diff"],
                         lower_x10 = scaled_tukey[, "lwr"], upper_x10 = scaled_tukey[, "upr"],
                         p_x10 = scaled_tukey[, "p adj"])
save_table(unit_check, "unit_check.csv")
stopifnot(max(abs(unit_check$difference_x10 - pairs$difference * 10)) < 1e-7)
stopifnot(max(abs(unit_check$lower_x10 - pairs$lower * 10)) < 1e-6)
stopifnot(max(abs(unit_check$p_x10 - pairs$p_adj)) < 1e-8)
cat("改单位后：差值和区间乘10，q和P保持不变。\n")

# 9. q分布图（Q-Q诊断已在ANOVA前完成） ----
png(file.path(out, "q_tail.png"), width = 950, height = 700, res = 130)
x <- seq(0, 7, length.out = 141)
plot(x, ptukey(x, 2, nu, lower.tail = FALSE), type = "l", col = 1,
     xlab = "q", ylab = "Pr(Q > q)", ylim = c(0, 1))
for (j in 2:3) lines(x, ptukey(x, k_values[j], nu, lower.tail = FALSE), col = j)
abline(h = 0.05, lty = 2)
legend("topright", legend = paste("k =", k_values), col = 1:3, lty = 1)
dev.off()
writeLines(c(R.version.string, "stats随R安装；所有核对通过。"), file.path(out, "version.txt"))

# 10. 全局H0下的整组误报教学模拟（每轮三组同时判断） ----
set.seed(20261102)
sim_trials <- 10000
sim_n <- 10
sim_mean <- 720
sim_variance <- 1842
sim_df <- k * (sim_n - 1)
sim_values <- array(rnorm(sim_trials * k * sim_n, mean = sim_mean,
                          sd = sqrt(sim_variance)), dim = c(sim_n, k, sim_trials))
sim_means <- matrix(0, nrow = sim_trials, ncol = k)
sim_mse <- numeric(sim_trials)
sim_max_q <- numeric(sim_trials)
sim_t <- matrix(0, sim_trials, 3)
sim_p <- matrix(0, sim_trials, 3)
for (b in seq_len(sim_trials)) {
  draw <- sim_values[, , b]
  sim_means[b, ] <- colMeans(draw)
  sim_sse <- sum((sweep(draw, 2, sim_means[b, ], "-"))^2)
  sim_mse[b] <- sim_sse / sim_df
  difference_sim <- c(sim_means[b,2]-sim_means[b,1],
                      sim_means[b,3]-sim_means[b,1],
                      sim_means[b,3]-sim_means[b,2])
  sim_t[b,] <- difference_sim / sqrt(2*sim_mse[b]/sim_n)
  sim_p[b,] <- 2*pt(-abs(sim_t[b,]), df=sim_df)
  sim_max_q[b] <- max(abs(outer(sim_means[b, ], sim_means[b, ], "-"))) /
    sqrt(sim_mse[b] / sim_n)
}
sim_thresholds <- c("Unadjusted pairs" = sqrt(2) * qt(0.975, sim_df),
                    "Tukey family" = qtukey(0.95, k, sim_df))
sim_false <- vapply(sim_thresholds, function(threshold) sum(sim_max_q > threshold), integer(1))
null_family_simulation <- data.frame(
  method = names(sim_thresholds), trials = sim_trials,
  false_positive_trials = as.integer(sim_false), rate = as.numeric(sim_false / sim_trials)
)
# 各对算双侧P，再按每轮“有没有至少一对P<.05”计数。
sim_t_any <- rowSums(sim_p < .05) > 0
sim_tukey_any <- sim_max_q > sim_thresholds[2]
stopifnot(sum(sim_t_any) == sim_false[1], sum(sim_tukey_any) == sim_false[2])
save_table(data.frame(trial=seq_len(sim_trials),
                     t_any=sim_t_any, tukey_any=sim_tukey_any), "null_trial_flags.csv")
example_trial <- which(sim_t_any)[1]
example_draw <- sim_values[,,example_trial]
colnames(example_draw) <- c("品种1","品种2","品种3")
write.csv(example_draw, "数据/模拟示例.csv", row.names=FALSE, fileEncoding="UTF-8")
example_result <- data.frame(trial=example_trial, pair=c("2-1","3-1","3-2"),
  t=sim_t[example_trial,], P=sim_p[example_trial,],
  reject_t=sim_p[example_trial,] < .05, MSE=sim_mse[example_trial])
save_table(example_result,"null_example.csv")
print(example_result)
save_table(null_family_simulation, "null_family_simulation.csv")
print(null_family_simulation)

# SSA演示：保留三组均值，改变各组人数后，重新算加权总均值及平方和。
demo_means <- c(756, 666, 738)
demo_ns <- list(c(10, 10, 10), c(5, 10, 15))
ssa_demo <- do.call(rbind, lapply(seq_along(demo_ns), function(i) {
  sizes <- demo_ns[[i]]
  grand <- sum(sizes * demo_means) / sum(sizes)
  data.frame(design = c("equal_n_10_10_10", "unequal_n_5_10_15")[i],
    grand_mean = grand, ssa = sum(sizes * (demo_means - grand)^2))
}))
save_table(ssa_demo, "ssa_demo.csv")
stopifnot(all(ssa_demo$grand_mean == c(720,717)), all(ssa_demo$ssa == c(45360,40230)))

# 11. 图05采用本脚本的R模拟结果；坐标与PNG同步保存。
false_n <- null_family_simulation$false_positive_trials
true_n <- sim_trials-false_n
plot_coords <- data.frame(method=rep(null_family_simulation$method,each=2),
                         status=rep(c("有误报","无误报"),2),
                         count=as.vector(rbind(false_n,true_n)))
write.csv(plot_coords,"文章配图/05-整组误报模拟.csv",row.names=FALSE,fileEncoding="UTF-8")
png("文章配图/05-整组误报模拟.png",width=1560,height=1536,res=240,
    pointsize=16,family="PingFang SC",type="cairo")
par(mar=c(.3,.3,.3,.3))
plot.new()
plot.window(xlim=c(-600,10600),ylim=c(0,11),xaxs="i",yaxs="i")
text(0,10.5,"10000轮试验，有多少轮误报？",adj=0,cex=1.18,font=2)
text(0,9.85,"R教学模拟：三组总体均值均为720",adj=0,cex=.9)
text(0,9.35,"共同方差1842，每组10个独立正态观测",adj=0,cex=.9)
for(j in 1:2){
  yy <- c(7.6,4.7)[j]
  text(0,yy+.75,c("逐对t检验","Tukey检验")[j],adj=0,cex=1,font=2)
  text(10000,yy+.75,sprintf("误报%d轮（%.2f%%）",false_n[j],100*false_n[j]/sim_trials),
       adj=1,cex=.9,col="#b34b46")
  rect(0,yy-.35,false_n[j],yy+.35,col="#d95f59",border=NA)
  rect(false_n[j],yy-.35,10000,yy+.35,col="#a8d5ba",border=NA)
  text(false_n[j]+true_n[j]/2,yy,sprintf("%d轮未误报",true_n[j]),cex=.95,col="#244b3a")
}
segments(0,3.7,10000,3.7)
for(xx in c(0,5000,10000)){
  segments(xx,3.7,xx,3.58)
  text(xx,3.28,format(xx,scientific=FALSE),cex=.9)
}
text(5000,2.65,"试验轮数（每轮同时比较三对）",cex=.9)
rect(0,1.65,450,1.95,col="#d95f59",border=NA)
text(650,1.8,"至少一对误报",adj=0,cex=.9)
rect(0,.9,450,1.2,col="#a8d5ba",border=NA)
text(650,1.05,"三对均未误报",adj=0,cex=.9)
dev.off()
article_fig <- file.path("..","文章配图")
if(dir.exists(article_fig)) file.copy("文章配图/05-整组误报模拟.png",article_fig,overwrite=TRUE)
