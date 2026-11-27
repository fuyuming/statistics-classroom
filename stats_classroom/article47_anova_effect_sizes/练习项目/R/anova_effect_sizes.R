# 1 读题与共享数据 ----
# RStudio先打开项目根统计课堂47.Rproj；从本节顺序运行，或Source整文件。
# 四种孔结构各5个独立试件，每件一次有效压缩模量MPa；独立、正态、共同方差模型。
# 20数字来自李春喜等《生物统计学》6版例7-1，只改情境单位，非真实支架实验。
# 本题求η²/ω²/偏η²、R²及A1−A3的全部6对Tukey同时95%区间。
# effectsize为必需依赖，安装及路径见README。默认看图；SAVE_RESULTS=1时保存CSV/PNG。
root <- getwd()
stopifnot(file.exists(file.path(root, "data", "scaffold.csv")))
save_results <- identical(Sys.getenv("SAVE_RESULTS"), "1")
out <- file.path(root, "运行结果", "R")
# read.csv读取相对路径CSV，返回数据框：group为孔结构，modulus为MPa。
# 下方将group转为分类因子；模量缺失会立即报错，不自动删试件。
scaffold <- read.csv(file.path(root, "data", "scaffold.csv"))
scaffold$group <- factor(scaffold$group)
print(names(scaffold))
print(head(scaffold))
stopifnot(!anyNA(scaffold$modulus))
# split把模量按组拆开；sd使用n−1分母的样本标准差。
groups <- split(scaffold$modulus, scaffold$group)
summary_table <- data.frame(group = names(groups),
  n = sapply(groups, length), mean = sapply(groups, mean), sd = sapply(groups, sd))
print(summary_table)

# 2 手算：SS、MSE、F及效应量 ----
ns <- sapply(groups, length)
means <- sapply(groups, mean)
N <- nrow(scaffold)
k <- length(groups)
grand <- mean(scaffold$modulus)
SSB <- sum(ns * (means - grand)^2)
SSW <- 0
for (x in groups) {
  SSW <- SSW + sum((x - mean(x))^2)
}
SST <- SSB + SSW
MSE <- SSW / (N - k)
f_stat <- (SSB / (k - 1)) / MSE
# pf返回累计F概率；lower.tail=FALSE取右尾，不额外乘2。
p_value <- pf(f_stat, k - 1, N - k, lower.tail = FALSE)
eta2 <- SSB / SST
omega2_raw <- (SSB - (k - 1) * MSE) / (SST + MSE)
omega2_report <- max(0, omega2_raw)
partial_eta2 <- SSB / (SSB + SSW)
main <- data.frame(N, k, SSB, SSW, SST, MSE, F = f_stat, P = p_value,
  eta2, omega2_raw, omega2_report, partial_eta2)
print(main, digits = 12)
cat(sprintf("样本组间平方和占比%.2f%%；ω²是总体比例的修正点估计。\n",100*eta2))
cat("整体均数检验P<0.05：", p_value < .05, "\n")

# 3 模型与effectsize：核对相同统计对象 ----
# aov默认带截距；modulus~group为四组分类模型，返回拟合对象fit。
fit <- aov(modulus ~ group, data = scaffold)
print(anova(fit))
lm_fit <- lm(modulus ~ group, data = scaffold)
# summary返回列表，r.squared与adj.r.squared是其中两个分量。
main$R2 <- summary(lm_fit)$r.squared
main$adjusted_R2 <- summary(lm_fit)$adj.r.squared
adjusted_manual <- 1 - (SSW/(N-k)) / (SST/(N-1))
stopifnot(isTRUE(all.equal(main$R2, eta2)))
stopifnot(isTRUE(all.equal(main$adjusted_R2, adjusted_manual)))
# FALSE指定普通效应量；ci=NULL只求点估计。返回data.frame中的Eta2、Omega2列。
eta_obj <- effectsize::eta_squared(fit, partial = FALSE, ci = NULL)
partial_obj <- effectsize::eta_squared(fit, partial = TRUE, ci = NULL)
omega_obj <- effectsize::omega_squared(fit, partial = FALSE, ci = NULL)
# 本普通单因素模型partial=TRUE会提示相等，并仍返回Eta2列。
stopifnot(isTRUE(all.equal(eta_obj$Eta2, eta2)))
stopifnot(isTRUE(all.equal(partial_obj$Eta2, partial_eta2)))
stopifnot(isTRUE(all.equal(omega_obj$Omega2, omega2_report)))
print(c(eta2_package = eta_obj$Eta2, omega2_package = omega_obj$Omega2), digits = 12)
print(effectsize::interpret_eta_squared(eta2, rules = "field2013"))
cat("R²与调整R²：", main$R2, main$adjusted_R2, "\n")

# 4 A1−A3：Tukey同时95%区间（全部6对组成比较族） ----
# TukeyHSD输入aov对象，返回列表；group分量每行含diff/lwr/upr/p adj。
tukey <- TukeyHSD(fit, "group", conf.level = .95)$group
# 按因子顺序原始行A3-A1方向相反。取负时上下界需交换：[-upr,-lwr]。
reverse_row <- tukey["A3-A1", ]
comparison <- data.frame(difference = -unname(reverse_row["diff"]),
  lower95 = -unname(reverse_row["upr"]), upper95 = -unname(reverse_row["lwr"]),
  p_adj = unname(reverse_row["p adj"]))
stopifnot(isTRUE(all.equal(comparison$difference, unname(means["A1"]-means["A3"]))))
cat("A1−A3（MPa），全部6对的同时95%区间：\n")
print(comparison, digits = 12)

# 5 辅助完整小题：2×2分母，不展开二因素检验 ----
# 孔结构D1/D2×温度low/high，各格2个独立试件，每件测一次模量MPa。
# 全部人工数D1low[9,11] D2low[11,13] D1high[13,15] D2high[15,17]。
# 独立、正态、共同方差；拟合含交互模型，仅求孔结构η²与偏η²。
fac <- read.csv(file.path(root, "data", "factorial_demo.csv"))
fac$design <- factor(fac$design)
fac$temperature <- factor(fac$temperature)
print(fac)
fac_fit <- aov(modulus ~ design * temperature, data = fac)
# *包含两个主效应及交互作用；本平衡设计Type I平方和正交分解。
fac_tab <- anova(fac_fit)
ssa <- fac_tab["design", "Sum Sq"]
ssb <- fac_tab["temperature", "Sum Sq"]
ssab <- fac_tab["design:temperature", "Sum Sq"]
sse <- fac_tab["Residuals", "Sum Sq"]
fac_sst <- sum((fac$modulus - mean(fac$modulus))^2)
factorial <- data.frame(SS_A=ssa, SS_B=ssb, SS_AB=ssab, SSE=sse, SST=fac_sst,
  A_eta2=ssa/fac_sst, A_partial_eta2=ssa/(ssa+sse))
stopifnot(isTRUE(all.equal(ssa+ssb+ssab+sse, fac_sst)))
print(factorial, digits = 12)

# 6 看图与改单位：全数MPa×1000为kPa，重新拟合 ----
# 横向轻微错开仅避免遮挡；短线为组均数。RStudio图窗显示，不默认另存。
if (save_results) {
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  png(file.path(out,"scaffold_groups.png"), width=1248, height=1080, res=240)
}
plot(rep(seq_along(groups), ns), scaffold$modulus, type="n", xaxt="n",
  xlab="Pore design", ylab="Compressive modulus (MPa)", xlim=c(.6,k+.4))
axis(1, at=seq_along(groups), labels=names(groups))
for (index in seq_along(groups)) {
  x <- groups[[index]]
  points(index + seq(-.075,.075,length.out=length(x)), x, pch=19, col="#3366cc")
  segments(index-.23, mean(x), index+.23, mean(x), lwd=2)
}
if (save_results) {
  dev.off()
}
exercise <- scaffold
exercise$modulus <- exercise$modulus * 1000
ex_fit <- aov(modulus ~ group, data=exercise)
ex_tab <- anova(ex_fit)
ex_b <- ex_tab["group", "Sum Sq"]
ex_w <- ex_tab["Residuals", "Sum Sq"]
ex_t <- ex_b + ex_w
ex_mse <- ex_w / (N-k)
ex_lm <- summary(lm(modulus ~ group, data=exercise))
exercise_result <- data.frame(eta2=ex_b/ex_t,
  omega2_raw=(ex_b-(k-1)*ex_mse)/(ex_t+ex_mse),
  R2=ex_lm$r.squared, adjusted_R2=ex_lm$adj.r.squared)
ex_reverse <- TukeyHSD(ex_fit, "group", conf.level=.95)$group["A3-A1",]
ex_comparison <- data.frame(difference=-unname(ex_reverse["diff"]),
  lower95=-unname(ex_reverse["upr"]), upper95=-unname(ex_reverse["lwr"]),
  p_adj=unname(ex_reverse["p adj"]))
for (name in names(exercise_result)) {
  stopifnot(isTRUE(all.equal(exercise_result[[name]], main[[name]])))
}
for (name in c("difference", "lower95", "upper95")) {
  stopifnot(isTRUE(all.equal(ex_comparison[[name]], 1000*comparison[[name]])))
}
print(exercise_result, digits=12)
print(ex_comparison, digits=12)
cat("单位改变使差值和区间乘1000，无单位效应量与R²保持不变。\n")
if (save_results) {
  write.csv(main, file.path(out,"main_results.csv"), row.names=FALSE)
  write.csv(comparison, file.path(out,"tukey_A1_A3.csv"), row.names=FALSE)
  write.csv(factorial, file.path(out,"factorial_results.csv"), row.names=FALSE)
  write.csv(exercise_result, file.path(out,"exercise_results.csv"), row.names=FALSE)
  write.csv(ex_comparison, file.path(out,"exercise_tukey_A1_A3.csv"), row.names=FALSE)
}
