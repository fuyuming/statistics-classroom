# 统计课堂40：假设心电干扰案例，不是实测数据。UTF-8。
# 目标：趋势/得分检验，Logistic拟合，OR/Wald区间，三种检验及预测。
# 原题：假设200名不同患者各提供一段心电记录，每段有100个已标注QRS。
# 将记录随机分为4组，每段仅接受一种干扰强度，各患者记录相互独立。
# 固定算法处理全部记录；漏检超过5个QRS判为失败，不超过5个判为成功。
# 4组各50人，50 Hz干扰峰峰值为0/.1/.2/.3 mV，预评分为0:3，失败数为5/10/18/27。
# 问题：是否有有序趋势？趋势统计量怎样对应Logistic模型？
# 每0.1 mV的失败胜算比及双侧95% Wald区间是多少，各组预测失败概率是多少？
# 为什么趋势检验与模型默认P值不同？采用双侧alpha=.05。
# 打开本课Rproj，从01.1逐节运行。使用R自带stats/graphics及lmtest；安装说明见README。

# 01.1 准备与输入 ----
cat("题目：检查共享数据，并准备失败事件、总数和观察失败率。\n")
if (!requireNamespace("lmtest", quietly = TRUE)) {
  stop("请先在R控制台运行 install.packages(\"lmtest\")，然后重新开始本节。")
}
# Rproj默认工作目录为项目根；也支持从“代码”目录运行。
root_dir <- getwd()
if (!file.exists(file.path(root_dir, "数据", "共享教学数据.csv"))) {
  root_dir <- normalizePath(file.path(root_dir, ".."))
}
data_path <- file.path(root_dir, "数据", "共享教学数据.csv")
stopifnot(file.exists(data_path))
result_dir <- file.path(root_dir, "运行结果", "R")
# 默认不保存；也可在运行前设save_outputs <- TRUE，或使用环境变量。
if (!exists("save_outputs")) {
  save_outputs <- identical(Sys.getenv("ARTICLE40_SAVE_OUTPUTS"), "1")
}
d <- read.csv(data_path, stringsAsFactors = FALSE)
print(d)
required <- c("score", "amplitude_mV", "failed", "n")
stopifnot(all(required %in% names(d)))
stopifnot(all(is.finite(unlist(d[required]))), all(d$n > 0),
          all(d$n == floor(d$n)), all(d$failed == floor(d$failed)),
          all(d$failed >= 0), all(d$failed <= d$n))
mV_per_score <- 0.1
stopifnot(length(mV_per_score) == 1, is.finite(mV_per_score), mV_per_score > 0)
# 教学单位参数：改为1即可把score按mV使用；括号必留，默认换算因子为1。
d$score <- d$score * (0.1 / mV_per_score)
d$passed <- d$n - d$failed
d$observed_rate <- d$failed / d$n
N <- sum(d$n)

# 01.2 二项Logistic拟合、coef、OR及95% Wald区间 ----
cat("问题：先看模型结果：失败概率如何随评分变化，斜率的OR及95% Wald区间是多少？\n")
# cbind第一列为事件数failed，第二列为非事件数passed。
# binomial选择二项计数；logit将失败概率连到评分直线。
# epsilon/maxit设迭代精度和最多迭代次数，便于后面复算。
ctrl <- glm.control(epsilon = 1e-12, maxit = 100)
fit <- glm(cbind(failed, passed) ~ score, data = d,
           family = binomial(link = "logit"), control = ctrl)
co <- coef(summary(fit))
# [行名, 列名]读取评分系数；不同于只有截距的共同率模型。
beta <- unname(coef(fit)["score"])
beta_se <- co["score", "Std. Error"]
z_wald <- co["score", "z value"]
P_Wald <- co["score", "Pr(>|z|)"]
OR <- exp(beta)
# qnorm输入累计概率：双侧95%区间两端各留2.5%，取0.975分位。
OR_CI <- exp(beta + c(-1, 1) * qnorm(.975) * beta_se)
print(co)
print(c(OR = OR, low = OR_CI[1], high = OR_CI[2]))
OR_per_0.1mV <- exp(beta * 0.1 / mV_per_score)
cat(sprintf("OR按每增加1评分（%.6g mV）解释；双侧95%% Wald区间=[%.10f, %.10f]；换算为每0.1 mV OR=%.10f。\n",
            mV_per_score, OR_CI[1], OR_CI[2], OR_per_0.1mV))

# 01.3 四组预测与曲线图 ----
cat("问题：用刚才的模型预测四组失败概率，和观察失败率如何对应？\n")
# response返回概率尺度；默认的link则返回logit尺度。
d$fitted_rate <- predict(fit, newdata = d, type = "response")
print(d)
# 横轴改为mV；预测输入除以教学单位，换回原模型使用的评分。
amplitude_grid <- seq(0, .3, length.out = 200)
curve_probability <- predict(fit,
  newdata = data.frame(score = amplitude_grid / mV_per_score), type = "response")
plot(d$amplitude_mV, d$observed_rate, pch = 19, ylim = c(0, .65),
     xlab = "50 Hz interference amplitude (mV, peak-to-peak)",
     ylab = "Failure probability")
lines(amplitude_grid, curve_probability, col = "#1f77b4", lwd = 2)
legend("topleft", c("Observed rate", "Logistic fit"),
       pch = c(19, NA), lty = c(NA, 1), col = c("black", "#1f77b4"), bty = "n")

# 01.4 趋势U/I_eff/Q/P与Wald、Rao、LRT比较 ----
cat("问题：趋势Q如何核对模型检验？这里I_eff作为得分U的方差估计使用。\n")
p0 <- sum(d$failed) / N
score_mean <- sum(d$n * d$score) / N
t_score <- d$score - score_mean
U <- sum(t_score * (d$failed - d$n * p0))
I_eff <- p0 * (1 - p0) * sum(d$n * t_score^2)
Q_manual <- U^2 / I_eff
trend <- prop.trend.test(d$failed, d$n, score = d$score)
Q_trend <- unname(trend$statistic)
P_trend <- trend$p.value
null <- glm(cbind(failed, passed) ~ 1, data = d,
            family = binomial(link = "logit"), control = ctrl)
score_test <- anova(null, fit, test = "Rao")
lrt_test <- anova(null, fit, test = "LRT")
Q_Rao <- score_test$Rao[2]
P_Rao <- score_test$`Pr(>Chi)`[2]
Q_LRT <- lrt_test$Deviance[2]
P_LRT <- lrt_test$`Pr(>Chi)`[2]
logLik_null <- as.numeric(logLik(null))
logLik_fit <- as.numeric(logLik(fit))
Q_from_logLik <- 2 * (logLik_fit - logLik_null)
course_lrt <- lmtest::lrtest(null, fit)
print(data.frame(N, p0, score_mean, U, I_eff, Q_manual, Q_trend, P_trend))
print(course_lrt)
stopifnot(abs(course_lrt$Chisq[2] - Q_LRT) < 1e-10)
cat("LogLik两行之差乘2就是似然比Q；不是Rao得分Q。三者都检验beta=0，有限样本P不必相同。\n")
print(data.frame(method = c("Rao", "Wald", "LRT"),
                 P = c(P_Rao, P_Wald, P_LRT)))

# 01.5 结果与可选保存 ----
cat("报告：以下对象给出样本量、趋势检验、单位OR及拟合结果。\n")
result <- data.frame(N, p0, score_mean, U, I_eff, Q_manual, Q_trend, Q_Rao,
 P_trend, P_Rao, Q_LRT, P_LRT, logLik_null, logLik_fit, Q_from_logLik, beta, beta_se, OR,
 OR_low = OR_CI[1], OR_high = OR_CI[2], Z_Wald = z_wald, P_Wald, mV_per_score)
print(result, digits = 12)
cat(sprintf("可读报告：本次%d段记录；未校正趋势Q=%.10f（df=1），双侧P=%.6g；二项Logistic每%.6g mV的OR=%.10f，双侧95%% Wald区间=[%.10f, %.10f]。\n",
            N, Q_trend, P_trend, mV_per_score, OR, OR_CI[1], OR_CI[2]))
if (save_outputs) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(result_dir, "R_summary.csv"), row.names = FALSE)
  write.csv(d, file.path(result_dir, "R_group_predictions.csv"), row.names = FALSE)
  png(file.path(result_dir, "R_observed_fitted.png"), width = 900, height = 650, res = 120)
  plot(d$amplitude_mV, d$observed_rate, pch = 19, ylim = c(0, .65),
       xlab = "50 Hz interference amplitude (mV, peak-to-peak)",
       ylab = "Failure probability")
  lines(amplitude_grid, curve_probability, col = "#1f77b4", lwd = 2)
  legend("topleft", c("Observed rate", "Logistic fit"),
         pch = c(19, NA), lty = c(NA, 1), col = c("black", "#1f77b4"), bty = "n")
  dev.off()
}
