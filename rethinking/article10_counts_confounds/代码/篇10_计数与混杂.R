# ============================================================
# 精读 10｜计数里的混杂：零膨胀、暴露量与调整
# 对应：McElreath 2023 第 10 讲 Counts and Confounds；教材第 11、12 章
# 用法： Rscript 代码/篇10_计数与混杂.R
#
# 约定：与前几篇同一套工程约定——确定性计算、三语言逐位一致。
# 本篇的主要拟合工具是**显式泊松 IRLS**（迭代重加权最小二乘，固定迭代次数、不提前停止），
# 这样批次系数会被真正估出来，而不是当"已知偏移"塞进模型。
# ============================================================
suppressPackageStartupMessages({ library(tidyverse) })

.args <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .args, value = TRUE)
root <- if (length(.f) > 0) dirname(dirname(normalizePath(sub("^--file=", "", .f[1])))) else getwd()
OUT <- file.path(root, "文章配图"); RES <- file.path(root, "运行结果")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(name, plot, w = 7.4, h = 4.4) {
  ggsave(file.path(OUT, name), plot, width = w, height = h, dpi = 150,
         device = ragg::agg_png, bg = "white")
}
num <- function(x, d = 6) {
  s <- formatC(round(as.numeric(x), d), format = "f", digits = d)
  s <- sub("0+$", "", s); s <- sub("\\.$", "", s)
  ifelse(s == "" | s == "-0", "0", s)
}
write_csv_fixed <- function(df, path) {
  hdr <- paste0('"', names(df), '"', collapse = ",")
  cols <- lapply(df, function(cl) if (is.numeric(cl)) num(cl) else as.character(cl))
  writeLines(c(hdr, do.call(paste, c(cols, sep = ","))), path, useBytes = TRUE)
}
r4 <- function(x) round(as.numeric(x), 4)

# 泊松 IRLS（固定 12 次迭代，不做提前停止，保证三语言一致）
# off＝偏移项（log 暴露量），系数固定为 1，不参与估计
pois_irls <- function(X, y, off = rep(0, nrow(X)), n_iter = 12) {
  b <- rep(0, ncol(X))
  for (it in 1:n_iter) {
    eta <- as.numeric(X %*% b) + off
    mu <- exp(eta)
    W <- mu
    z <- eta + (y - mu) / mu - off        # 工作响应里要去掉偏移
    b <- solve(t(X) %*% (W * X), t(X) %*% (W * z))
  }
  list(beta = b, XtWX = t(X) %*% (W * X))
}

# ---------- 1) 数据：16 个培养单元 ----------
# 处理（0 对照 / 1 接种）、批次（1 = 不同批，与处理**不平衡**：对照组 5 个、处理组 3 个）、
# 接种体积（10 或 20 微升）、是否接种成功（记录里能看到）、每皿菌落数
unit <- data.frame(
  单元 = 1:16,
  处理 = c(0,0,0,0,0, 0,0,0, 1,1,1, 1,1,1,1,1),
  批次 = c(1,1,1,1,1, 0,0,0, 1,1,1, 0,0,0,0,0),
  体积微升 = c(10,20,10,20,10, 20,10,20, 10,20,10, 20,10,20,10,20),
  接种成功 = c(1,1,1,0,1, 1,1,1, 1,1,0, 1,1,1,1,1),
  菌落数 = c(6,13,5,0,7, 6,3,5, 8,19,0, 9,4,8,5,7)
)
write_csv_fixed(unit, file.path(RES, "01_计数数据.csv"))

offs <- log(unit$体积微升 / 10)
y <- unit$菌落数

# 模型 A：只放处理（含体积偏移）—— 不调整批次
XA <- cbind(1, unit$处理)
fitA <- pois_irls(XA, y, offs)
# 模型 B：处理 + 批次 + 偏移 —— 调整批次
XB <- cbind(1, unit$处理, unit$批次)
fitB <- pois_irls(XB, y, offs)

ratio <- function(fit, idx = 2) exp(fit$beta[idx])
se_log <- function(fit, idx = 2) sqrt(solve(fit$XtWX)[idx, idx])
r_tab <- data.frame(
  模型 = c("只放处理（含体积偏移）", "再加批次（混杂）"),
  率比 = r4(c(ratio(fitA), ratio(fitB))),
  下界89 = r4(c(exp(log(ratio(fitA)) - 1.598 * se_log(fitA)), exp(log(ratio(fitB)) - 1.598 * se_log(fitB)))),
  上界89 = r4(c(exp(log(ratio(fitA)) + 1.598 * se_log(fitA)), exp(log(ratio(fitB)) + 1.598 * se_log(fitB)))),
  批次系数 = c("（未估计）", num(r4(fitB$beta[3])))
)
write_csv_fixed(r_tab, file.path(RES, "02_率比调整前后.csv"))
print(r_tab)

# ---------- 2) 零的比例：三种说法 ----------
# ① 观测到的零比例；② 单一泊松在拟合值下预测的零比例；③ 只用"接种成功"的单元拟合后的零比例
muA <- exp(as.numeric(XA %*% fitA$beta) + offs)
p_zero_pois <- mean(exp(-muA))
keep <- unit$接种成功 == 1
fitC <- pois_irls(XB[keep, , drop = FALSE], y[keep], offs[keep])
muC <- exp(as.numeric(XB[keep, , drop = FALSE] %*% fitC$beta) + offs[keep])
zero_tab <- data.frame(
  量 = c("观测到的零比例", "单一泊松预测的零比例（含全部单元）",
         "只看接种成功单元的率比", "只看接种成功单元：预测的零比例"),
  值 = r4(c(mean(y == 0), p_zero_pois, exp(fitC$beta[2]), mean(exp(-muC)))),
  说明 = c("16 个单元里有 2 个计数为 0（都是没接种上的）",
           "泊松把零当成「恰好长了 0 个」，于是严重低估零的出现",
           "把没接种上的单元排除后重估",
           "排除结构性零之后，泊松对「真零」的预测就合理了")
)
write_csv_fixed(zero_tab, file.path(RES, "03_零的比例.csv"))
print(zero_tab)

# ---------- 3) 过离散：泊松标准误是否可信 ----------
muB <- exp(as.numeric(XB %*% fitB$beta) + offs)
pearson <- sum((y - muB)^2 / muB)
df_res <- nrow(unit) - ncol(XB)
phi <- pearson / df_res
se_pois <- se_log(fitB)
se_quasi <- se_pois * sqrt(phi)
od_tab <- data.frame(
  量 = c("Pearson 卡方", "残差自由度", "离散系数 φ（卡方/自由度）",
         "泊松假设下的 log 率比标准误", "按 φ 校正后的标准误", "校正倍数"),
  值 = r4(c(pearson, df_res, phi, se_pois, se_quasi, sqrt(phi)))
)
write_csv_fixed(od_tab, file.path(RES, "04_过离散.csv"))
print(od_tab)

# 过离散校正后的率比区间（用 φ 放大的标准误）
b_adj <- log(ratio(fitB))
rr_quasi <- data.frame(
  量 = c("调整批次后的率比", "泊松假设下的 89% 区间下", "泊松假设下的 89% 区间上",
         "按 φ 校正后的 89% 区间下", "按 φ 校正后的 89% 区间上"),
  值 = r4(c(exp(b_adj),
            exp(b_adj - 1.598 * se_pois), exp(b_adj + 1.598 * se_pois),
            exp(b_adj - 1.598 * se_quasi), exp(b_adj + 1.598 * se_quasi)))
)
write_csv_fixed(rr_quasi, file.path(RES, "06_校正后的率比.csv"))
print(rr_quasi)

write_csv_fixed(data.frame(单元 = 1:16, 处理 = unit$处理, 批次 = unit$批次, 体积微升 = unit$体积微升,
                           接种成功 = unit$接种成功, 菌落数 = y,
                           模型B预测均值 = r4(muB)),
                file.path(RES, "05_拟合值.csv"))

# ---------- 图 ----------
p1 <- ggplot(unit, aes(factor(体积微升), 菌落数, colour = factor(批次))) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.85) +
  labs(title = "先看清三件事：体积、批次、以及两个空皿",
       subtitle = "大体积皿计数更多（要放偏移）；批次 1 的计数也偏高（要考虑是不是混杂）；两个 0 是没接种上的结构性零",
       x = "接种体积（微升）", y = "每皿菌落数", colour = "批次")
save_fig("01_数据总览.png", p1, 7.4, 4.4)

p2 <- ggplot(r_tab, aes(模型, 率比)) +
  geom_col(fill = c("#e8b4a8", "#167d80"), width = 0.5) +
  geom_errorbar(aes(ymin = 下界89, ymax = 上界89), width = 0.12, colour = "#163d48") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#9aa8a8") +
  geom_text(aes(label = sprintf("%.4f", 率比)), vjust = -0.9, size = 4) +
  labs(title = "不调整批次，率比就偏了",
       subtitle = paste0("批次系数 = ", r_tab$批次系数[2], "（明显不为 0）；把它放进模型后，处理率比从 ",
                         r_tab$率比[1], " 变成 ", r_tab$率比[2]),
       x = NULL, y = "率比（处理 vs 对照）")
save_fig("02_率比调整.png", p2, 7.4, 4.4)

p3 <- ggplot(zero_tab[1:2, ], aes(量, 值)) +
  geom_col(fill = c("#c1462c", "#167d80"), width = 0.45) +
  geom_text(aes(label = sprintf("%.4f", 值)), vjust = -0.6, size = 4) +
  labs(title = "泊松解释不了「零太多」",
       subtitle = paste0("观测到 ", r4(mean(y == 0)), " 的零；单一泊松只预测 ", r4(p_zero_pois),
                         "——差距就是「没长出来」那一层"),
       x = NULL, y = "零的比例") +
  ylim(0, max(zero_tab$值[1:2]) * 1.3)
save_fig("03_零太多.png", p3, 7.4, 4.2)

p4 <- ggplot(data.frame(写法 = c("泊松标准误", "按离散系数校正"), 标准误 = c(se_pois, se_quasi)),
             aes(写法, 标准误)) +
  geom_col(fill = c("#e8b4a8", "#c1462c"), width = 0.45) +
  geom_text(aes(label = sprintf("%.4f", 标准误)), vjust = -0.6, size = 4) +
  labs(title = "过离散会让区间假窄",
       subtitle = paste0("本例离散系数 φ = ", r4(phi), "；把标准误乘 √φ = ", r4(sqrt(phi)),
                         " 之后，率比区间才是可信的宽度"),
       x = NULL, y = "log 率比的标准误")
save_fig("04_过离散.png", p4, 7.4, 4.2)

p5 <- ggplot(unit, aes(factor(处理), 菌落数)) +
  geom_boxplot(fill = "#167d80", alpha = 0.45, width = 0.4) +
  geom_jitter(aes(colour = factor(批次)), width = 0.1, size = 2.6) +
  labs(title = "调整之后再看这组差异",
       subtitle = "箱线是原始计数（含体积与批次的影响）；调整后的率比见第 2 张图",
       x = "处理（0 对照 / 1 接种）", y = "每皿菌落数", colour = "批次")
save_fig("05_处理对比.png", p5, 7.4, 4.2)

cat("\n完成：5 张图 + 5 份 CSV 已写入\n")
