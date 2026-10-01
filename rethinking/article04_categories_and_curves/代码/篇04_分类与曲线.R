# ============================================================
# 精读 04｜同一批数据，换个问题就该换个模型：分类变量、中心化、曲线
# 对应：McElreath 2023 第 04 讲 Categories and Curves；教材第 4、5 章
# 用法： Rscript 代码/篇04_分类与曲线.R
#
# 约定：
#   · 主线案例＝我们自己的场景（三个处理组各 5 个重复的终点厚度；另一批体系的饱和曲线）
#   · 对照案例＝教材 Howell1（性别→体重；年龄→体重的曲线）
#   · 全部确定性：解析后验（σ 固定为合并组内标准差）、显式最小二乘（指数/样条基自己搭）
#     → 三语言输出逐位一致
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(rethinking)
})

# 目录自动定位：脚本放在 <篇目录>/代码/ 下
.args <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .args, value = TRUE)
root <- if (length(.f) > 0) dirname(dirname(normalizePath(sub("^--file=", "", .f[1])))) else getwd()
OUT  <- file.path(root, "文章配图")
RES  <- file.path(root, "运行结果")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(name, plot, w = 7.4, h = 4.2) {
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
# 显式最小二乘：解正规方程（三语言同一套运算顺序，保证逐位一致）
ols_ls <- function(X, y) {
  XtX <- t(X) %*% X; Xty <- t(X) %*% y
  b <- solve(XtX, Xty)
  res <- y - X %*% b
  list(coef = as.numeric(b), sigma_hat = sqrt(sum(res^2) / (length(y) - ncol(X))))
}

# ---------- 0) 教材对照数据 ----------
data(Howell1)
d <- Howell1[Howell1$age >= 18, c("height", "weight", "male")]
d <- d[order(d$height), ]
write_csv_fixed(data.frame(height = d$height, weight = d$weight, male = d$male),
                file.path(RES, "00_教材数据_Howell1成年人.csv"))

# 两个年龄端点的曲线数据（教材第 4 章末的经典图：体重随年龄非线性）
k <- Howell1
write_csv_fixed(data.frame(age = k$age, weight = k$weight), file.path(RES, "00_教材数据_Howell1全部.csv"))

# ---------- 1) 我们的案例（一）：三个处理组的终点厚度 ----------
groups <- tibble(
  处理 = c(rep("对照", 5), rep("低剂量", 5), rep("高剂量", 5)),
  厚度 = c(41.2, 43.0, 39.8, 42.5, 40.9,
           46.5, 48.1, 45.2, 47.4, 46.0,
           52.8, 54.6, 51.9, 53.7, 52.2)
)
groups$处理 <- factor(groups$处理, levels = c("对照", "低剂量", "高剂量"))
write_csv_fixed(data.frame(处理 = as.character(groups$处理), 厚度 = groups$厚度),
                file.path(RES, "01_我们的案例_三组终点厚度.csv"))

# 合并组内标准差（σ 固定，用作解析后验的抽样标准差）
g_mean <- tapply(groups$厚度, groups$处理, mean)
g_n    <- tapply(groups$厚度, groups$处理, length)
ss_within <- sum(sapply(levels(groups$处理), function(g) {
  v <- groups$厚度[groups$处理 == g]; sum((v - mean(v))^2)
}))
df_within <- nrow(groups) - length(levels(groups$处理))
sigma_pool <- sqrt(ss_within / df_within)

# 索引编码：每个组的均值后验 = N(样本均值, sigma_pool/sqrt(n))（σ 固定）
post_sd <- sigma_pool / sqrt(as.numeric(g_n))
z89 <- qnorm(0.945)
means_df <- tibble(
  组 = levels(groups$处理),
  样本均值 = round(as.numeric(g_mean), 4),
  后验标准差 = round(post_sd, 4),
  区间89下 = round(as.numeric(g_mean) - z89 * post_sd, 4),
  区间89上 = round(as.numeric(g_mean) + z89 * post_sd, 4)
)
write_csv_fixed(means_df, file.path(RES, "02_组均值与后验区间.csv"))
print(means_df)

# 对比：三组两两之差（解析后验；P(>0) 用正态尾部）
pairs_df <- tibble(
  对比 = c("低剂量 - 对照", "高剂量 - 对照", "高剂量 - 低剂量"),
  A = c("低剂量", "高剂量", "高剂量"), B = c("对照", "对照", "低剂量")
) %>%
  mutate(均值差 = round(as.numeric(g_mean[A]) - as.numeric(g_mean[B]), 4),
         差的标准误 = round(sqrt(post_sd[match(A, levels(groups$处理))]^2 +
                                   post_sd[match(B, levels(groups$处理))]^2), 4),
         区间89下 = round(均值差 - z89 * 差的标准误, 4),
         区间89上 = round(均值差 + z89 * 差的标准误, 4),
         P大于0 = round(pnorm(均值差 / 差的标准误), 4)) %>%
  select(对比, 均值差, 差的标准误, 区间89下, 区间89上, P大于0)
write_csv_fixed(pairs_df, file.path(RES, "03_组间对比.csv"))
print(pairs_df)

# ---------- 2) 我们的案例（二）：饱和曲线 ----------
t_hours <- 1:12
resid_fix <- c(0.5, -1.1, 0.9, -0.6, 1.2, -0.4, 0.7, -1.0, 0.3, 1.1, -0.8, 0.6)
satu <- tibble(时间 = t_hours, 厚度 = 58 * (1 - exp(-0.32 * t_hours)) + resid_fix)
write_csv_fixed(data.frame(时间 = satu$时间, 厚度 = round(satu$厚度, 6)),
                file.path(RES, "04_饱和曲线数据.csv"))

# 三种拟合：直线、三次多项式、线性样条（结点 4 与 8）
X_lin  <- cbind(1, satu$时间)
X_poly <- cbind(1, satu$时间, satu$时间^2, satu$时间^3)
X_spl  <- cbind(1, satu$时间, pmax(satu$时间 - 4, 0), pmax(satu$时间 - 8, 0))
f_lin <- ols_ls(X_lin,  satu$厚度)
f_poly <- ols_ls(X_poly, satu$厚度)
f_spl <- ols_ls(X_spl,  satu$厚度)

fit_df <- tibble(时间 = t_hours, 观测 = round(satu$厚度, 6),
                 直线 = round(as.numeric(X_lin %*% f_lin$coef), 4),
                 三次多项式 = round(as.numeric(X_poly %*% f_poly$coef), 4),
                 线性样条 = round(as.numeric(X_spl %*% f_spl$coef), 4))
write_csv_fixed(fit_df, file.path(RES, "05_三种拟合.csv"))

fit_summary <- tibble(
  模型 = c("直线", "三次多项式", "线性样条（结点 4、8）"),
  参数个数 = c(2, 4, 4),
  残差标准差 = round(c(f_lin$sigma_hat, f_poly$sigma_hat, f_spl$sigma_hat), 4)
)
write_csv_fixed(fit_summary, file.path(RES, "06_拟合优度.csv"))
print(fit_summary)

# 外推：把时间推到 16 小时，看三种模型怎么走
t_ext <- 13:16
ext_df <- tibble(时间 = t_ext,
                 直线 = round(as.numeric(cbind(1, t_ext) %*% f_lin$coef), 4),
                 三次多项式 = round(as.numeric(cbind(1, t_ext, t_ext^2, t_ext^3) %*% f_poly$coef), 4),
                 线性样条 = round(as.numeric(cbind(1, t_ext, pmax(t_ext - 4, 0), pmax(t_ext - 8, 0)) %*% f_spl$coef), 4))
write_csv_fixed(ext_df, file.path(RES, "07_外推四小时.csv"))
print(ext_df)

# 教材对照：年龄→体重，同样三种拟合（结点 10、15、25、40）
age <- as.numeric(k$age); wt <- as.numeric(k$weight)
Xl <- cbind(1, age); Xp <- cbind(1, age, age^2, age^3)
Xs <- cbind(1, age, pmax(age - 10, 0), pmax(age - 15, 0), pmax(age - 25, 0), pmax(age - 40, 0))
hw_fit <- tibble(
  模型 = c("直线", "三次多项式", "线性样条（结点 10/15/25/40）"),
  参数个数 = c(2, 4, 6),
  残差标准差 = round(c(ols_ls(Xl, wt)$sigma_hat, ols_ls(Xp, wt)$sigma_hat, ols_ls(Xs, wt)$sigma_hat), 4)
)
write_csv_fixed(hw_fit, file.path(RES, "08_教材对照_拟合优度.csv"))
print(hw_fit)

# ---------- 图 ----------
p1 <- ggplot(groups, aes(处理, 厚度)) +
  geom_point(size = 2.4, colour = "#c1462c", position = position_jitter(width = 0.08, seed = 1)) +
  geom_errorbar(data = means_df, aes(x = 组, ymin = 区间89下, ymax = 区间89上),
                inherit.aes = FALSE, width = 0.18, colour = "#163d48", linewidth = 0.9) +
  geom_point(data = means_df, aes(x = 组, y = 样本均值), inherit.aes = FALSE,
             size = 4, colour = "#167d80") +
  labs(title = "三个处理组的终点厚度：先把「组均值」的后验画出来",
       subtitle = "点：5 个重复；深色点与竖线：组均值的后验均值与 89% 区间（σ 固定为合并组内标准差）",
       x = NULL, y = "终点厚度（微米）")
save_fig("01-三组终点厚度.png", p1, 6.6, 4.0)

dens_x <- seq(-12, 14, length.out = 400)
dens_df <- do.call(rbind, lapply(seq_len(nrow(pairs_df)), function(i) {
  m <- pairs_df$均值差[i]; s <- pairs_df$差的标准误[i]
  tibble(x = dens_x, 密度 = dnorm(dens_x, m, s), 对比 = pairs_df$对比[i])
}))
dens_df$对比 <- factor(dens_df$对比, levels = pairs_df$对比)
p2 <- ggplot(dens_df, aes(x, 密度, colour = 对比)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#9aa8a8") +
  labs(title = "报系数不如报对比：三组两两差多少",
       subtitle = paste0("三条曲线是差值的后验分布；P(差值 > 0) 分别是 ",
                         paste(pairs_df$P大于0, collapse = "、")),
       x = "厚度差（微米）", y = "后验密度", colour = NULL)
save_fig("02-组间对比后验.png", p2, 7.4, 4.2)

# 中心化：同一份数据、同一根直线，两种写法下"截距"回答的是两个不同问题
tbar <- mean(satu$时间); n_obs <- length(satu$时间)
Sxx <- sum((satu$时间 - tbar)^2)
sd_a_raw <- f_lin$sigma_hat * sqrt(1 / n_obs + tbar^2 / Sxx)      # 未中心化：截距=时间 0 处
sd_a_ctr <- f_lin$sigma_hat / sqrt(n_obs)                          # 中心化：截距=平均时间处
ctr_df <- tibble(
  写法 = c("未中心化：截距＝时间 0 处的厚度", "中心化：截距＝平均时间处的厚度"),
  截距后验均值 = round(c(f_lin$coef[1], f_lin$coef[1] + f_lin$coef[2] * tbar), 4),
  截距后验标准差 = round(c(sd_a_raw, sd_a_ctr), 4)
)
write_csv_fixed(ctr_df, file.path(RES, "09_中心化对照.csv"))
print(ctr_df)
dens2 <- do.call(rbind, lapply(seq_len(nrow(ctr_df)), function(i) {
  tibble(x = dens_x + 40, 密度 = dnorm(dens_x + 40, ctr_df$截距后验均值[i], ctr_df$截距后验标准差[i]),
         写法 = ctr_df$写法[i])
}))
p3 <- ggplot(dens2, aes(x, 密度, fill = 写法)) +
  geom_area(alpha = 0.5) +
  labs(title = "中心化换的是「截距在回答哪个问题」",
       subtitle = "未中心化：截距＝时间 0 处（外推远、后验宽）；中心化：截距＝平均时间处（后验窄）",
       x = "截距（微米）", y = "后验密度", fill = NULL)
save_fig("03-中心化.png", p3, 7.4, 4.0)

t_fine <- seq(1, 16, length.out = 200)
curve_df <- tibble(
  时间 = t_fine,
  直线 = as.numeric(cbind(1, t_fine) %*% f_lin$coef),
  三次多项式 = as.numeric(cbind(1, t_fine, t_fine^2, t_fine^3) %*% f_poly$coef),
  线性样条 = as.numeric(cbind(1, t_fine, pmax(t_fine - 4, 0), pmax(t_fine - 8, 0)) %*% f_spl$coef)
) %>% pivot_longer(-时间, names_to = "模型", values_to = "拟合")
curve_df$模型 <- factor(curve_df$模型, levels = c("直线", "三次多项式", "线性样条"))
p4 <- ggplot() +
  geom_point(data = satu, aes(时间, 厚度), size = 2.4, colour = "#c1462c") +
  geom_line(data = filter(curve_df, 时间 <= 12), aes(时间, 拟合, colour = 模型), linewidth = 1) +
  labs(title = "直线不够用的时候：多项式与样条",
       subtitle = "1–12 小时里三次多项式贴得最紧（0.996），样条次之（1.46）；出范围后的走向见下图",
       x = "培养时间（小时）", y = "生物膜厚度（微米）", colour = NULL)
save_fig("04-曲线拟合.png", p4, 7.4, 4.2)

p5 <- ggplot() +
  geom_point(data = satu, aes(时间, 厚度), size = 2.4, colour = "#c1462c") +
  geom_vline(xintercept = 12, linetype = "dashed", colour = "#9aa8a8") +
  geom_line(data = curve_df, aes(时间, 拟合, colour = 模型), linewidth = 1) +
  labs(title = "外推四小时，三种模型给出三种未来",
       subtitle = "虚线右侧没有数据：直线按恒定斜率延伸、三次多项式弯得更快、样条沿最后一段斜率延伸",
       x = "培养时间（小时）", y = "生物膜厚度（微米）", colour = NULL)
save_fig("05-外推对照.png", p5, 7.4, 4.2)

cat("\n完成：5 张图 + 9 份 CSV 已写入\n")
