# ============================================================
# 精读 08｜MCMC：网格走不动的时候怎么办
# 对应：McElreath 2023 第 08 讲 Markov Chain Monte Carlo；教材第 9 章
# 用法： Rscript 代码/篇08_MCMC.R
#
# 设计与前几篇的差别：MCMC 天生要用随机数，而本系列要求三个语言给出同一个答案。
# 这里用**确定性低差异序列**（Halton）代替伪随机数：提议步长由 qnorm(Halton) 给出，
# 接受判断用 Halton 的均匀序列。算法结构与真实 Metropolis 一致，但结果可复现、可核对。
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
wq <- function(v, w, p) {
  o <- order(v); v <- v[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  k <- which(cw >= p)[1]
  if (k == 1) return(v[1])
  v[k - 1] + (v[k] - v[k - 1]) * (p - cw[k - 1]) / (cw[k] - cw[k - 1])
}
# Halton 低差异序列（确定性"随机数"）
halton <- function(n, base) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    f <- 1; r <- 0; k <- i
    while (k > 0) { f <- f / base; r <- r + f * (k %% base); k <- k %/% base }
    out[i] <- r
  }
  out
}
orth_resid <- function(v, X) v - X %*% solve(t(X) %*% X, t(X) %*% v)

# ---------- 1) 数据：16 个观测，两参数模型 ----------
n <- 16
x <- 1:n
X0 <- cbind(1, x)
resid_y <- orth_resid(c(0.9, -1.1, 0.7, -0.8, 1.2, -0.4, 0.6, -1.0, 0.5, 1.1, -0.7, 0.8, -0.9, 0.4, -0.6, 1.0),
                      cbind(1, x))
y <- 12 + 2.5 * x + resid_y
write_csv_fixed(data.frame(单元 = 1:n, 时间 = x, 观测 = r4(y)), file.path(RES, "01_数据.csv"))

# 似然：σ 固定为残差标准差（与第 03、04 篇同一简化，便于对照）
s_y <- sqrt(sum((y - (cbind(1, x) %*% solve(t(X0) %*% X0, t(X0) %*% y)))^2) / (n - 2))
log_lik <- function(a, b) -0.5 * sum(((y - (a + b * x)) / s_y)^2)

# ---------- 2) 网格后验（两参数，供对照）----------
a_seq <- 0 + (0:299) * (24 / 299)
b_seq <- 0 + (0:299) * (5 / 299)
grid_df <- expand.grid(a = a_seq, b = b_seq)
lp <- apply(grid_df, 1, function(r) log_lik(r[1], r[2]))
grid_df$post <- exp(lp - max(lp)); grid_df$post <- grid_df$post / sum(grid_df$post)
grid_sum <- data.frame(
  方法 = "网格（300×300 = 9 万格点）",
  截距均值 = r4(weighted.mean(grid_df$a, grid_df$post)),
  截距89下 = r4(wq(grid_df$a, grid_df$post, 0.055)),
  截距89上 = r4(wq(grid_df$a, grid_df$post, 0.945)),
  斜率均值 = r4(weighted.mean(grid_df$b, grid_df$post)),
  斜率89下 = r4(wq(grid_df$b, grid_df$post, 0.055)),
  斜率89上 = r4(wq(grid_df$b, grid_df$post, 0.945))
)

# ---------- 3) Metropolis（确定性低差异序列）----------
metropolis <- function(start, n_step, step_a, step_b, offset) {
  chain <- matrix(NA_real_, n_step, 2)
  chain[1, ] <- start
  z1 <- qnorm(halton(n_step, 2)); z2 <- qnorm(halton(n_step, 3)); u <- halton(n_step, 5)
  lp_cur <- log_lik(chain[1, 1], chain[1, 2])
  n_accept <- 0
  for (i in 2:n_step) {
    j <- ((i + offset - 2) %% n_step) + 1
    prop <- chain[i - 1, ] + c(step_a * z1[j], step_b * z2[j])
    lp_prop <- log_lik(prop[1], prop[2])
    if (log(u[j]) < (lp_prop - lp_cur)) {
      chain[i, ] <- prop; lp_cur <- lp_prop; n_accept <- n_accept + 1
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
  }
  list(chain = chain, accept = n_accept / (n_step - 1))
}
n_step <- 20000; burn <- 2000
# 三组步长：太大（几乎全被拒绝）、合适、太小（漫游太慢）
ch1  <- metropolis(c(3, 0.5), n_step, 0.9, 0.9, 0)        # 步长太大
ch2  <- metropolis(c(3, 0.5), n_step, 0.25, 0.025, 0)     # 合适
ch3  <- metropolis(c(3, 0.5), n_step, 0.05, 0.005, 0)     # 太小
ch2b <- metropolis(c(20, 4.0), n_step, 0.25, 0.025, 777)  # 第二条链（起点很远）
keep1 <- ch1$chain[(burn + 1):n_step, ]
keep2 <- ch2$chain[(burn + 1):n_step, ]
keep3 <- ch3$chain[(burn + 1):n_step, ]
keep2b <- ch2b$chain[(burn + 1):n_step, ]
all_keep <- rbind(keep2, keep2b)          # 用"合适步长"的两条链做后验汇总

# 自相关（滞后 1–20）：用合并后的两条链
acf_of <- function(v, lag) {
  v <- v - mean(v)
  sum(v[(lag + 1):length(v)] * v[1:(length(v) - lag)]) /
    sum(v^2)
}
acf_tab <- data.frame(
  滞后 = 1:20,
  步长太大 = r4(sapply(1:20, function(l) acf_of(keep1[, 2], l))),
  合适步长 = r4(sapply(1:20, function(l) acf_of(keep2[, 2], l))),
  步长太小 = r4(sapply(1:20, function(l) acf_of(keep3[, 2], l)))
)

# 收敛：运行均值；简化 R-hat（两条"合适步长"链的方差比）
run_mean1 <- cumsum(keep2[, 2]) / seq_along(keep2[, 2])
run_mean2 <- cumsum(keep2b[, 2]) / seq_along(keep2b[, 2])
w <- mean(c(var(keep2[, 2]), var(keep2b[, 2])))
b <- var(c(mean(keep2[, 2]), mean(keep2b[, 2])))
r_hat <- sqrt((((nrow(keep2) - 1) / nrow(keep2)) * w + b / nrow(keep2)) / w)
n_eff <- function(v) length(v) * (1 - acf_of(v, 1)) / (1 + acf_of(v, 1))
mcmc_sum <- data.frame(
  方法 = "Metropolis（两条链，预热 2000 步，各留 18000 步）",
  截距均值 = r4(mean(all_keep[, 1])),
  截距89下 = r4(wq(all_keep[, 1], rep(1 / nrow(all_keep), nrow(all_keep)), 0.055)),
  截距89上 = r4(wq(all_keep[, 1], rep(1 / nrow(all_keep), nrow(all_keep)), 0.945)),
  斜率均值 = r4(mean(all_keep[, 2])),
  斜率89下 = r4(wq(all_keep[, 2], rep(1 / nrow(all_keep), nrow(all_keep)), 0.055)),
  斜率89上 = r4(wq(all_keep[, 2], rep(1 / nrow(all_keep), nrow(all_keep)), 0.945))
)
cmp_tab <- rbind(grid_sum, mcmc_sum)
write_csv_fixed(cmp_tab, file.path(RES, "02_网格与MCMC对照.csv"))
print(cmp_tab)

diag_tab <- data.frame(
  量 = c("步长太大：接受率", "合适步长：接受率", "步长太小：接受率",
         "步长太大：滞后 1 自相关", "合适步长：滞后 1 自相关", "步长太小：滞后 1 自相关",
         "合适步长：近似有效样本量（3.6 万步）", "简化 R-hat（斜率，两条合适步长的链）",
         "斜率后验均值（链 1）", "斜率后验均值（链 2）"),
  值 = r4(c(ch1$accept, ch2$accept, ch3$accept,
            acf_tab$步长太大[1], acf_tab$合适步长[1], acf_tab$步长太小[1],
            n_eff(keep2[, 2]), r_hat, mean(keep2[, 2]), mean(keep2b[, 2])))
)
write_csv_fixed(diag_tab, file.path(RES, "03_诊断.csv"))
print(diag_tab)

# 维度诅咒：每维 100 格，格点数随参数个数增长
dim_tab <- data.frame(参数个数 = 1:8, 格点数 = r4(100^(1:8)))
write_csv_fixed(dim_tab, file.path(RES, "04_维度诅咒.csv"))

run_tab <- data.frame(步数 = seq(100, nrow(keep1), by = 500),
                      链1运行均值 = r4(run_mean1[seq(100, nrow(keep1), by = 500)]),
                      链2运行均值 = r4(run_mean2[seq(100, nrow(keep2), by = 500)]))
write_csv_fixed(run_tab, file.path(RES, "05_运行均值.csv"))

# ---------- 图 ----------
p1 <- ggplot() +
  geom_raster(data = grid_df, aes(a, b, fill = post), alpha = 0.85, interpolate = TRUE) +
  geom_point(data = data.frame(a = all_keep[seq(1, nrow(all_keep), by = 20), 1],
                               b = all_keep[seq(1, nrow(all_keep), by = 20), 2]),
             aes(a, b), colour = "#c1462c", size = 0.5, alpha = 0.5) +
  scale_fill_gradient(low = "#eef5f5", high = "#167d80", guide = "none") +
  labs(title = "MCMC 做的事：在网格画出的后验上走一圈",
       subtitle = paste0("底色＝网格算出的后验（对照用）；红点＝确定性 Metropolis 链抽样 1/20 的点（",
                         n_step * 2 - 2 * burn, " 步中的一部分）"),
       x = "截距 a", y = "斜率 b")
save_fig("01_MCMC与网格.png", p1, 7.4, 4.6)

trace_df <- rbind(
  data.frame(步数 = 1:1200, 链 = "链 1（起点 a=3, b=0.5）", 斜率 = ch2$chain[1:1200, 2]),
  data.frame(步数 = 1:1200, 链 = "链 2（起点 a=20, b=4）", 斜率 = ch2b$chain[1:1200, 2])
)
p2 <- ggplot(trace_df, aes(步数, 斜率, colour = 链)) +
  geom_line(linewidth = 0.45) +
  geom_vline(xintercept = burn, linetype = "dashed", colour = "#9aa8a8") +
  labs(title = "采样轨迹：两条链从很远的起点走到一起",
       subtitle = "虚线是预热结束的位置；看轨迹是不是在同一片区域上下抖动、有没有长期趋势",
       x = "步数", y = "斜率 b", colour = NULL)
save_fig("02_轨迹.png", p2, 7.4, 4.4)

p3 <- ggplot(acf_tab %>% pivot_longer(-滞后, names_to = "步长", values_to = "自相关"),
             aes(滞后, 自相关, colour = 步长)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 0, colour = "#9aa8a8") +
  labs(title = "自相关：相邻样本有多像",
       subtitle = "步长太大（几乎全被拒绝）时自相关居高不下；合适步长才让链更快「忘记」起点",
       x = "滞后", y = "自相关", colour = NULL) +
  ylim(0, 1)
save_fig("03_自相关.png", p3, 7.4, 4.4)

p4 <- ggplot(run_tab %>% pivot_longer(-步数, names_to = "链", values_to = "运行均值"),
             aes(步数, 运行均值, colour = 链)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = mcmc_sum$斜率均值, linetype = "dashed", colour = "#c1462c") +
  labs(title = "运行均值：两条链最后收敛到同一个数",
       subtitle = "红线＝合并两条链后的后验均值；两条链前期差得多，后期贴合",
       x = "步数", y = "斜率的运行均值", colour = NULL)
save_fig("04_运行均值.png", p4, 7.4, 4.4)

p5 <- ggplot(dim_tab, aes(factor(参数个数), 格点数)) +
  geom_col(fill = "#c1462c", width = 0.6) +
  scale_y_log10() +
  geom_text(aes(label = format(格点数, scientific = TRUE, digits = 2)), vjust = -0.4, size = 3.2) +
  labs(title = "为什么网格走不动：每维 100 格，格点数按参数个数指数增长",
       subtitle = "纵轴取对数；6 个参数就已经是 10¹² 个格点，8 个参数是 10¹⁶",
       x = "参数个数", y = "格点数（对数刻度）") +
  ylim(1, 1e17)
save_fig("05_维度诅咒.png", p5, 7.4, 4.4)

cat("\n完成：5 张图 + 5 份 CSV 已写入\n")
