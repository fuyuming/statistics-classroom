# 统计课堂35｜关联强度与区间。输入为第11课的假设工程试验汇总计数。
# 300台独立记录器随机分3组各100台：A硅胶、B环氧、C热缩封装。
# 弯折100次后接模拟器记录10min，无持续≥1s数据缺口判合格（课堂阈值）。
# 主表按A/B/C行、合格/不合格列排列为[[90,10],[80,20],[65,35]]。
# 在项目根目录按节顺序运行；默认不写文件，ARTICLE35_SAVE_RESULTS=1才保存。

# 01.1 准备与输入 ----
raw <- read.csv("data/device_pass.csv", stringsAsFactors = FALSE)
print(raw)
tab <- as.matrix(raw[, c("pass", "fail")])
stopifnot(all(rowSums(tab) > 0), all(colSums(tab) > 0))
n <- sum(tab)
r <- nrow(tab)
s <- ncol(tab)

# 01.2 手算预期频数、Q、P、C与V ----
# outer计算各行合计与各列合计的乘积，除以总数得到预期频数。
expected <- outer(rowSums(tab), colSums(tab)) / n
Q <- sum((tab - expected)^2 / expected)
df <- (r - 1) * (s - 1)
m <- min(r - 1, s - 1)
p_value <- stats::pchisq(Q, df, lower.tail = FALSE)  # 读取χ²右尾
C <- sqrt(Q / (Q + n))
Cmax <- sqrt(m / (m + 1))
V <- sqrt(Q / (n*m))
point <- data.frame(n=n,r=r,s=s,df=df,m=m,Q=Q,p=p_value,C=C,Cmax=Cmax,V=V)
print(expected)
print(point)

# 01.3 A/C子表：有符号phi与0/1相关 ----
ac <- tab[c(1, 3), ]                  # R从1计数，第1和第3行是A/C
n_sub <- sum(ac)
a <- ac[1, 1]
b <- ac[1, 2]
c_count <- ac[2, 1]
d <- ac[2, 2]
phi <- (a*d - b*c_count) / sqrt(prod(rowSums(ac))*prod(colSums(ac)))
# 按原四格人数展开类别标签；用于复核0/1相关的代数等价。
x <- rep(c(1,1,0,0), c(a,b,c_count,d))
y <- rep(c(1,0,1,0), c(a,b,c_count,d))
phi_pearson <- stats::cor(x, y)
flipped <- ac[2:1, ]
phi_flipped <- det(flipped) / sqrt(prod(rowSums(flipped))*prod(colSums(flipped)))
Q_sub <- unname(stats::chisq.test(ac, correct=FALSE)$statistic)
V_sub <- sqrt(Q_sub/n_sub)
phi_result <- data.frame(a=a,b=b,c=c_count,d=d,n_sub=n_sub,Q_sub=Q_sub,
                         phi=phi,phi_pearson=phi_pearson,
                         phi_row_flipped=phi_flipped,V_sub=V_sub)
stopifnot(abs(phi-phi_pearson)<1e-12, abs(phi_flipped+phi)<1e-12)
print(phi_result)

# 01.4 相同比例与不同样本量的假设试验 ----
scaling <- data.frame()
for (k in c(1, 2, 4)) {
  enlarged <- k * tab
  nn <- sum(enlarged)
  ee <- outer(rowSums(enlarged), colSums(enlarged)) / nn
  qq <- sum((enlarged-ee)^2/ee)
  row <- data.frame(multiplier=k,n=nn,Q=qq,
                    p=stats::pchisq(qq,df,lower.tail=FALSE),Q_over_n=qq/nn,
                    C=sqrt(qq/(qq+nn)),V=sqrt(qq/(nn*m)))
  scaling <- rbind(scaling, row)
}
stopifnot(all(abs(scaling$V-V)<1e-12))
print(scaling)

# 01.5 完整软件核对与可选保存 ----
# chisq.test提供Q/P；effectsize分别提供Pearson C与Cramér V；ci=NULL仅返回点估计。
fit <- stats::chisq.test(tab, correct=FALSE)
# adjust=FALSE关闭Cramér V的小样本偏差修正，保持课堂公式；ci=NULL仅返回点估计。
fit_c <- effectsize::pearsons_c(tab, ci=NULL)
fit_v_point <- effectsize::cramers_v(tab, adjust=FALSE, ci=NULL)
function_comparison <- data.frame(quantity=c("Q","P","C","V"),
  manual=c(Q,p_value,C,V),
  software=c(unname(fit$statistic),fit$p.value,
             fit_c$Pearsons_c,fit_v_point$Cramers_v))
stopifnot(max(abs(function_comparison$manual-function_comparison$software))<1e-10)
print(function_comparison)
# rstatix独立核对V；correct=FALSE关闭Yates连续性修正。
V_rstatix <- rstatix::cramer_v(tab, correct=FALSE)
stopifnot(abs(V_rstatix-V)<1e-10)
print(V_rstatix)
if (identical(Sys.getenv("ARTICLE35_SAVE_RESULTS"), "1")) {
  dir.create("results", showWarnings=FALSE)
  write.csv(point, "results/R_point.csv", row.names=FALSE)
  write.csv(phi_result, "results/R_phi.csv", row.names=FALSE)
  write.csv(scaling, "results/R_scaling.csv", row.names=FALSE)
}
