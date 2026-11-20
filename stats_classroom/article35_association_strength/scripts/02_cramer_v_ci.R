# 统计课堂35｜关联强度与区间。输入为第11课的假设工程试验汇总计数。
# 300台独立记录器随机分3组各100台：A硅胶、B环氧、C热缩封装。
# 弯折100次后接模拟器记录10min，无持续≥1s数据缺口判合格（课堂阈值）。
# 主表按A/B/C行、合格/不合格列排列为[[90,10],[80,20],[65,35]]。
# 在项目根目录按节顺序运行；默认不写文件，ARTICLE35_SAVE_RESULTS=1才保存。

# 02.1 准备输入 ----
raw <- read.csv("data/device_pass.csv", stringsAsFactors=FALSE)
print(raw)
tab <- as.matrix(raw[, c("pass", "fail")])
stopifnot(all(rowSums(tab)>0), all(colSums(tab)>0))

# 02.2 调用函数计算V及双侧95%区间 ----
# cramers_v返回点估计与双侧区间；adjust=FALSE保留课堂未修正V。
fit_v <- as.data.frame(effectsize::cramers_v(tab, adjust=FALSE, ci=.95,
                                             alternative="two.sided"))

# 02.3 读取结果 ----
# V为整体关联强度，V_L/V_U为区间下限/上限，confidence记录置信水平。
V <- fit_v$Cramers_v
V_L <- fit_v$CI_low
V_U <- fit_v$CI_high
ci <- data.frame(V=V, CI_low=V_L, CI_high=V_U, confidence=.95)
print(ci)

# 02.4 效应分级 ----
interpretation <- data.frame(
  component=c("lower", "estimate", "upper"),
  V=c(V_L, V, V_U),
  label=as.character(effectsize::interpret_cramers_v(c(V_L, V, V_U),
                                                      rules="funder2019")),
  rule=rep("funder2019", 3), stringsAsFactors=FALSE)
stopifnot(all(interpretation$V >= 0 & interpretation$V <= 1))
print(interpretation)

# 02.5 本例解释与报告 ----
# chisq.test创建检验对象；report生成文字，report_table返回结构化结果。
chi_fit <- stats::chisq.test(tab, correct=FALSE)
report_result <- report::report(chi_fit, adjust=FALSE, alternative="two.sided", ci=.95)
report_table_result <- report::report_table(chi_fit, adjust=FALSE,
                                            alternative="two.sided", ci=.95)
report_v <- unname(report_table_result$Cramers_v)
report_low <- unname(report_table_result$Cramers_v_CI_low)
report_high <- unname(report_table_result$Cramers_v_CI_high)
stopifnot(abs(report_v-V)<1e-10, abs(report_low-V_L)<2e-7,
          abs(report_high-V_U)<2e-7)
cat(sprintf("本例Cramér V=%.6f，双侧95%%区间=[%.6f, %.6f]；区间表示估计不确定性。\n",
            V, V_L, V_U))
# 第1、2、3行分别是下限、点估计、上限的等级。
labels <- interpretation$label
cat(sprintf("按Funder–Ozer规则，点估计为%s；区间两端等级为%s/%s。\n",
            labels[2], labels[1], labels[3]))
print(report_result)

# 02.6 可选保存 ----
if (identical(Sys.getenv("ARTICLE35_SAVE_RESULTS"), "1")) {
  dir.create("results", showWarnings=FALSE)
  write.csv(ci, "results/R_ci.csv", row.names=FALSE)
  write.csv(interpretation, "results/R_interpretation.csv", row.names=FALSE)
}
