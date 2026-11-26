# 一、读入两份独立教学数据；安装依赖：install.packages(c("car","nortest")) ----
options(stringsAsFactors = FALSE, digits = 17)
root <- getwd()
if (basename(root) == "代码") root <- dirname(root)
source(file.path(root,"代码","studentized_qq.R"))
x <- read.csv(file.path(root,"数据","normality_example.csv"))$value
d <- read.csv(file.path(root,"数据","anova_example.csv"))
d$group <- factor(d$group, levels = c("A1","A2","A3"))
print(x)
print(head(d))
alpha <- 0.05
save_outputs <- identical(Sys.getenv("RESIDPLOT_SAVE_OUTPUTS"),"1")
outdir <- file.path(root,"运行结果","R")
if (save_outputs) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 二、正态性：原始独立样本x；AD在nortest中要求至少8个观测。 ----
sw <- shapiro.test(x)
ad <- nortest::ad.test(x)
print(sw)
print(ad)
if (save_outputs) png(file.path(outdir,"normal_qq.png"),width = 1100,height = 800)
car::qqPlot(x, distribution = "norm", line = "quartiles", envelope = FALSE,
            id = FALSE, main = "Independent sample: normal Q-Q plot")
if (save_outputs) dev.off()

# 三、方差齐性：median中心的Levene即本课Brown-Forsythe修正。 ----
bartlett <- bartlett.test(value ~ group,data = d)
levene <- car::leveneTest(value ~ group,data = d,center = median)
statistics <- c(unname(sw$statistic),unname(ad$statistic),unname(bartlett$statistic),levene[1,"F value"])
pvalues <- c(sw$p.value,ad$p.value,bartlett$p.value,levene[1,"Pr(>F)"])
tests <- data.frame(method = c("Shapiro-Wilk","Anderson-Darling","Bartlett","Levene-median"),
                    sample = c("independent20","independent20","groups15","groups15"),
                    statistic = statistics,p_value = pvalues,alpha = alpha,reject = pvalues < alpha)
print(tests,row.names = FALSE)

# 四、car::qqPlot模型方法：外部学生化残差，t参照，关闭包络统一三语。 ----
fit <- lm(value ~ group,data = d)
if (save_outputs) png(file.path(outdir,"studentized_qq.png"),width = 1100,height = 800)
qq <- studentized_qq(fit)
if (save_outputs) dev.off()
if (save_outputs) png(file.path(outdir,"group_boxplot.png"),width = 1100,height = 800)
boxplot(value ~ group,data = d,xlab = "Group",ylab = "Compressive modulus (MPa)",
        main = "Group spread and observations",col = "grey88")
stripchart(value ~ group,data = d,vertical = TRUE,method = "stack",add = TRUE,pch = 19,col = "grey25")
if (save_outputs) {
  dev.off()
  write.csv(tests,file.path(outdir,"diagnostics_tests.csv"),row.names = FALSE,quote = FALSE)
  write.csv(qq$points,file.path(outdir,"studentized_qq.csv"),row.names = FALSE,quote = FALSE)
}
