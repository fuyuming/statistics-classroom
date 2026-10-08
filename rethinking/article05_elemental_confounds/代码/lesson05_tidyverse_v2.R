# 【课堂问题】结婚率较高的地区，离婚率也较高。这种关联是否混入了结婚年龄的影响？
# 【数据】作者WaffleDivorce数据，49州和哥伦比亚特区，共50行。
# 【计算思路】①统一尺度；②写出似然与先验；③近似联合后验；④解释相同年龄下的斜率。
# 【第二个问题】为什么控制第三个变量，有时消除关联，有时反而制造关联？
# 【核对办法】按四种结构的生成规则枚举0/1状态，比较全体与分组后的相关。
# 建议先运行第一部分并读后验结果，再进入第二部分；数值输出和图形要结合当前问题阅读。

# 精读05 tidyverse版：RStudio打开精读05_v2.Rproj，再Source本文件。
# dplyr/tidyr/readr处理数据和结果，ggplot2绘图，rethinking::quap拟合贝叶斯模型。
# 作者PPT31页的quap模型；输入已附，无须先跑其他语言。
# 二元结构精确枚举，区别于PPT中1000次随机模拟的频数。
# 第一部分对应正文二至四：真实婚姻数据的贝叶斯回归。
# 第二部分对应四种因果结构：已知生成规则下的概率计算，不需要拟合后验。
library(rethinking)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
# Source时使用脚本位置；项目根目录中运行也可。
source_path <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
args <- commandArgs(FALSE)
arg <- grep("^--file=", args, value=TRUE)
if (length(arg)) source_path <- sub("^--file=", "", arg[1])
root <- if (!is.null(source_path)) dirname(dirname(normalizePath(source_path))) else getwd()
out <- file.path(root,"运行结果/v2/R_tidyverse")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read_csv(file.path(root,"数据/v2/WaffleDivorce.csv"), show_col_types=FALSE)
# 先读数据：每行一个州或特区；Divorce、Marriage单位为每千名成年人，
# MedianAgeMarriage是结婚年龄中位数（岁）。本项目CSV列名应保持一致。
# 换数据时先核对列名、单位与缺失值，不能把个体记录直接套入本例解释。
print(head(d))
# standardize减均值除样本标准差，与PPT一致。
# transmute只保留模型要用的三列，并把数据字段映射为公式符号。
dat <- d |> transmute(D=standardize(Divorce), M=standardize(Marriage),
                     A=standardize(MedianAgeMarriage)) |> as.list()
# D：离婚率，M：结婚率，A：结婚年龄，三者均已标准化。
# 似然 D ~ Normal(mu,sigma)：描述同一条件均值周围的数据波动。
# mu中的bM比较结婚年龄相同时的结婚率差异，bA反过来比较。
# a、bM、bA的正态先验约束合理尺度；sigma的指数先验只支持正数。
# quap寻找后验峰值，再用峰顶曲率给出联合正态近似；不是MCMC抽样。
# 参数怎么填：dnorm(均值,标准差)，不是填方差；0.2和0.5对应标准化尺度。
# dexp(1)的1是速率，均值为1/速率；不能把它理解为固定sigma=1。
# data=dat把公式中的D/M/A与数据对应。start只是优化起点，
# 依次给a、bM、bA、sigma初值；sigma必须为正，初值不决定先验均值。
# 调整先验应改dnorm/dexp中的参数，并重新检查先验预测；不要只改start。
# ① 写出模型：D的波动与其条件均值分别放在两行 --------------------
fit <- quap(
  alist(
    D ~ dnorm(mu, sigma),       # 似然：mu为各地区均值，sigma为剩余波动
    mu <- a + bM*M + bA*A,      # 固定A后，bM描述M变化对应的平均D变化
    a ~ dnorm(0, 0.2),          # 截距先验：均值0，标准差0.2
    bM ~ dnorm(0, 0.5),         # 结婚率斜率先验：均值0，标准差0.5
    bA ~ dnorm(0, 0.5),         # 结婚年龄斜率先验：同样尺度
    sigma ~ dexp(1)             # 指数先验：速率1，sigma仍是未知参数
  ),
  data = dat,                  # 使用上面整理好的D、M、A
  start = list(a=0, bM=0, bA=-0.5, sigma=0.8)  # 优化从这里开始
)

# ② 从联合后验近似中提取均值与不确定性 --------------------------
# 近似分布以峰值为中心：此处mean是近似后验均值。
# vcov保留参数间协方差，区间同时反映先验与数据的信息。
means <- coef(fit)[c("a","bM","bA","sigma")]
cv <- vcov(fit)[names(means),names(means)]
sds <- sqrt(diag(cv))  # 对角线是各参数的后验方差，开平方得到后验SD
# 区间怎么填：中央概率为c时，上分位点是(1+c)/2。
# 本例c=0.89，所以填0.945；若改95%，填0.975，并同步改列名和图注。
q <- qnorm(.945) # 中间89%，两端各5.5%
result <- tibble(parameter=names(means), mean=unname(means), sd=unname(sds)) |>
 mutate(lower89=mean-q*sd, upper89=mean+q*sd)
write_csv(select(result,-parameter),file.path(out,"posterior.csv"))
write_csv(result,file.path(out,"posterior_named.csv"))
# 结果怎么看：posterior.csv四行依次为a、bM、bA、sigma，sd是参数的后验标准差。
# 其中sigma本身是数据的残差标准差；sigma那一行的sd是对sigma估计的不确定性。
# bM约-0.065，89%区间约[-0.306,0.176]：年龄相同后，方向仍不确定；
# 区间跨0不证明作用为0。bA约-0.614，表示每增加1个年龄标准差，
# 平均离婚率降低约0.614个离婚率标准差，前提是结婚率相同。
print(result)  # 先看bM的区间，再看bA；对应正文图4
cat("bM是在年龄相同条件下的斜率；因果解释仍依赖图和模型假设。\n")
# 给出所有0/1组合，而不是用随机次数估计概率。
# 下面切换到已知概率的教学模型，X/Y/Z/A是0或1状态；这里的A不是结婚年龄。
# bern(value,p)的p填“该变量为1的概率”，应在0到1之间。
# 0.1+0.8*z：z=0时概率0.1，z=1时概率0.9；乘法来自生成规则的条件独立。
# 这只是枚举所有数据状态，不是在参数网格上近似后验。
# ③ 切换问题：生成规则已知，比较不同分组下的关联 ----------------
bern <- function(value,p) ifelse(value==1,p,1-p)
# model：1叉、2管、3对撞、4中间变量的后代。
# 前三种的A只是独立占位变量，乘0.5后求和即可消去。
# 每一行的概率沿箭头相乘，全部行加总为1。
# crossing列出全部可能状态，case_when按结构选用对应联合概率。
rows <- crossing(model=1:4,X=0:1,Z=0:1,Y=0:1,A=0:1) |>
 mutate(probability=case_when(
  model==1 ~ .5*bern(X,.1+.8*Z)*bern(Y,.1+.8*Z)*.5,
  model==2 ~ .5*bern(Z,.1+.8*X)*bern(Y,.1+.8*Z)*.5,
  model==3 ~ .25*bern(Z,if_else(X+Y>0,.9,.2))*.5,
  model==4 ~ .5*bern(Z,.1+.8*X)*bern(Y,.1+.8*Z)*bern(A,.1+.8*Z)))
# 复制为五种比较：全体、Z=0、Z=1、A=0、A=1。
# 分组筛选之后必须重新归一化，不能沿用全体的概率权重。
answers <- crossing(rows, group=0:4) |>
 filter(group==0 | (group==1 & Z==0) | (group==2 & Z==1) |
        (group==3 & A==0) | (group==4 & A==1)) |>
 group_by(model,group) |>
 mutate(w=probability/sum(probability)) |>
 summarise(ex=sum(w*X), ey=sum(w*Y), exy=sum(w*X*Y), .groups="drop") |>
 transmute(model,group,correlation=(exy-ex*ey)/sqrt(ex*(1-ex)*ey*(1-ey)))
stopifnot(all(abs((rows |> group_by(model) |> summarise(p=sum(probability)))$p-1)<1e-12))
write_csv(rows,file.path(out,"joint.csv"))
# 枚举输出怎么看：model=1/2/3/4分别为叉/管/对撞/后代。
# group=0全体、1固定Z=0、2固定Z=1、3固定A=0、4固定A=1。
# 后代模型全体相关约0.64，固定Z后为0，固定A后约0.390。
# 练习：只把后代模型最后的bern(A,0.1+0.8*Z)概率改成0.5
# （原版变量名是小写a、z）；固定A后的相关应回到0.64，其余机制保持不变。
write_csv(answers,file.path(out,"associations.csv"))
# 原书171页生成规则的期望计算，不是一次样本拟合。
# 植物参数：treatment的0/1分别是不处理/处理；0.5和0.1是出现真菌的概率。
# 5与2是无真菌/有真菌时的平均增长。按两种状态的概率加权，
# 得到3.5与4.7，两者之差1.2是生成规则给出的理论总效应。
plant <- tibble(treatment=0:1,fungus_probability=c(.5,.1)) |>
 mutate(mean_growth=(1-fungus_probability)*5+fungus_probability*2)
write_csv(plant,file.path(out,"plant_expectation.csv"))
print(plant)

# ggplot对象需要显式print，Source运行时才会在RStudio Plots中显示。
# ragg负责保存，字体检测避免生成中文方框；不会修改用户的图形后端。
library(systemfonts)
library(ragg)
fonts <- unique(system_fonts()$family)
font <- intersect(c("PingFang SC","Microsoft YaHei","Noto Sans CJK SC","Heiti SC"),fonts)[1]
if(is.na(font)) stop("请安装Noto Sans CJK SC中文字体。")
style <- theme_minimal(base_family=font,base_size=13)
# 图怎么读：点是近似后验均值，横线是89%后验区间，竖虚线为零。
# 区间宽度反映参数不确定性，不是地区之间的离婚率波动。
model_plot <- result |> filter(parameter %in% c("bM","bA")) |>
 mutate(label=if_else(parameter=="bM","结婚率","结婚年龄")) |>
 ggplot(aes(x=mean,y=label)) +
 geom_vline(xintercept=0,linetype=2,color="grey60") +
 geom_segment(aes(x=lower89,xend=upper89,yend=label),linewidth=1,color="#167e83") +
 geom_point(size=3,color="#167e83") +
 labs(x="标准化斜率及89%后验区间",y=NULL,title="年龄相同以后，结婚率的斜率接近零") + style
# 图10：理论相关按生成规则计算，不是婚姻模型的贝叶斯后验。
descendant_plot <- answers |> filter(model==4,group %in% c(0,1,3)) |>
 mutate(correlation=if_else(abs(correlation)<1e-12,0,correlation),
        label=factor(group,levels=c(0,1,3),labels=c("全体","固定中间变量Z","固定后代A"))) |>
 ggplot(aes(x=label,y=correlation)) + geom_col(fill="#167e83") +
 geom_text(aes(label=sprintf("%.3f",pmax(correlation,0))),vjust=-.6) +
 scale_y_continuous(limits=c(0,.75)) +
 labs(x=NULL,y="X与Y的理论相关",title="后代携带中间变量的信息") + style
for(name in c("model","descendant")) {
 plot <- get(paste0(name,"_plot"))
 if(interactive()) print(plot)
 ggsave(file.path(out,paste0(name,".png")),plot,device=ragg::agg_png,
        width=8,height=5,dpi=180,bg="white")
}
cat("已保存后验表、概率枚举表和两张中文图。\n")

# ④ 从斜率含义到两种情境的预测 ----------------------------------
# 同一行参数与年龄用于两种情境；均值差用于解释bM，观测差另加随机波动。
set.seed(20261009)
post <- extract.samples(fit, n = 10000)
post <- post[post$sigma > 0, ]
age_same <- sample(dat$A, nrow(post), replace = TRUE)
mu0 <- with(post, a + bM * 0 + bA * age_same)
mu1 <- with(post, a + bM * 1 + bA * age_same)
delta_mean <- mu1 - mu0
stopifnot(isTRUE(all.equal(delta_mean, post$bM)))
D0 <- rnorm(nrow(post), mean = mu0, sd = post$sigma)
D1 <- rnorm(nrow(post), mean = mu1, sd = post$sigma)
delta_observed <- D1 - D0
print(quantile(delta_mean, c(0.055, 0.5, 0.945)))
print(c(均值差的后验SD=sd(delta_mean), 两次观测之差的SD=sd(delta_observed)))
write.csv(data.frame(age_standardized=age_same,mu_M0=mu0,mu_M1=mu1,
                    mean_change=delta_mean,observed_change=delta_observed),
          file.path(out,"intervention_pairs.csv"),row.names=FALSE)
