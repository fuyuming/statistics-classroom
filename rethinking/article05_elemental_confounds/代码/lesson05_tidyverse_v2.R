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
fit <- quap(alist(D ~ dnorm(mu,sigma), mu <- a+bM*M+bA*A,
                 a ~ dnorm(0,0.2), bM ~ dnorm(0,0.5), bA ~ dnorm(0,0.5),
                 sigma ~ dexp(1)), data=dat, start=list(a=0,bM=0,bA=-0.5,sigma=0.8))
# 近似分布以峰值为中心：此处mean是近似后验均值。
# vcov保留参数间协方差，区间同时反映先验与数据的信息。
means <- coef(fit)[c("a","bM","bA","sigma")]
cv <- vcov(fit)[names(means),names(means)]
sds <- sqrt(diag(cv))
q <- qnorm(.945) # 中间89%，两端各5.5%
result <- tibble(parameter=names(means), mean=unname(means), sd=unname(sds)) |>
 mutate(lower89=mean-q*sd, upper89=mean+q*sd)
write_csv(select(result,-parameter),file.path(out,"posterior.csv"))
write_csv(result,file.path(out,"posterior_named.csv"))
print(result)
cat("bM是在年龄相同条件下的斜率；因果解释仍依赖图和模型假设。\n")
# 给出所有0/1组合，而不是用随机次数估计概率。
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
write_csv(answers,file.path(out,"associations.csv"))
# 原书171页生成规则的期望计算，不是一次样本拟合。
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
