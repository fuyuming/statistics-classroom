# 明哥的微生物世界 · 统计课堂 06
# 指数分布与威布尔分布：逐题计算、解释与配图
# 所有参数均为教学设定；不含真实寿命、微生物或学生个体数据。
# 文件保存为 UTF-8。RStudio 可用 File > Reopen with Encoding > UTF-8。
# 运行方式：打开本文件，点 Source；或 Rscript 代码/exponential_weibull.R。
# 计算仅用基础 R；绘图需 install.packages(c("ggplot2", "showtext"))。
# 在下面设置项目文件夹（包含“代码”文件夹），仅自动定位失败时需要。
project_dir <- ""

# ---------- 0. 定位输出目录 ----------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!nzchar(project_dir)) {
  if (length(file_arg)) project_dir <- dirname(dirname(normalizePath(file_arg[1])))
  else if (!is.null(source_file)) project_dir <- dirname(dirname(normalizePath(source_file)))
  else if (dir.exists("代码")) project_dir <- getwd()
  else if (basename(getwd()) == "代码") project_dir <- dirname(getwd())
  else stop("请用 Source 运行完整文件，或把 project_dir 设置为本课文件夹的完整路径。")
}
fig_dir <- file.path(project_dir, "文章配图")
out_dir <- file.path(project_dir, "运行结果")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# split=TRUE 同时把结果显示在控制台并存为文本；不上传任何个人信息。
sink(file.path(out_dir, "逐题计算结果.txt"), split = TRUE)
tryCatch({
section <- function(title) cat("\n", strrep("=", 48), "\n", title, "\n", sep = "")
answer <- function(label, value) cat(sprintf("%s：%.6f（%.2f%%）\n", label, value, 100 * value))

section("1. 泊松数次数，指数量等待：同一个事件")
lambda <- 0.5  # 每小时平均 0.5 通；这是率，单位为 1/小时。
t <- 4
cat("设定：来电服从齐次泊松过程，恒定发生率 0.5/小时。\n")
cat(sprintf("4 小时内来电数的期望 = 0.5 × 4 = %.1f 通。\n", lambda*t))
p_zero <- dpois(0, lambda = lambda*t)
# lower.tail=FALSE 直接算右尾 P(T>t)；默认 TRUE 算 P(T<=t)。
p_wait <- pexp(t, rate = lambda, lower.tail = FALSE)
answer("4 小时零来电 P(N(4)=0)", p_zero)
answer("等待超过 4 小时 P(T>4)", p_wait)
cat("解释：这两句话描述同一个事件，所以概率相等。平均等待为 2 小时，不是固定每 2 小时来电。\n")

section("2. 无记忆性：已经等 2 小时，再等 1 小时")
s <- 2; u <- 1
num <- pexp(s+u, lambda, lower.tail = FALSE)
den <- pexp(s, lambda, lower.tail = FALSE)
p_cond <- num/den
answer("从起点算，等超过 3 小时", num)
answer("分母：已经等超过 2 小时", den)
answer("已等 2 小时后，还要再等超过 1 小时", p_cond)
answer("刚开始时，等超过 1 小时", pexp(u, lambda, lower.tail=FALSE))
cat("解释：0.606531 是在已等过 2 小时的情形中计算的比例；0.223130 是从起点计算的比例。\n")

section("3. 教材猫寿命例：先接受题设，再做条件概率")
cat("题设：总寿命 X 从出生起算，服从 rate=0.1/年的指数分布。\n")
cat("这是用于练习的模型假设，不是真实猫寿命的研究结论。\n")
lambda_cat <- 0.1
cat(sprintf("模型均值：%.2f 年；中位数：%.2f 年。\n", 1/lambda_cat, qexp(.5,lambda_cat)))
answer("一年内发生事件的概率（并非直接用 0.1）", pexp(1,lambda_cat))
answer("从出生算，活过 15 年", pexp(15,lambda_cat,lower.tail=FALSE))
p_cat <- pexp(18,lambda_cat,lower.tail=FALSE)/pexp(10,lambda_cat,lower.tail=FALSE)
answer("已活过 10 年，再活过 8 年 = S(18)/S(10)", p_cat)
answer("用于比较：从出生算，活过 8 年", pexp(8,lambda_cat,lower.tail=FALSE))
answer("另一事件：从出生算，在 10 到 18 年间死亡", pexp(18,lambda_cat)-pexp(10,lambda_cat))
cat("解释：最后一项不是题目所求。年龄 10 岁是条件，不能拿它来估计 rate。\n")

section("4. 威布尔：让风险率随时间改变")
eta <- 10  # 尺度，单位为年；除 shape=1 外一般不是均值。
k <- c(0.5, 1, 2)
# 用对数生存概率相减，再指数化，避免尾部很小的概率直接相除。
conditional_survival <- function(age, extra, shape, scale) {
  exp(pweibull(age+extra, shape, scale, lower.tail=FALSE, log.p=TRUE) -
        pweibull(age, shape, scale, lower.tail=FALSE, log.p=TRUE))
}
comparison <- data.frame(
  shape=k, scale_years=eta,
  mean_years=eta*gamma(1+1/k),
  survive_8_from_start=pweibull(8,k,eta,lower.tail=FALSE),
  survive_8_after_age10=conditional_survival(10,8,k,eta)
)
print(comparison, row.names=FALSE, digits=6)
cat("解释：固定的是尺度 10 年，三个模型均值不同。这是比较模型性质，未用数据拟合。\n")
cat("shape<1：风险率下降；shape=1：恒定；shape>1：上升。\n")
cat("尺度为 10 年表示到第 10 年累计事件概率约 63.21%，不表示平均寿命必为 10 年。\n")
write.csv(comparison,file.path(out_dir,"weibull_comparison.csv"),row.names=FALSE,fileEncoding="UTF-8")

section("5. 小练习：试着先手算，再核对")
answer("发生率 2/小时：等待超过半小时", pexp(.5,2,lower.tail=FALSE))
answer("发生率 2/小时：已经等 1 小时，再等超过半小时", pexp(1.5,2,lower.tail=FALSE)/pexp(1,2,lower.tail=FALSE))
answer("威布尔 shape=2、scale=10 年：已用 10 年，再用超过 2 年", conditional_survival(10,2,2,10))

# 定义核对：避免只运行代码而没有判断结果。
grid_t <- seq(0,30,length.out=301)
stopifnot(abs(p_zero-p_wait)<1e-12,
          abs(p_cond-exp(-lambda*u))<1e-12,
          abs(p_cat-exp(-.8))<1e-12,
          max(abs(pweibull(grid_t,1,10)-pexp(grid_t,.1)))<1e-12,
          comparison$survive_8_after_age10[1] > comparison$survive_8_from_start[1],
          comparison$survive_8_after_age10[3] < comparison$survive_8_from_start[3])
cat("\n公式、条件概率及 Weibull(shape=1)=Exponential 的核对全部通过。\n")
}, finally = sink())

# ---------- 6. 可选绘图：中文乱码与字体 ----------
# 编码管文字是否读对；字体管图里是否能画出中文字。这是两件事。
if (!requireNamespace("ggplot2",quietly=TRUE) || !requireNamespace("showtext",quietly=TRUE)) {
  message("计算已完成；安装 ggplot2、showtext 后重新运行即可生成中文配图。")
} else {
  # 可自行指定中文字体文件，例如 C:/Windows/Fonts/msyh.ttc。
  font_file <- ""
  candidates <- c("C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/simhei.ttf",
                  path.expand("~/Library/Fonts/MSYH.TTC"),
                  "/System/Library/Fonts/PingFang.ttc",
                  Sys.glob("/System/Library/AssetsV2/com_apple_MobileAsset_Font*/*/AssetData/PingFang.ttc"),
                  "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
  if (!nzchar(font_file)) font_file <- candidates[file.exists(candidates)][1]
  if (is.na(font_file) || !file.exists(font_file)) {
    message("计算已完成；未找到中文字体。请安装思源黑体/Noto Sans CJK，并填写 font_file 后重跑绘图。")
  } else {
    library(ggplot2)
    sysfonts::font_add("chinese", regular=font_file)
    showtext::showtext_opts(dpi=180)
    showtext::showtext_auto()
    cols <- c("#C8793A", "#187F7B", "#5363AD")
    theme_set(theme_minimal(base_size=15, base_family="chinese") +
      theme(plot.background=element_rect(fill="#FAF9F5",colour=NA),
            panel.grid.minor=element_blank(),
            plot.title=element_text(face="bold",size=20),
            plot.subtitle=element_text(size=12,margin=margin(b=15)),
            plot.caption=element_text(size=10,hjust=0),
            legend.position="top",legend.title=element_blank(),
            plot.margin=margin(20,24,20,20)))
    save_plot <- function(name,p,height=5.4) ggsave(file.path(fig_dir,name),p,width=9,height=height,dpi=180)
    tt <- seq(0,10,length.out=501)
    d <- data.frame(t=tt,density=dexp(tt,.5),survival=pexp(tt,.5,lower.tail=FALSE))
    p1 <- ggplot(d,aes(t,density)) +
      geom_area(data=subset(d,t>=4),fill="#C8793A",alpha=.7) +
      geom_line(colour="#187F7B",linewidth=1.2) +
      geom_vline(xintercept=4,linetype="dashed",colour="#C8793A") +
      annotate("text",x=6.8,y=.29,label="等超过 4 小时\n= 前 4 小时零来电\n= exp(−0.5 × 4) ≈ 13.53%",family="chinese",size=5,colour="#263B44")+
      labs(title="数次数与量等待，怎样连起来？",subtitle="齐次泊松过程：每小时平均 0.5 通来电",x="从开始等待算起的时间（小时）",y="概率密度（1/小时）",
           caption="橙色是右尾面积（右侧延伸至无穷），不是曲线高度。图中仅显示到第 10 小时。\n参数为教学设定；明哥的微生物世界 · 统计课堂 06")
    save_plot("01_次数与等待.png",p1)
    d2 <- rbind(data.frame(t=tt,p=exp(-.5*tt),curve="刚开始：P(T > u)"),
                data.frame(t=tt,p=exp(-.5*(2+tt))/exp(-.5*2),curve="已等 2 小时：P(T > 2+u | T > 2)"))
    p2 <- ggplot(d2,aes(t,p,colour=curve,linetype=curve))+
      geom_line(linewidth=1.25)+scale_colour_manual(values=c("#187F7B","#C8793A"))+
      scale_linetype_manual(values=c("solid","dashed"))+
      geom_point(data=data.frame(t=1,p=exp(-.5)),aes(t,p),inherit.aes=FALSE,size=3,colour="#263B44")+
      annotate("text",x=4.5,y=.65,label="两条曲线处处重合\n再等超过 1 小时的概率均为 60.65%",family="chinese",size=5)+
      labs(title="无记忆性，看的是剩余等待时间",subtitle="先限定“已经等过 2 小时”，再比较未来还要等多久",x="从此刻起，还要等的时间 u（小时）",y="剩余等待超过 u 的概率",
           caption="无记忆性来自指数模型的恒定风险率；不是所有等待都具有此性质。")
    save_plot("02_无记忆性.png",p2)
    tw <- seq(.1,25,length.out=500)
    d3 <- do.call(rbind,lapply(k,function(kk) data.frame(t=tw,shape=paste0("k = ",kk),hazard=kk/eta*(tw/eta)^(kk-1),survival=pweibull(tw,kk,eta,lower.tail=FALSE))))
    dd <- rbind(data.frame(t=d3$t,shape=d3$shape,value=d3$hazard,quantity="风险率 h(t)（1/年）"),data.frame(t=d3$t,shape=d3$shape,value=d3$survival,quantity="生存概率 S(t)"))
    dd$quantity <- factor(dd$quantity, levels=c("风险率 h(t)（1/年）", "生存概率 S(t)"))
    p3 <- ggplot(dd,aes(t,value,colour=shape))+geom_line(linewidth=1.1)+
      facet_wrap(~quantity,ncol=1,scales="free_y")+scale_colour_manual(values=cols)+
      labs(title="威布尔多一个形状参数，风险就可以变化",subtitle="k < 1 下降；k = 1 恒定；k > 1 上升——说的是风险率",x="从起点算起的时间（年）",y=NULL,
           caption="均固定尺度 η = 10 年，均值并不相同；全部为教学模型。\nk < 1 时风险率在 t → 0 时趋于无穷，图从 0.1 年开始绘制。")
    save_plot("03_威布尔风险率与生存.png",p3,height=8)
    ages <- seq(0,20,length.out=200)
    d4 <- do.call(rbind,lapply(k,function(kk) data.frame(age=ages,p=conditional_survival(ages,8,kk,eta),shape=paste0("k = ",kk))))
    labels <- data.frame(age=10,p=comparison$survive_8_after_age10,shape=paste0("k = ",k),label=sprintf("%.1f%%",comparison$survive_8_after_age10*100))
    p4 <- ggplot(d4,aes(age,p,colour=shape))+geom_line(linewidth=1.2)+scale_colour_manual(values=cols)+
      geom_vline(xintercept=10,linetype="dotted",colour="#777777")+
      geom_point(data=labels,size=3)+geom_text(data=labels,aes(label=label),show.legend=FALSE,hjust=-.3,vjust=-.7,family="chinese",size=4.5)+
      coord_cartesian(ylim=c(0,1))+
      labs(title="已经用得更久，再用 8 年的机会怎样变？",subtitle="同一模型内比较不同已使用时间；尺度 η = 10 年",x="已经使用且尚未发生事件的时间 s（年）",y="再保持超过 8 年无事件的条件概率",
           caption="计算的是 S(s+8)/S(s)。只有 k = 1 的水平线具有无记忆性。\n模型间均值不同；曲线描述假设，不能单凭形状判定实际老化或筛选机制。")
    save_plot("04_剩余时间的变化.png",p4)
    cat("\n4 张中文配图已生成：",fig_dir,"\n")
  }
}
