# Statistical Rethinking 精读04：按作者PPT修订 v3

先打开 `分类与曲线_手机预览_v3.html` 阅读；正文源文件是同名 v3 Markdown。各节标出对应PPT页序。公式均为PNG图片，SVG和生成代码一并提供。

本版以2023年第04讲PPT为核心：因果问题 → 模拟检验 → 索引编码 → 总体对比与个体预测 → 同身高处对比 → 多项式与局部样条 → 年龄—身高 → Bonus联合模型。

## 完整下载后怎么运行

保留本目录结构，尤其不要只下载主脚本。数据与样条设计矩阵已提供，三种语言无需依次运行。

### R / RStudio

打开 `精读04_v3.Rproj`，再打开 `代码/lesson04_v3.R`，点击 Source。结果写入 `运行结果/v3/`，六张中文图写入 `运行结果/v3/R_figures/`，交互环境的 Plots 面板可以前后翻看。脚本内部调用 `lesson04_plots_v3.R`。

依赖作者的 `rethinking` 包及其依赖，还有 `systemfonts`、`ragg`；`splines` 随R提供。安装作者包请遵循其官方说明：https://github.com/rmcelreath/rethinking 。其他依赖可在R中用 `install.packages(c("systemfonts","ragg"))` 安装。

R额外实现PPT第97页的完整联合模型，输出 `R_full.csv` 及 `R_full_contrasts.csv`。后者分别记录总的平均效应、两个独立模拟个体的差，不能混同为个体反事实效应。

### Python / Spyder

在实际使用的Python环境安装 `requirements_v3.txt` 中的 numpy、scipy、matplotlib。在Spyder中打开 `代码/lesson04_v3.py`，按 F5 运行完整文件。默认采用环境的交互绘图后端，不强制Agg；保存的六张公众号统计图位于 `文章配图_v3/`。

输出图片但没有窗口时，检查Spyder图形后端：可在Plots面板查看内嵌图，或选择自动/Qt后端查看独立窗口。无图形服务器只能保存图片。代码中 `plt.show()` 已保留。

`代码/lesson04_teaching_v3.py` 另行生成因果图、模拟规则示意、多项式形状图。这三张是辅助教学图；多项式采用最小二乘作函数形状演示，不宣称精确复现PPT中的后验。

### MATLAB

打开 `代码/lesson04_v3.m` 并 Run。脚本用自身位置定位项目，自动把绘图函数目录加入搜索路径。需要基础MATLAB，不使用统计工具箱；六张Figure保留，并导出到 `运行结果/v3/MATLAB_figures/`。

### 中文与公式

三种语言按本机可用字体选取 PingFang SC、Microsoft YaHei、Noto Sans CJK SC 等；缺少中文字体会明确报错。安装 Noto Sans CJK SC 后重启运行环境即可。源代码用UTF-8。

`代码/render_formula_v3.py` 是单独的公式导出工具，输出13张PNG和SVG。它使用无界面后端只为生成公式文件，不影响三语言统计绘图的GUI行为。

## 模型对应与范围

| 模型 | 对应来源 | 关键设定 |
|---|---|---|
| m1：每组平均体重 | PPT38—39 | α各自Normal(60,10)，σUniform(0,10) |
| m2：每组身高—体重直线 | PPT59 | α同上，β各自Uniform(0,1)，σ同上 |
| ms：年龄—身高样条 | PPT84—92；作者样条脚本后半段 | 基准Normal(120,1)，权重Normal(0,25)，logσ Normal(0,0.5) |
| full：身高与体重联合模型 | PPT97 | m2加组平均身高Normal(160,10)和τ Uniform(0,10) |

m1、m2、ms有R/Python/MATLAB三份独立计算；full目前由R实现。Python、MATLAB实现相同模型的二次近似，没有调用R quap。并非使用网格枚举全部参数。

R使用向量先验与显式初值以避免部分rethinking版本重设索引参数初值；数学模型与PPT的逐组先验相同。

总体和新个体差采用男性减女性；同身高处的差采用PPT61页女性减男性。不要跨行混读正负号。

## 数据与计算

`数据/v3/Howell1_all.csv`：作者rethinking包的544人；字段height（厘米）、weight（千克）、age（岁）、male（原数据0/1编码）。成年人文件筛选age≥18，共352人。分类模型重新编码S=male+1，即女性1、男性2。

年龄样条使用全部544人，20个等间距节点位置，按作者bs调用约定得到23列。`age_basis*.csv` 由R的 `splines::bs` 导出，三个语言共享同一矩阵，避免节点端点和截距约定差异。两端包含于knots的写法是作者脚本设定，不能擅自改成原书的15分位数节点。

R可重建这些CSV。模拟文件 `synthetic_case*.csv` 由教学图脚本生成，不是真实实验数据。随机预测不要求跨语言逐行一致。

已实跑三种语言并核对参数，详见 `VALIDATION_v3.md`。图形界面可用性依据原生绘图调用和保留的Figure检查；本次未逐个手动操作RStudio、Spyder、MATLAB桌面窗口。

## 原始资料

PPT：https://speakerdeck.com/rmcelreath/statistical-rethinking-2023-lecture-04

课程代码：https://github.com/rmcelreath/stat_rethinking_2023

数据与作者包：https://github.com/rmcelreath/rethinking

逐段页码对应与重绘差异：`PPT对应表_v3.md`。原书全文不包含在项目包中。
