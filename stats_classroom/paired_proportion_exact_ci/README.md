# 统计课堂28：配对率差的四种区间

三语入口使用同一份20人配对资料（before/after均为0/1）：两次阳性8、阴转阳1、阳转阴5、两次阴性6，`n=20`，后测减前测 `D=-.20`。

## 从哪里开始

打开项目后从各入口第一节开始运行：准备节只定位项目与设置保存开关，问题1读取字段、前3行并计算配对计数和`D`，问题2计算当前`confidence`下的区间，Python/MATLAB的问题3比较固定的90%、95%、99%并画图；R的问题3展示当前置信水平下的四种区间。R入口是完整多法比较：`ratesci::rdpairci`同时给出 Wald、MOVER-NW、AS(Tango)，`ExactCIdiff::PairedCI`另算 Wang；Python/MATLAB只用各自原生函数复算 Wald。三者都读取共享CSV，图默认显示且不写盘；`SAVE_RESULTS`/`SAVE_RESULTS <- TRUE`/`true`才保存CSV。

- R/RStudio：打开`配对率差.Rproj`，运行`代码/lesson28.R`。需安装`ratesci`和`ExactCIdiff`。
- Python/Spyder：打开`代码/lesson28.py`，按`# %%`分节运行。需numpy、matplotlib；从项目根或代码目录均可定位数据，非交互运行也不依赖`__file__`。
- MATLAB：打开`代码/lesson28.m`，按`%%`分节运行。只需基础函数；Current Folder可为项目根或代码目录。

## 核对与改参

基准Wald标准误为`0.1140175425`，95% Wald区间`[-.4234703,.02347028]`。R基准四法为：MOVER-NW `[-.4068100,.0428730]`、AS(Tango) `[-.4369083,.04912859]`、Wang `[-.44923,.06303]`，点估计均为`-.20`。R的`x=c(a=8,b=1,c=5,d=6)`按后测行×前测列重排；文章列表采用前测行×后测列，代码注释说明了方向。练习：将问题2的`confidence=.95`改成`.90`后重跑；R会同步重算全部四法，Python/MATLAB会显示当前Wald变窄，而问题3仍保留90/95/99%固定比较标签。

完整四法只在R入口完成；Python/MATLAB是Wald及置信水平变化的复算，不表示其他生态实现。

## 依赖与保存

R首次准备可运行`install.packages(c("ratesci", "ExactCIdiff"))`；Python在Spyder实际解释器内安装requirements.txt。脚本不自动安装依赖。Rproj采用UTF-8，不恢复或保存旧工作空间。

保存开关默认关闭。开启后，R写入`运行结果/R/四种配对率差区间.csv`（4行）；Python和MATLAB分别写入各语言运行结果目录的`Wald区间.csv`（固定90/95/99%共3行）。图默认在运行环境查看，可另存。

## 本篇下载与现有结果

代码目录：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/paired_proportion_exact_ci

本篇ZIP：https://raw.githubusercontent.com/fuyuming/statistics-classroom/main/stats_classroom/paired_proportion_exact_ci/paired_proportion_ci_project.zip

运行结果目录附默认20人、95%设置下已验证的三语言CSV。R表含四种方法，Python/MATLAB表含90/95/99%的Wald比较；改参而未开启保存时，以当前Console与图为准。

## 对应文章

[公众号源稿](文章/配对率差的精确区间_公众号稿.md) · [手机预览HTML](文章/配对率差的精确区间_手机预览.html)
