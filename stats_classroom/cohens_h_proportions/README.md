# 统计课堂30：Cohen h——两个率的角变换

统计课堂30学生项目。共享数据 data/lesson10_independent.csv 为新方案 50/100、常规对照 30/100，比较方向固定为“新方案−对照”，目标是理解 Cohen h 的角变换、方差稳定化及分级，并以原始 RD 作参照。本篇练习围绕h定义、角变换、估计精度及Cohen参照展开。

## 入口

- R：RStudio打开[项目](统计课堂30_独立两率.Rproj)，再打开[R脚本](R/independent_two_proportions.R)，从开头依次分节执行，再按需Source运行整文件。
- Python：Spyder先打开本项目文件夹，再打开[Python脚本](Python/independent_two_proportions.py)；从开头按`# %%`分节运行，F5运行整文件。项目根或`Python/`均可启动，不依赖`__file__`。
- MATLAB：Current Folder设为项目根，打开[MATLAB脚本](MATLAB/independent_two_proportions.m)；先运行初始化小节，再依`%%`逐节执行，Run可运行整文件。

数据入口：[lesson10_independent.csv](data/lesson10_independent.csv)；数值结果：[R](results/R/)、[Python](results/Python/)、[MATLAB](results/MATLAB/)。

脚本按 1→5 顺序执行：读表并算 RD=.2 与 h；hSE、双侧95% CI及 effectsize::cohens_h(t(tab),ci=.95)；固定 RD=.2 的 h 表与 |h|=.2/.5/.8 参照；n=20/100/1000 的 h CI、无连续性校正 Pearson 双侧 P；p=.1/.3/.5、n=100 的原 SE、g 斜率、近似 SE/方差及 epsilon=.01 的 exact/linear 增量。默认只打印。需要保存时，把脚本中的`SAVE_OUTPUTS`赋值行改为R的`TRUE`、Python的`True`或MATLAB的`true`；也可在运行前设环境变量`SAVE_OUTPUTS=true`。输出分别写入[results/R](results/R/)、[results/Python](results/Python/)、[results/MATLAB](results/MATLAB/)。

当前基准为新方案50/100、常规对照30/100，比较方向固定为“新方案−对照”：h=0.411516846067488，hSE=sqrt(1/n1+1/n0)，95% CI=[0.134336081198,0.688697610937]。验证状态见[VALIDATION.md](VALIDATION.md)与[validation.json](validation.json)；历史证据留在工作区，不纳入本学生包。

依赖：R的effectsize（本次1.0.3）；Python的numpy、pandas、scipy（见requirements.txt）；MATLAB使用基础函数，无统计工具箱依赖。h的角变换按弧度计算。当前四类CSV均有列名，可先与results中的参考输出核对。完整资料包另附文章科学图的Python图源。


## 配套文章

- [公众号正文](独立两率的效应大小_公众号稿.md)
- [手机预览HTML](独立两率的效应大小_手机预览.html)（与文章配图目录一起下载后打开）
- [本篇三语练习项目ZIP](三语练习项目.zip)
- [配图源代码](figure_source/绘制文章配图.py)
