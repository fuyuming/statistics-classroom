# 统计课堂：Bootstrap 与置换检验

本目录是《Bootstrap与置换检验：反复抽样，究竟在模拟什么？》的配套练习。三种语言使用同一份大鼠数据和同一组固定 Bootstrap 索引，因此可以逐值核对。数据是教材例题的教学复现，不含学生个人资料。

## 从哪里开始

完整下载并解压仓库 ZIP 后，进入 `stats_classroom/bootstrap_permutation`。先打开一种熟悉的语言，完整运行该语言的主题入口：

- R：`代码/bootstrap_permutation.R`
- Python / Spyder：`代码/bootstrap_permutation.py`
- MATLAB：`代码/bootstrap_permutation.m`

三份脚本都按“看数据 → 普通组内 Bootstrap 估计 → 中心化、学生化的 H0 参照 → 置换检验 → 小样本率例”分节。第一次先核对 `运行结果/` 中的结果表，再改一个参数重新运行。

## 数据与固定索引

`data/rats.csv` 是两组大鼠增重（g）：

- A 组：134, 146, 106, 119, 124, 161, 107, 83, 113, 129, 97, 123（12只）
- B 组：70, 118, 101, 85, 107, 132, 94（7只）

数据来自李春喜：《生物统计学》第6版，例4-8，pp.97–98；本文用 A−B 表示均数差。`bootstrap_indices_A.csv` 为 10000×12、`bootstrap_indices_B.csv` 为 10000×7，每行是一轮组内有放回抽样的索引。索引从1开始；Python读入后减1，R和MATLAB原样使用。固定索引不代表增加了真实动物数，只是让三套脚本对同一批课堂演示轮次。

`data/rate_examples.csv` 用于边界教学例：9/10 与 0/10。它们不是本次大鼠数据的另一项研究。

## R / RStudio

双击 `统计课堂-Bootstrap与置换检验.Rproj`，打开 `代码/bootstrap_permutation.R`，点击 **Source**。也可在项目目录运行：

```r
source("代码/bootstrap_permutation.R", encoding = "UTF-8")
```

第0.1节直接调用 `boot.pval::boot_t_test()`，首次使用先运行 `install.packages("boot.pval")`。后面的原理演示使用基础R，绘图采用基础图形设备。脚本不恢复旧工作空间，也不依赖上次会话的变量。

## Python / Spyder

在 Spyder 中选择 **Projects → New Project → Existing directory**，选择本目录，打开 `代码/bootstrap_permutation.py` 后按 F5。需要 NumPy、pandas、SciPy、Matplotlib；在运行脚本的同一解释器中安装：

```sh
python -m pip install -r requirements.txt
python 代码/bootstrap_permutation.py
```

脚本使用 `__file__` 定位项目目录；命令行无图形界面时可用 `MPLBACKEND=Agg`。运行结果目录保留英文通用图；文章配图优先使用系统中文字体，没有中文字体时自动改画英文。

## MATLAB

需要 Statistics and Machine Learning Toolbox。把 Current Folder 切到本目录，打开 `代码/bootstrap_permutation.m` 并点击 Run，或运行：

```matlab
run('代码/bootstrap_permutation.m')
```

Current Folder 与搜索路径是两回事；本脚本按自身位置查找 `data/`，无需把整个目录永久加入 MATLAB 搜索路径。

## 运行结果和图

三语言分别写入 `运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`，避免互相覆盖；`文章配图/` 保存 Python 入口生成的两张文章图。关键结果应为：

- 观察均数差 A−B = 19.1666667 g，观察 T = 1.931930817；
- 普通 Bootstrap SE 约 9.270038 g，百分位区间 [1.047619, 37.369345] g（只作估计演示）；
- 中心化、学生化 H0 参照：左 9636、右 364、等尾 P = .0728；
- 与该 T 分布配套的 Bootstrap-t 95% CI [-2.28522435, 40.92199759] g；
- 19个观测混合后按 12/7 重分的置换检验，直接枚举 C(19,12)=50388 个分组，三语言共用同一精确置换分布，双侧 P = 0.0735294；
- 9/10 的普通百分位率区间为 [.7, 1]，0/10 退化为 [0, 0]，CP 95% 区间为 [0, .308497]。真实 p=.9 的全成功概率为 .9^10=.34867844；对指定区间算法、给定教学模型作无限 Bootstrap 分布精确枚举，覆盖率为 .6496866225（不是模拟结果）。

`文章配图/01_两种重抽样与统计量.png` 对照估计 Bootstrap 与 H0 参照；`文章配图/02_置换分布与观察差值.png` 解释置换检验如何混合19个观测后重分。两张图都由 Python 入口脚本生成；R 和 MATLAB 也会在各自结果目录生成同名图，供跨语言查看。

## 参考资料

- Hesterberg (2015), *What Teachers Should Know About the Bootstrap*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4784504/
- R `boot` 文档：https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.html
- SciPy `permutation_test` 文档：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html

## 练习

1. 在读入固定索引后取前2000行，观察较少重抽样时区间端点和尾部计数的模拟波动；不要把 B 当成样本量。
2. 只把两组均数差的定义改成 B−A，检查观察差值和区间端点如何变，双侧 P 是否改变。
3. 为什么50388个分组全枚举后不再有有限重抽样次数的随机波动，而10000轮Bootstrap仍有？比较两种P值的计算分母。
4. 将 9/10 改成 0/10，先预测普通 Bootstrap 与 CP 区间，再运行核对。

练习1请在计算重抽样轮数之前截取两份索引的前2000行，并只运行到估计与检验部分。末节的固定数值断言用于核验原始10000轮案例，修改数据或轮数后不应继续用原数值作为通过标准。练习2先在纸上判断符号与区间的变化，再逐节修改，保持所有统计量的差值方向一致。

## 中文与编码

脚本和数据使用 UTF-8。Python 绘图会查找系统中文字体（微软雅黑、黑体、Noto Sans CJK SC、苹方等）；若找不到则改用英文图，数值照常输出，避免生成方框字。安装任一上述字体后重启 Python。请保持整个项目解压后的目录结构，输出均留在本项目内。

## 专用包快捷调用

R入口第0.1节使用boot.pval 0.7.0、seed=20261023、R=10000、type="stud"、var.equal=FALSE：差19.1667g、95%CI[-2.319064,40.927598]、P=.0730。结果另存`运行结果/R/boot_pval_package_result.csv`。包调用与共享索引演示不是同一批抽样，不要求末位一致。
