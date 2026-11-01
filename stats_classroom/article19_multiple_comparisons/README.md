# 统计课堂19：多重比较与Bonferroni、Holm

数据取自课堂/教材三种疗法例题的汇总频数，未获得患者个体原始记录。`数据/therapies.csv`包含group（组别）、effective（有效人数）、ineffective（无效人数）；physical为物理疗法，drug为药物，topical为膏药。三种语言读取同一文件。

## 怎样运行

完整解压配套，保留目录层级，任选一种语言：

- R：打开本目录Rproj，再打开`代码/lesson19_multiple_comparisons.R`，Source完整运行。只用基础R。
- Python：在Spyder以本目录创建现有文件夹项目，打开`代码/lesson19_multiple_comparisons.py`，F5运行。先在Spyder实际解释器中安装requirements.txt所列依赖。主脚本自动调用同目录`plot_multiple_comparisons.py`生成两图及坐标，不能只复制主文件。
- MATLAB：将Current Folder设为`代码`目录，运行`lesson19_multiple_comparisons.m`。它按自身位置寻找数据，使用基础数组运算及gammainc，无需统计工具箱。Current Folder决定当前目录，与设置函数搜索路径不同。

结果分别保存在`运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`。从读数据开始，按章节运行；不要从中间复制一行而跳过前面的变量定义。

## 结果怎样读

整体Pearson检验不作Yates校正。每张2×2子表重新计算期望频数，再作Yates连续性校正得到原始P；此处“原始”仅指尚未作多重比较调整。不同患者独立是研究设计前提。

- `overall_result.csv`：整体χ²≈21.037651、df=2、P≈0.00002702。
- `pair_results.csv`：原顺序物理—药物、物理—膏药、药物—膏药；原始P约0.016754307、0.000009344616、0.047763774；另列校正P与子表最小期望数。
- `holm_intermediate.csv`：排序后原始P、乘数3/2/1、乘积及累计最大值。累计最大值就是本行及之前各行乘积中的最大值，随后结果被放回原来的比较位置。
- `holm_demo.csv`：人工算法例子.010/.011/.012，乘积.030/.022/.012，累计最大值均.030；不是疗法数据。
- `alpha001_exercise.csv`：把整族α改为.01，两种方法都只保留物理—膏药这一对。校正P本身不随α改变，改变的是判断门槛。

0.05标准下，Bonferroni只有物理—膏药显著，Holm三对显著。第一对Bonferroni P=.0502629197，不能先舍入成.05再判断。每个脚本还置换输入P顺序并断言还原后结果不变。

图1为各原假设全真、检验独立且单次误报率恰为5%的理论示意，不是本例三次相关检验的真实FWER。图2比较相同原始P经两种方法的结果，两图均保留逐点坐标。

本篇完整项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article19_multiple_comparisons

完整仓库：https://github.com/fuyuming/statistics-classroom

仓库ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

解压后进入`statistics-classroom-main/stats_classroom/article19_multiple_comparisons/`，先读本README。

## 补充：效应量与置信区间

三语言末节新增effect_ci.csv：物理—药物率差6.492052个百分点，近似单项95%区间0.981438至12.002665个百分点。R使用prop.test(correct=TRUE)；Python/MATLAB按其两样本区间算法复现。该区间没有作多重比较调整，不能用是否跨0代替Bonferroni的整族判断。
