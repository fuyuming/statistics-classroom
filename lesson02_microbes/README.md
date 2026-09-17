# 统计课堂 · Python与MATLAB快速入门：微生物数据练习

本包对应第二篇公众号。第一篇的专业构成与抛硬币例子保留在仓库根目录；本篇请进入 **lesson02_microbes** 再运行，不要混用两课的数据路径。

## 数据是什么

`data/microbes_od600_demo.csv` 是人工构造的模拟教学数据，不来自真实实验，不含个人信息。

设同一非致病微生物在培养基A、B中各有8个独立培养物，在同一固定终点记录扣除空白后的OD600；每组人为设置1个缺失。每一行代表一个独立培养物，不是重复读数，也不是生长曲线的一个时间点。

| 字段 | 含义 |
|---|---|
| culture_id | 虚构培养物编号A01—A08、B01—B08 |
| group | medium_A、medium_B |
| od600 | 模拟OD600，空白代表缺失 |

OD600是无量纲浊度读数，不能直接解释为活菌数或生长速率；模拟数据也不能证明培养基的实际效果。缺失不填零；课堂计算先报告缺失，再使用有效观测。

## 下载与运行

在仓库首页点击 Code → Download ZIP，完整解压后进入lesson02_microbes。直接下载链接：

https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

Python：打开 `02_python_quickstart.py`，在Spyder中完整运行。需要pandas、numpy、matplotlib；也可在本文件夹的终端输入 `python 02_python_quickstart.py`。脚本使用自身路径定位数据，完整运行时才有 `__file__`。

MATLAB：打开 `matlab_quickstart.m` 并点击Run。仅使用基础MATLAB，文件名须以字母开头，文件采用UTF-8。已在R2025a执行。

R对照：打开 `course_demo.Rproj`，运行 `source("01_R_quickstart.R", encoding="UTF-8")`。此脚本仅需基础R。

每种语言输出到对应的 `output/Python`、`output/MATLAB`、`output/R`：summary.csv、group_plot.png、histogram.png及版本记录。重新运行会覆盖对应练习结果，不修改输入数据。

## 核对统计结果

| group | 总记录 | 有效数 | 缺失 | mean | sd | median | q1 | q3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| medium_A | 8 | 7 | 1 | 0.414286 | 0.028200 | 0.410 | 0.395 | 0.430 |
| medium_B | 8 | 7 | 1 | 0.514286 | 0.028200 | 0.510 | 0.495 | 0.530 |

SD使用n−1分母。四分位数用位置1+(n−1)p的线性插值，三种语言口径一致。两组教学数值特意保留相同的离散程度，方便区分平均水平与离散程度。散点图横线表示样本均值，没有显著性检验或误差区间含义。

直方图采用相同的分箱边界；小样本与边界归箱会影响形状，不能单凭形状判断分布。结果图使用英文标签配合文章中文图注，避免依赖某一操作系统的字体。

## 编程基础与课后练习

`R_basics.R`、`python_basics.py`、`matlab_basics.m` 是独立的小型语法示例，不是另一份微生物实验数据。

1. 先解释数据类型、结构、索引、条件、循环和函数。
2. 选一种语言运行完整OD600流程，解释每组有效数与缺失数。
3. 比较两组平均水平和离散程度，并说明为什么不能由模拟结果宣称某培养基更好。
4. 复制CSV修改一个值，再改脚本输入文件名，观察均值、中位数与SD的变化。
5. 重启软件，从头运行，确认不依赖旧会话变量。

## 中文乱码排查

脚本注释乱码：先按正确编码重新打开，确认后另存UTF-8，不要覆盖乱码文本。

数据乱码：本CSV使用UTF-8 BOM；Python用utf-8-sig，R用UTF-8-BOM，MATLAB使用readtable的Encoding='UTF-8'。其他来源的CSV应核对其实际编码，不要忽略解码错误。

图上中文方框：检查绘图字体。Python用matplotlib.font_manager查询已安装字体；MATLAB用listfonts。选择本机真实存在的中文字体，再重新绘图、导出。修改字体不能修复已经读错的文字。

## 软件入口

- Spyder：https://www.spyder-ide.org/download
- Spyder安装：https://docs.spyder-ide.org/current/installation.html
- MATLAB试用：https://www.mathworks.com/campaigns/products/trials.html
- R：https://cran.r-project.org/

验证：三种语言的微生物统计表已核对，全部数值列误差小于1e-10；输出图已检查。测试环境为macOS，未逐一测试Windows软件安装。
