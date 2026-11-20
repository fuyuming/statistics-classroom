# 统计课堂33：独立与配对

数据是用于教学的假设例题，不是真实临床记录。项目根目录包含 `data/`、`scripts/`、`results/`；相对路径可移植到RStudio、Spyder或MATLAB。

## 六个分析脚本

`01_paired_means_R.R`、`01_paired_means.py`、`m01_paired_means.m`：8名不同患者，每人一段固定30秒心电记录，A/B在同一记录估计心率（次/分），比较B−A均值。`02_paired_binary_R.R`、`02_paired_binary.py`、`m02_paired_binary.m`：100名患者，A/B分别标记一段心电中是否至少出现一次室性早搏；无金标准，不比较诊断准确性。MATLAB脚本加`m`前缀是因为MATLAB文件名不能以数字开头。

二分类表实际读取 `data/ecg_2x2.csv`：行是A阳性/阴性，列是B阳性/阴性；本例无缺失，每人贡献一个格子。连续数据读取 `data/ecg_rates.csv`，8个ID一一对应且无缺失。

## 运行方式

RStudio打开 `article33_independent_paired.Rproj`，再打开`scripts/01_paired_means_R.R`，从“问题与输入”节运行到最后；之后同样运行`02_paired_binary_R.R`。选中整节后用Ctrl+Enter（Mac为Cmd+Enter）执行，也可Source整份文件。

Spyder以本项目作为工作目录，打开`scripts/`中的对应Python文件，按`# %%`从第一单元逐节执行，也可按F5运行整文件。Python只按当前项目目录读取数据，不依赖`__file__`。命令行备选为`python scripts/01_paired_means.py`。

MATLAB把Current Folder设为本项目根目录，先在命令窗口执行`addpath('scripts')`，再打开`m01_paired_means.m`或`m02_paired_binary.m`，从第一节开始按`%%`逐节运行。若提示切换目录或加入路径，选择加入路径，保持Current Folder在项目根目录。整文件备选是在命令窗口输入脚本名，如`m01_paired_means`。

R使用基础包；Python分析需要`numpy`、`pandas`、`scipy`和`statsmodels`，绘图另需`matplotlib`、`seaborn`，依赖列在requirements.txt；MATLAB需要Statistics and Machine Learning Toolbox。Python可在自己使用的解释器环境安装`pip install -r requirements.txt`。运行会将每种语言的两张结果表写入results/，互不覆盖。

每个脚本都先读入数据并显示输入，再分别显示手算对象与函数核对结果。连续例先看n=8、平均差−2、差值SD约0.9258、t约−6.1101、df=7、双侧P约0.0004862和95%区间［−2.7740，−1.2260］次/分；二分类例先看n=100、m=30、率差0.20、校正Q约12.0333、P约0.0005226、精确P约0.0003249和普通配对Wald区间约［0.09956，0.30044］。率差区间乘100才是百分点；它并非校正检验的反演区间。

可练习的改动：改变A/B一对数值并观察均值差、SD、SE和CI；交换二分类表的b/c并观察率差方向；把ID、A、B一起整行重排，结果应不变；只倒序B会改变差值波动。均值差不说明算法准确性，阳性判定率更高也不说明性能更好。

`scripts/plot_pairs.py`可从项目根目录运行，重建8对心率连线图到`results/paired_plot.png`。文中概率独立的定义是联合分布性质；配对数据结构保留同对内联系，简单配对推断仍要求不同对子满足相应独立性前提。
