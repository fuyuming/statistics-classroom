# 精读03 v5：同一个模型，网格与 quap 两种计算

先读 `一条直线_手机预览_v5.html`，正文源文件为 `一条直线_公众号稿_v5.md`。旧版本保留，但本版数字仅对应 `运行结果/v5/`。

## 下载与打开

项目：https://github.com/fuyuming/statistics-classroom/tree/main/rethinking/article03_geocentric_line

ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

下载解压后进入 `rethinking/article03_geocentric_line`。

- RStudio：打开 `一条直线.Rproj`，打开 `代码/howell_v5.R`，点击 Source。依赖 `rethinking`，安装以作者 README 为准：https://github.com/rmcelreath/rethinking 。本脚本使用 quap，不需要运行 Stan 采样。
- Spyder：打开本课文件夹项目，打开 `代码/howell_v5.py`，F5。实际解释器需要 numpy、scipy、matplotlib（命令行也可 `python 代码/howell_v5.py`）。它自动调用同目录 `howell_plots_v5.py`；两个文件都要保留。完整网格检查需几百 MB 内存。
- MATLAB：打开本课目录，将其设为 Current Folder，打开 `代码/howell_v5.m`，Run。仅基础 MATLAB，不需要统计工具箱。

要复现全部资料，可先运行 R，再运行 Python，最后运行 MATLAB。Python 不依赖 R 才能计算；如果存在本次 R 输出，会核对后优先使用真实 quap 的均值和协方差。

## 模型与作者对应

作者脚本：https://github.com/rmcelreath/stat_rethinking_2023/blob/main/scripts/03_howell_new_weight_model.r

定位 `m_adults`，而非同文件后面的男女分组或多项式模型。使用 352 位年龄≥18岁的 Howell1 成年人；共享 CSV 已按身高排序，字段为 height（厘米）、height_c（中心化厘米）、weight（千克），来源由旧 R 数据导出步骤记录。CSV 数值保留六位小数。

给定身高和参数，观测条件独立：

- W_i ~ Normal(mu_i, sigma)，第二个参数为标准差。
- mu_i = a + b × (H_i − mean(H))。
- a ~ Normal(60,10)，b ~ LogNormal(0,1)，sigma ~ Uniform(0,10)，先验独立。

两种方法都估计三个参数。R 真正调用 `rethinking::quap`；Python/MATLAB 的独立二次近似实现不是其他语言的 quap 函数。联合近似保留协方差，没有把边际后验独立相乘。

## 计算与精度

三维网格：a∈[43,47]、b∈[0.45,0.8]、sigma∈[3.3,5.5]，各101点。这是后验积分窗口，不是改变先验。网格体积相同，归一化时共同约掉；采用对数权重避免直接连乘下溢。

Python 还运行各151点的加密网格，以及 a∈[42,48]、b∈[0.4,0.85]、sigma∈[3,6] 的扩大网格。后验均值变化小于 1e-5；外层网格点总质量约 4.66e-9。这是边缘诊断，不是对窗口外质量的严格上界。数值窗口针对本数据核验过，更换数据要重新检查。

边际等尾区间使用网格中心累计权重插值，受步长限制。预测采用262144个扰乱 Sobol 低差异点（固定种子20261005）；网格按条件分布依次取 sigma、b、a，保留依赖；quap 从联合正态变换生成。新观测额外使用独立维度的正态分位数。输出保留较多位仅为复核，正文用约数。半数与全数积分点检查端点变化，阈值0.02千克；R 的10万次普通随机模拟是独立思路复核，区间不会逐位相同。

本版两类区间均展示89%，而作者脚本平均趋势使用99%、个体预测使用89%，这是为便于比较作出的展示调整。

## 结果在哪里

- `comparison.csv`：网格与 quap（或同原理近似）的均值和89%边际区间。
- `R_quap_parameters.csv`、`R_quap_covariance.csv`：R quap 实际输出。
- `R_grid_means.csv`、`Python_grid_means.csv`、`MATLAB_comparison.csv`：三语言对照。
- `grid_posterior.npz`：三维网格坐标与全部归一化权重，用 numpy.load 打开。
- `grid_predictions.csv`、`quap_predictions.csv`：各身高的平均体重及两种89%区间。
- `R_quap_simulation_160.csv`：R后验抽样与新个体模拟结果。
- `checks.json`：网格范围、密度与预测积分检查。
- `文章配图_v5/`：6张统计图，另有公式图和横/方封面。文件名前缀沿用旧图资产编号，以图内与正文图号为准。

图4只展示均值线，不是完整先验/后验预测分布；先验允许不合理直线的问题在正文明确展示。所有统计图由脚本计算，封面和公式图沿用已有图像资产。

## 练习

第一次核对 b≈0.6287，grid sigma≈4.257，quap sigma≈4.230。第二次加密网格，观察结果变化与组合数量。更换先验或数据时，两种算法应同步修改，并重新检查计算窗口，不能继续套用旧数值。
