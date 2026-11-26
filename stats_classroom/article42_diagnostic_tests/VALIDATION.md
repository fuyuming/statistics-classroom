# 三语项目验证记录

状态：成功。R、Python、MATLAB 完整脚本与分节入口均已运行。制作过程日志由作者保留。

## 实际环境

|环境|版本|
|---|---|
|Python|3.12.12|
|NumPy|2.4.3|
|pandas|3.0.1|
|SciPy|1.17.1|
|matplotlib|3.10.8|
|statsmodels|0.14.6|
|R|4.6.1|
|MATLAB|R2025a|

|检验|统计量|R P|Python P|MATLAB P|
|---|---:|---:|---:|---:|
|Shapiro–Wilk|0.9823344752|0.9605572776|0.9605572782|0.9605572776|
|Anderson–Darling|0.1449997844|0.9614828348|0.9614828348|0.9734209193|
|Bartlett|5.23369015899|0.0730329128|0.0730329128|0.0730329128|
|Levene（中位数）|1.27983810952|0.3134564230|0.3134564230|0.3134564230|

独立比较三语 `diagnostics_tests.csv` 的四项统计量/P 值及 `studentized_qq.csv` 全部坐标。SW、Bartlett、Levene 差均 <1e-7；AD 统计量相同，MATLAB P 值为原生 `adtest` 校准，故单列。Q-Q 坐标差均 <1e-7。

残差图算法与外部学生化残差定义已核验；`nbreaks=3/20`、`draw=false` 时三语面积为 1 且不新建图。MATLAB Shapiro–Wilk 为自编实现，许可见 `代码/LICENSE_shapiro_wilk.txt`；AD 使用 Statistics and Machine Learning Toolbox 原生函数。

ZIP解压后，R、Python、MATLAB的两个完整入口分别重新运行成功；新生成的30份CSV与当前项目结果逐列核对一致。
