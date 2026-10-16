# 配套教学算例 v2

本目录提供四个可复核的、完全假设数据的课堂算例。四份脚本共同读取 `input_data.csv`，并把长表结果写到 `results/`。脚本只做确定性枚举、加权平均和普通最小二乘（OLS），不生成图形。

这些数字用于说明控制、标准化和条件化的算术关系。它们与原课中的 `quap`、贝叶斯后验以及重复模拟不同；算例 3 的四行数据也只用于有限行的系数演示，不输出显著性或标准误（SE）。

## 算例与应看到的数字

1. 培训标准化：原始培训/未培训均值为 `84/56`，差 `28`；把高、低分档各设为 50% 后为 `75/65`，差 `10`。
2. 朋友筛选：四格（都运动、仅甲、仅乙、都不运动）筛选前为 `25,25,25,25`，筛选后为 `20,5,5,20`。按正文方向计算 `P(乙运动 | 甲运动)` 与 `P(乙运动 | 甲不运动)`，结果由 `50%/50%` 变成 `80%/20%`。
3. 偏倚放大：`Z=[0,0,1,1]`、`U=[0,1,0,1]`、`X=[90,110,190,210]`、`Y=[10,14,10,14]`。`Y~X` 的 X 系数为 `1/130 = 0.0076923...`，`Y~X+Z` 的 X 系数为 `0.2`。这里把真效应设为零，Y 的差异由 U 造成。
4. Table 2：经验 1 年未训/已训均值 `50/60`，经验 2 年 `52/62`，培训率 `.2/.6`；两年的总均值为 `52/58`，差 `6`；固定培训状态下经验差均为 `2`，固定经验下培训差均为 `10`。

## 运行

在本目录执行：

```bash
python3 teaching_examples.py
Rscript teaching_examples_base.R
Rscript teaching_examples_tidyverse.R
```

MATLAB（命令行或 MATLAB 界面）运行 `teaching_examples.m`。脚本用自身所在目录定位 CSV，因此从其他工作目录调用也可以。四个脚本会分别生成：

```text
results/python_results.csv
results/r_base_results.csv
results/r_tidyverse_results.csv
results/matlab_results.csv
```

Python 版只用标准库；R 基础版只用 R 自带的 `stats::lm`；R tidyverse 版需要已经安装 `tidyverse`，脚本不会自动安装依赖；MATLAB 版只用 `readtable`、反斜杠 OLS 与 `writetable`。若当前环境没有 tidyverse 或 MATLAB，请做静态检查，确认路径、字段名、公式和预期数值，不要把缺少环境当成数字验证通过。

## 输入字段

`input_data.csv` 是四个算例的共享输入。培训算例使用 `stratum/trained/mean_value/count`；朋友筛选使用 `phase/friend_cell/count`；偏倚放大使用 `x/y/z/u`；Table 2 使用 `experience_years/trained/mean_value/training_rate`。空字段是其他算例的占位，不代表缺失数据需要插补。
