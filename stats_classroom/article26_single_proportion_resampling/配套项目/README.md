# 统计课堂26配套：非参数重抽样与蒙特卡洛

本配套项目对应“非参数重抽样与蒙特卡洛”示例。R脚本是Bootstrap、固定率原假设模拟和Permutation的主要演示；Python与MATLAB保留同题的单样本率基础练习及绘图/官方接口对照。

## 运行顺序

1. 在RStudio打开 `article26_single_proportion_resampling.Rproj`。
2. 先运行 `R/infer_walkthrough.R`，查看观察率、Bootstrap区间和固定率 H0 模拟。
3. 运行 `R/article26_v8_additions.R` 查看率差区间与单组换序。
4. 运行 `R/permutation_walkthrough.R` 查看置换分布、P值和70种精确分组。
5. `R/boot_examples.R`、`R/package_examples.R` 是可选的其他接口示例；`R/single_proportion_resampling.R` 展示基础R实现。

R依赖：R 4.1+；`infer`、`ggplot2`、`boot`，以及可选的 `confintr`、`ragg`、`systemfonts`。Python依赖见 `requirements.txt`；MATLAB基础脚本需要 Statistics and Machine Learning Toolbox。三种语言均使用种子20261110，但随机数发生器不同，数值允许有小差异。

## 主要结果

- 单样本：180人中117人达标，观察率0.65；Bootstrap双侧95%百分位区间 `[0.5777778, 0.7166667]`。
- H0: p=0.45 的10000轮模拟命中为0；精确右尾二项 P=`5.23826971452729e-8`。
- 两组置换：A=3/4、B=1/4；10000轮命中2462次（0.2462），70种精确分组中17次达到观察差值。

单组换序示例输出三次0.65至 `data/infer_v4_shuffle_stats.csv`。数据和运行记录在 `data/`、`results/`；生成图在上一级 `文章配图/`。`VALIDATION.md`记录本轮环境与验证结果。
