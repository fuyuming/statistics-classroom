# 运行核对记录

本记录只写实际运行过的环境与关键输出，不把文件存在当作计算通过。

## 实跑环境

- R：R 4.6.1（2026-06-24）；使用 `Rscript 代码/t_tests_as_regression.R` 运行成功。
- Python：Python 3.12.12；使用 `MPLBACKEND=Agg python 代码/t_tests_as_regression.py` 运行成功。自动核对时设置了 `MPLBACKEND=Agg`，脚本在 Spyder 中仍保留 `plt.show()`。
- MATLAB：MATLAB R2025a；使用 `matlab -batch` 运行 `代码/t_tests_as_regression.m` 成功。需要 Statistics and Machine Learning Toolbox。

## 跨语言确定性结果

三种语言的 `运行结果/<语言>/all_key_results.csv` 相互一致（末位浮点舍入除外）：

| example | estimate | df | t | P | 95% CI |
| --- | ---: | ---: | ---: | ---: | ---: |
| single_B_vs_100 | 1.000000 | 6 | 0.1282875 | 0.9021134 | [−18.0736590, 20.0736590] |
| two_independent_groups_pooled_OLS / `group01` | 19.166667 | 17 | 1.9156630 | 0.0723832 | [−1.9425426, 40.2758759] |
| paired_new_minus_old | 1.300000 | 9 | 2.7508479 | 0.0224457 | [0.2309462, 2.3690538] |

两独立样本的 Welch 对照为 `t=1.9319308`、`df=13.0158821`、`P=0.0754315`、95%区间 `[−2.2636714, 40.5970047]`。普通 OLS 的整体 `F=3.6697649`、`df1=1`、`df2=17`、`P=0.0723832`，等于组差系数的 `t²`；它不是总体方差比 F 检验。

R、Python、MATLAB 均生成了两张 PNG：组均数与 OLS 截距、新旧款配对连线。图片分别保存在项目 `文章配图/`，并在各语言的 `运行结果/` 中留有一份。

## 解释边界

- B组与100 g的区间是“均数−100 g”的区间；换回均数尺度为 `[81.9263, 120.0737] g`。
- A−B 的区间跨过0，P>0.05只表示这组数据在该模型下没有拒绝零差异，不能说两组相同或已经等效。
- 配对分析使用同一人的 `new-old`，不能把 old、new 打散后当成两组独立观测。
