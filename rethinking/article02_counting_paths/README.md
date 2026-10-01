# 精读 02｜数路径：贝叶斯推断为什么可以先不背公式

《Statistical Rethinking》精读系列的配套代码。对应公众号「明哥的微生物世界 · Statistical Rethinking 精读 02」。

- 课程：Statistical Rethinking 2023 第 02 讲 *Garden of Forking Data*
- 教材：McElreath, R. (2020). *Statistical Rethinking: A Bayesian Course with Examples in R and Stan*, 2nd ed. 第 2 章（并用到第 3 章少量内容）

## 这一篇做了什么

1. **四面地球仪数路径**：观测序列 W L W，五个候选水面比例的路径数是 0、3、8、9、0，除以总数 20 就是后验（0、0.15、0.40、0.45、0）。
2. **20 点网格 + 三种先验**：平坦先验、`exp(-5·|p-0.5|)`、"p<0.5 时为 0"，看后验被拉动多少。
3. **解析解对照**：Beta(7, 4) 的均值 0.636364 与 20 点网格的 0.636382 相差 0.000018。
4. **从后验到预测**：beta-二项给出"未来两次取点里有几次是水"的分布（0.151515 / 0.424242 / 0.424242）。
5. **测量误差**：判定会看错时，观测到"水"有 3×2 + 1×1 = 7 条路（观测到"陆"5 条）。

## 运行

三个语言的脚本各自独立可跑，输出应当完全一致（本篇连 CSV 文件本身都逐字节一致）。

```bash
# R（主）
Rscript 代码/篇02_数路径.R

# Python
pip install -r requirements.txt
python 代码/篇02_数路径.py

# MATLAB（-batch 传不了中文路径，脚本里自带 cd）
matlab -batch "run('代码/篇02_数路径.m')"
```

输出：

- `运行结果/01–05*.csv` —— 路径计数、三种先验汇总、解析对照、预测分布、误分类路径计数
- `文章配图/*.png` —— 正文五张配图（R 版）；`06-`、`07-` 是 Python、MATLAB 的同内容对照图

## 练习开关（正文第七节）

- **练习 1**：把脚本里的 `n_faces` 从 `4` 改成 `6`（六面地球仪），路径数变成 5、16、27、32、25（合计 105）。
- **练习 2**：把 `judge_correct / judge_wrong` 从 `2 / 1` 改成 `9 / 1`（判定有 10 种取法），观测到"水"的路径数变成 3×9 + 1×1 = 28。

## 环境

- R 4.6.1：ggplot2、ragg、patchwork、dplyr／tidyr／tibble（R 代码按 tidyverse 语法写）
- Python 3.9+：numpy、scipy、matplotlib
- MATLAB R2025a：Statistics and Machine Learning Toolbox（只用 `binopdf`／`betainv`）
- 中文字体：macOS 用 PingFang SC，Windows 换 Microsoft YaHei，Linux 换 Noto Sans CJK SC

## 正文

- `数路径_公众号稿.md` —— 本篇公众号正文《9 次里 6 次是水，水面比例就是三分之二吗？——从数路径理解贝叶斯》
