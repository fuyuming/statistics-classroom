# 精读06｜控制了谁，我们就在比较谁？

对应“明哥的微生物世界” Statistical Rethinking 2023 第06讲《Good & Bad Controls》。以原课讲义与讲解为主，原书第二版第6章为补充。

- [阅读全文](好控制与坏控制_公众号稿_v7.md)
- [配套脚本、输入数据与运行说明](配套/README.md)
- [下载四组教学算例](teaching_examples.zip)

## 当前文章的算例

| 算例 | 关键结果 |
|---|---|
| 自愿报名培训，按入学成绩分层并加权 | 原始差28分；按高低档各50%重算后差10分。调整差异不自动等于培训因果效应 |
| 筛选朋友 | 两人的运动爱好独立时条件比例均50%；只留下朋友后变为80%与20% |
| 浇水与偏倚放大 | 未调整斜率约0.077克/10毫升，调整档位后为2克/10毫升；设定真实水量效应为0，两种回归均有混杂 |
| 同一回归表的系数 | 培训差10分，固定培训的经验差2分，按设定计算经验总差6分 |

## 运行

在`配套/`目录执行：

```sh
python3 teaching_examples.py
Rscript teaching_examples_base.R
Rscript teaching_examples_tidyverse.R
```

MATLAB运行`配套/teaching_examples.m`。四份脚本共享`input_data.csv`，输出保存在`配套/results/`。Python只需标准库，R基础版不需额外包，tidyverse版需tidyverse，MATLAB版使用基本运算。

这些是完全假设数据的确定性枚举、加权平均与普通最小二乘演示，不是原课quap贝叶斯拟合或重复模拟的全量复现。四语言结果此前已核验一致；本次同步复跑Python并比对已保存的其他语言结果。

## 资料

- 官方课程目录：https://github.com/rmcelreath/stat_rethinking_2023
- 原讲义：https://speakerdeck.com/rmcelreath/statistical-rethinking-2023-lecture-06
- McElreath (2020), Statistical Rethinking, 2nd edition, Chapter 6。

当前目录只保留与正文一致的四组教学算例。仓库不包含公众号密钥、草稿回执或原书全文。
