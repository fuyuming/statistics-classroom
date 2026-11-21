# 统计课堂36：P值调整（主练习）

本项目的默认入口只有 `01_p_adjust`：读取共享 CSV，练习 Bonferroni、Holm、BH 和 BY 四种 P 值调整，并复算第19篇的三疗法例题。R、Python、MATLAB 三选一即可。

## 先打开项目，再运行主练习

共享输入已随项目提供，三种语言都直接读取 `data/screening_pvalues.csv` 和 `data/therapies.csv`，无需先运行 R 生成数据。

| 语言 | 工作目录与分节运行 | 主入口 |
| --- | --- | --- |
| R | 用 RStudio 打开本目录的 `article36_p_adjust_qvalue.Rproj`；从第一节用 Ctrl/Cmd+Enter 逐节运行，或用 Source 整文件运行 | `scripts/01_p_adjust.R` |
| Python | 在 Spyder 中打开本目录，并将本目录设为工作目录；从首个单元用 Shift+Enter 逐单元运行，或用 F5 整文件运行 | `scripts/01_p_adjust.py` |
| MATLAB | 将 Current Folder 设为本目录，执行 `addpath('scripts')`；在 Editor 按 `%%` 顺序运行，或在命令窗口运行 `m01_p_adjust` | `scripts/m01_p_adjust.m` |

MATLAB 的 Current Folder 决定相对数据路径，搜索路径负责找到脚本和函数。保持 Current Folder 为项目根目录；若 Editor 询问处理脚本位置，选择加入路径，不要切换到 `scripts` 目录。脚本会自行加入随附的 `helpers` 目录。

命令行备选（均在项目根目录运行）：

```text
Rscript scripts/01_p_adjust.R
python scripts/01_p_adjust.py
```

Python 使用 `requirements.txt` 中列出的常规依赖。MATLAB 核心计算使用随附的 `p_adjust_helpers.m`，无额外工具箱要求；若已安装 Bioinformatics Toolbox，01 中的 `mafdr(...,'BHFDR',true)` 仅用于可选 BH 核对。

## 输入与统计对象

`data/screening_pvalues.csv` 含 1000 个模拟 P 值：800 项代表没有差异，200 项代表有差异。它们是 P 值层面的教学输入，不是患者原始记录；调整算法只读取 `p` 列，不读取模拟时附带的 `is_null` 标签。

- `id`：指标编号；正文的转录组筛选情境可将其对应为基因编号。
- `p`：算法读取的 P 值。
- `is_null`：模拟时已知的真相标签，仅用于说明输入，不参与调整。

生成来源为 `data/generate_pvalues.R`，R 种子为 20261121。三种语言读取同一 CSV，便于逐项比较算法输出。基因例有 1000 项（`m=1000`）；01 另有两组各自单独构成一个三项 Holm 比较族（均为 `m=3`）：正文停止例 `0.020、0.023、0.040`，以及累计值递增例 `0.020、0.040、0.100`。两组不合并为一个比较族。

`data/therapies.csv` 是第19篇的课堂汇总计数：物理疗法199人有效、7人无效，药物治疗164/18，外用膏药118/26。三次两两比较分别取2×2表，采用Yates连续性校正，再将三个P值一起作Bonferroni或Holm调整。这里的 `m=3` 数的是三次两两比较。病例来源见[统计课堂19](https://mp.weixin.qq.com/s/Ox3tXlD1XchvU6GSmtiuWg)。

## 怎样读主练习

01 先展开排序、位置、候选乘数、累计处理和还原顺序，再用完整函数核对。Holm 取前向累计最大，BH 取后向累计最小，BY 在 BH 基础上乘以调和数。原始 P 值与调整值均与目标水平 0.05 比较。

固定输入的发现数为：原始/Bonferroni/Holm/BH/BY = **149/23/23/67/38**。Holm停止例的调整值为 **0.060/0.060/0.060**；递增例的调整值为 **0.060/0.080/0.100**。这些是当前输入的发现数，不等同于重复研究定义的 FWER/FDR。疗法原题在0.05水平下，Bonferroni只判定物理与膏药有差异；Holm判定三对均有差异。

## 再改一次水平，观察名单怎样变化

跑通后，试着把01.1中的`alpha`从0.05改为0.01，保留同一份1000项P值，重新从头运行。记录Bonferroni、Holm、BH和BY的发现数，比较名单怎样变化。疗法原题仍按它原先设定的0.05水平核对。

## 可选保存与结果文件

默认只在会话中查看对象，运行不会改写 `results`。需要保存时，先设置相应语言的开关，再运行主入口：

```text
R:       Sys.setenv(ARTICLE36_SAVE_RESULTS="1")
Python:
import os
os.environ["ARTICLE36_SAVE_RESULTS"] = "1"
MATLAB:  setenv('ARTICLE36_SAVE_RESULTS','1')
```

主练习 01 的结果 CSV 使用各语言的 `*_01_p_adjust.csv`、`*_01_summary.csv`、`*_01_mini_example.csv` 和 `*_01_therapies.csv` 文件名；不设置开关则不写文件。可将开关设为 `0` 恢复仅在会话查看。

## 方法来源

- R 调整方法：https://search.r-project.org/R/refmans/stats/help/p.adjust.html
- Python 调整方法：https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
