# 统计课堂40｜趋势检验到 Logistic 回归

Logistic 回归用评分等输入描述二分类事件发生的概率；本项目的事件是检测失败。这是同一份课堂假设数据的 R、Python、MATLAB 初学者练习：200 名独立患者各有一段、每段 100 个已标注 QRS，固定算法漏检超过 5 个判失败、不超过 5 个判成功；记录随机分为 4 组，每段仅接受一种干扰强度，各患者记录独立，每组 50 人，50 Hz 干扰峰峰值为 0/.1/.2/.3 mV，预评分为 0:3，失败数为 5/10/18/27。数据在 `数据/共享教学数据.csv`，无缺失。

| 字段 | 含义与单位 |
| --- | --- |
| `score` | 预定评分 0–3；脚本默认每单位对应 0.1 mV |
| `amplitude_mV` | 50 Hz 干扰峰峰值，mV |
| `failed` | 漏检超过 5/100 个 QRS 的记录数 |
| `n` | 该组总记录数 |

## 从解压到第一次运行

1. 完整解压本目录，不要单独拖出脚本。
2. 选择一种语言：RStudio 打开 `.Rproj`；Python 在 Spyder 创建或打开本文件夹项目，将工作目录设为项目根，再打开 `.py`；MATLAB 打开 `代码/trend_logistic.m`，推荐把项目根目录设为 Current Folder（代码目录也兼容）。MATLAB 的 Current Folder 决定相对路径起点，搜索路径只用于查找函数，不能替代数据所在目录。
3. 在新会话按 `01.1`→`01.2`→`01.3`→`01.4`→`01.5` 顺序运行，或整文件运行：先准备题干数据与单位，再先看 Logistic 模型结果和 OR，再看预测图，最后用趋势统计量核对 Wald/Rao/LRT 与 logLik。每节使用前节创建的对象，脚本不要求每节独立运行；若重新开始，请从 `01.1` 起。

R 使用 `stats`/`graphics` 和 `lmtest`（在 R 控制台运行 `install.packages("lmtest")`）；Python 在 Spyder 所用解释器的终端运行 `python -m pip install -r requirements.txt`；MATLAB 需要 Statistics and Machine Learning Toolbox（`glmfit`、`glmval`）。终端备选（项目根目录）：`Rscript 代码/趋势检验_Logistic.R`、`python 代码/趋势检验_Logistic.py`、`matlab -batch "run('代码/trend_logistic.m')"`。Rscript 无图形窗口时可能创建默认 `Rplots.pdf`；默认不保存仅指 IDE 练习模式。

## 每节在回答什么

`01.1` 展示完整假设题并读取检查字段、准备失败/成功和单位参数；`01.2` 拟合二项 Logistic，先展示 coef、beta、SE、OR 与双侧 95% Wald CI；`01.3` 显示四组观察率和拟合曲线；`01.4` 计算 U、I_eff、趋势 Q/P，再比较得分/Rao、Wald、LRT、logLik（R 同时展示 `lmtest::lrtest` 接口）；`01.5` 汇总原有结果并提供动态报告/可选保存。这里把 I_eff 作为得分 U 的方差估计使用。R 的 Rao 为原生比较，Python/MATLAB 按 `U²/Ieff` 手算同一得分 Q，均未调用 R 后端。

## 结果显示、参考结果与保存

默认只在 IDE 显示数据、解释和图，不覆盖文件。`运行结果/` 中现有 CSV/PNG 是随包提供的参考输出，来自当前输入和同一计算流程在 2026-11-22 新手顺序重排后的默认整文件实际运行；这些随包文件是参考结果。同版脚本另已通过顺序分节和单位换算练习。新输入或改参数后，必须重新顺序运行；只有打开保存开关才会更新文件。

输出固定在项目根的 `运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`，不随 Current Folder 改变。保存前设置 R `save_outputs <- TRUE`、Python `save_outputs = True`、MATLAB `save_outputs = true`，或设环境变量 `ARTICLE40_SAVE_OUTPUTS=1`；关闭则用 FALSE/False/false。参考核对值：U=37，Ieff=52.5，Q=26.0761904762，P≈3.282064e-7；beta≈0.7832025，OR≈2.188470，LRT Q≈27.443894。

结果汇总CSV的 `mV_per_score` 列记录当次模型中一个评分单位对应多少mV，OR及其区间均按这一单位解释。完整预期字段与结果见 [`VALIDATION.md`](VALIDATION.md)。

## 第二遍练习

把某组 `failed` 改为相邻整数，观察 U、拟合率和三种检验变化。完成后恢复原始CSV，再做单位练习：在所选脚本 `01.1` 的唯一教学参数处把 `mV_per_score` 从 `0.1` 改为 `1`，然后从 `01.1` 按顺序运行；不改 CSV，也不扩大观察范围。此时模型中的 `score` 按 mV 表示，预测、趋势 Q/P 和换算回每 0.1 mV 的 OR 应保持不变，输出的“每一个评分单位 OR”会按新单位显示。

文章中的儿童肥胖真实研究用于说明研究目标和资料结构。本项目的计算仍只复现表1的假设心电数据；原研究数据、评分与多因素模型不在本项目中重算。
