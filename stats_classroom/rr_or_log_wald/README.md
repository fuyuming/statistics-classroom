# 统计学堂32：RR、OR与log-Wald区间

这是给初学者的三语配套项目。数据是人为设定的两组独立康复教学数据，不用于真实医学结论。推荐先打开项目目录，再按脚本分节从0到5运行；也可整文件运行。学生可任选一种语言。

## 运行入口

- R：在RStudio打开 `统计学堂32_RR_OR.Rproj`，确认Files面板位于项目根目录；从0节开始逐段运行（选中代码后Run），或Source整文件。需 `DescTools`；首次安装：`install.packages("DescTools")`。
- Python：在Spyder把项目根目录设为工作目录，确认解释器是已安装 `requirements.txt` 的环境；打开 `代码/rr_or_log_wald.py`，按 `# %%` 的顺序逐节运行或F5整文件。命令行可运行 `python 代码/rr_or_log_wald.py`；依赖见 `requirements.txt`，可在该解释器中执行 `python -m pip install -r requirements.txt`。
- MATLAB：在MATLAB把项目目录设为 Current Folder，打开 `代码/rr_or_log_wald.m`，按 `%%` 分节或Run整文件；脚本用基础 `erfcinv` 计算正态分位点，不需要工具箱。

脚本默认关闭 `SAVE_RESULTS`，只在控制台显示结果，不写结果文件。要导出相同schema的 `results.csv`，设置环境变量 `RR_OR_SAVE_RESULTS=1`（MATLAB 可在命令行执行 `setenv('RR_OR_SAVE_RESULTS','1')`），文件分别写入 `运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`；也可在脚本开头直接修改开关。

保存开关在每份脚本开头：R写`SAVE_RESULTS <- TRUE`，Python写`SAVE_RESULTS = True`，MATLAB写`SAVE_RESULTS = true`。首次学习保持关闭即可。MATLAB的搜索路径决定去哪里查找函数，Current Folder决定相对路径从哪里出发；本项目把项目根目录设为Current Folder，不需要另加搜索路径。

Python在Spyder中先运行文件顶部的导入与目录设置，再从第0节往下运行；可用“运行当前单元”逐节执行。R同样先运行顶部设置，再逐段Run。分节运行依赖前面已经算出的变量，因此请按顺序进行。

## 学习顺序

0 读取两组计数并检查正整数；1 计算RD、RR、odds与OR并解释对象；2 从二项方差和log/odds导数计算两个log尺度SE；3 在log尺度构造95%区间并exp回原尺度；4 交换两组，观察比值倒数、SE不变及端点 `[1/U,1/L]`；5 将两组人数同时乘4，观察点估计不变、SE减半。最后请自行解释：为什么RD不变而log-Wald区间变窄？

第2节中的log/odds导数具体指log p和logit(p)的导数；logit(p)=log[p/(1−p)]。

## 核对结果与练习

| 对象 | 本例预期结果 |
| --- | --- |
| RD | 0.20，即高20个百分点 |
| RR、OR | 1.66666667、2.33333333；分别是概率比与优势比 |
| 两个log尺度SE | 0.1825741858、0.2960051480 |
| RR的95% log-Wald区间 | [1.16530959，2.38372515] |
| OR的95% log-Wald区间 | [1.30622126，4.16808747] |
| 交换组序 | RR变为0.60、OR变为3/7，SE不变，区间为原端点按[1/U,1/L]变换 |
| 人数同时乘4 | 两个SE分别变为0.0912870929、0.1480025740，均减半 |

人数乘4的练习表示在相同概率条件下观察更多相互独立的患者；把已有每条患者记录复制四遍不会增加独立信息。两个样本率保持不变，点估计就保持不变；样本率方差中的分母变为4n，标准误因开平方而缩小为原来的1/2。再把脚本第5节的倍数改为9，先预测SE会怎样改变，再运行核对。

输入文件只接受new/control各唯一一行，success/failure是有限正整数，无缺失；本项目同时演示未经调整的RR与OR log-Wald，因此遇到零计数会在计算前给出说明，不自动加0.5。改变输入时保持每组患者数与教学问题一致。

R脚本另外演示 `DescTools::RelRisk`、`OddsRatio`、默认 `delta=.5` 与 `delta=0` 的差别，以及 `fisher.test` 的条件OR区间；Fisher结果不是对Wald区间的验证。Python/MATLAB只使用本语言公式，不调用R。

本篇完整项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/rr_or_log_wald

## 文章与配图

[文章源稿](RR与OR_公众号稿.md) · [手机预览HTML](手机预览.html) · [下载本篇完整项目ZIP](RR_OR_完整项目.zip) · [只下载三语言练习](RR_OR_三语言练习.zip)

完整项目ZIP解压后进入`rr_or_log_wald/`，先读本README；仅练习ZIP解压后进入`配套项目/`。`构建源码/build_assets.py`生成18幅公式和2幅科学示意图，`构建源码/build_preview.py`从正文生成预览。封面沿用第10课无姓名课堂封面。
