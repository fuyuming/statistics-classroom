# 统计课堂43｜方差分析的正态性、方差齐性：指什么，一定要检验吗？

从三组独立试件的完整例题出发，说明正态性、共同方差怎样支撑经典F检验，以及怎样结合已有依据、残差图与必要的正式检验判断模型条件。

- [阅读全文](公众号稿.md)
- [下载本篇R、Python与MATLAB配套项目](统计课堂43_ANOVA条件与残差_配套项目.zip)：解压后进入“配套项目”。
- [三语言运行说明](配套项目/README.md)
- [实际运行与三语言结果核对](配套项目/VALIDATION.md)
- [离线手机预览](手机预览.html)：下载后用浏览器打开，图片已嵌入。

## 运行练习

RStudio打开 `配套项目/ANOVA条件与残差.Rproj`，运行 `代码/anova_residuals.R`。
Python在Spyder将工作目录设为“配套项目”，按 `# %%` 运行 `代码/anova_residuals.py`。
MATLAB将Current Folder设为“配套项目”，按 `%%` 运行 `代码/anova_residuals.m`。

三种语言读取同一份 `数据/residual_cases.csv`，仅包含主例场景A：三组、每组5个独立观测。数据为教学模拟资料。依赖、图形和结果开关见配套README。

## 重建图与预览

安装 `构建工具/requirements.txt` 所列依赖，并准备支持中文的字体。在本目录执行：

```sh
python3 构建工具/绘图.py
Rscript 构建工具/绘制残差分布.R
python3 构建工具/构建预览.py
```

正文是唯一内容源，配图和预览由构建工具生成。项目保留按组Shapiro结果、Anderson–Darling未运行状态、中位数中心Levene及Bartlett结果，不把P＞0.05解释为模型前提已被证明。
