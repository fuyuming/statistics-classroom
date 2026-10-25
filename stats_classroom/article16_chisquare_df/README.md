# 统计课堂16：四格表为什么自由度是1

本篇项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article16_chisquare_df

完整下载：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

下载并完整解压，进入 `statistics-classroom-main/stats_classroom/article16_chisquare_df/`，选择熟悉的一种语言即可。脚本按四部分展开，从数据准备开始，按顺序选中运行；命名对象用于查看结果，保存结果另设可选节。

- R：双击本目录的Rproj，打开 `R/001_四格表自由度.R`，可以Source整文件或按注释分段运行。只用R基础包。
- Python：在Spyder中把工作目录设为本文件夹，打开 `Python/001_four_cell_df.py`，F5运行或从第一节按 `# %%` 顺序运行。需NumPy、SciPy，依赖见requirements.txt。控制台中不依赖 `__file__`。
- MATLAB：把Current Folder设为本文件夹，打开 `MATLAB/four_cell_df.m`，运行整份或从第一节依次运行；需要Statistics and Machine Learning Toolbox。

数据沿用课堂例6：行是女性、男性，列是收看、未收看，计数为225/275、128/272。`data/`保存同一教学表；脚本中直接写出四格人数便于新生逐行核对。

先核对E、df=1、未校正X²=15.75423及Yates校正X²=15.21361，再把x=225改为226并按合计约束改另外三格。输出保存在 `output/`，实际运行证据与版本说明见 `VALIDATION/README.md`。文中的卡方输出截图原始文本另保存在 `R运行输出/`。

率差为13个百分点；`ci_corrected`是本题与R prop.test默认校正同口径的近似95%区间[0.06457346,0.19542654]，不要与未校正区间混用。
