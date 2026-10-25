# 验证记录

日期：2026-10-24。项目主入口为 `配套练习/`，三语言均从此工作目录运行。

- R 4.6.1：`Rscript --vanilla R/001_四格表自由度.R`，通过；输出：`R运行输出/001_main.txt`，列联表完整核对：`R运行输出/列联表输出.txt`。
- Python 3.12.11（Spyder运行环境），NumPy 2.5.2、SciPy 1.18.0；脚本通过；结果：`output/Python/001_result.json`。按 `# %%` 顺序、无 `__file__` 的分段执行通过；结果：`VALIDATION/python_sections.json`。
- MATLAB R2025a 25.1.0.2943329：`-batch` 实跑通过；输出：`output/MATLAB/001_run.txt`。正式入口已改名为 `MATLAB/four_cell_df.m`，在本项目中文路径下直接运行通过；实测日志为 `output/MATLAB/正式文件运行.txt`，数值CSV为 `output/MATLAB/001_result.csv`。
- 三语言均包含 `ci_corrected`：本题近似公式与 R `prop.test` 默认校正一致，端点约 0.06457346、0.19542654；不是精确区间。
- 未运行项：无。
