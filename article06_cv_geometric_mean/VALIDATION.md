# 运行与公开内容检查

验证日期：2026-09-18；macOS 15.6.1。

## 环境

- R 4.6.1（`DescTools`、`EnvStats`、`moments`、`nortest`、`showtext` 等已安装；`2-9` 的 R 脚本会调用多个统计包）
- Python 3.12.12（numpy、pandas、scipy、matplotlib）
- MATLAB R2025a（Statistics and Machine Learning Toolbox）

## 运行结果

- R：`2-9-几何均数.R`、`2-10-变异系数.R` 运行通过。
- Python：`2-9-几何均数.py`、`2-10-变异系数.py` 运行通过（共用模块 `medstats_examples.py`）。
- MATLAB：`medstats_examples('2-9')`、`medstats_examples('2-10')` 运行通过。
- 三套语言关键数字一致：MIC 算术均数 0.225、中位数 0.0625、几何均数 0.10153155、几何标准差 3.70265996；身高 CV 8.0968196%、体重 CV 11.2828279%。
- 标准差统一使用 n−1 分母（R `sd()`、Python `std(ddof=1)`、MATLAB `std(x,0)`）；四分位数按位置线性插值（R 默认 type=7，MATLAB 以 `q7` 复现同一规则）。
- 图形导出：R 写入 `output/R/`（`show_fig()` 直接输出 PDF；早期 `dev.copy(pdf, …)` 在命令行下会生成 0 页空文件，已修正），Python 写入 `output/Python/`，MATLAB 写入 `output/MATLAB/`。

## 公开内容检查

- 数据为教学/模拟数据：10 个 MIC 值（二倍稀释档位）与 10 组新生儿身高体重；**不含姓名、学号、联系方式等身份字段**。
- 公开版把课堂使用的 `06.xlsx` 导出为 `data/newborn_2025.csv`（数值逐值比对一致）；MIC 数据沿用课堂的 `05.txt`。
- 脚本仅保留这两节（几何均数、变异系数）；课件页码引用已移除，未纳入课件文件、学生数据与个人路径。
- `output/` 为运行产物，未纳入版本记录（见 `.gitignore`）。

未在 Windows 上实机安装软件；不同系统的编码与字体问题请先阅读 README 的排查说明。
