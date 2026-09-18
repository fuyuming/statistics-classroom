# 运行与公开内容检查

验证日期：2026-09-18（本地）；macOS 15.6.1。

## 环境

- R 4.6.1（ggplot2、dplyr、flextable、gt、openxlsx、showtext、svglite、webshot2、eoffice 等已安装）
- Python 3.12.12（numpy、pandas、scipy、matplotlib）
- MATLAB R2025a（Statistics and Machine Learning Toolbox）

## 运行结果

- R：10 个脚本逐一 `Rscript` 运行通过；Python：10 个入口脚本逐一运行通过；MATLAB：2-1、2-3、2-9、2-10、2-11、2-12、2-13 抽查运行通过。
- 三套语言关键数字一致（逐项对照见 README 的"核对统计结果"表）：例如 MIC 几何均数 0.1015、几何标准差 3.703；身高 CV 8.0968%、体重 CV 11.2828%。
- 标准差统一使用 n−1 分母（R `sd()`、Python `std(ddof=1)`、MATLAB `std(x,0)`）；四分位数用位置线性插值。
- 图形导出：R 经 `show_fig()` 直接写入 `output/R/*.pdf`（早期用 `dev.copy(pdf, …)` 在命令行下会生成 0 页空文件，已修正）；Python/MATLAB 分别写入 `output/Python`、`output/MATLAB`。

## 公开内容检查

- 输入数据均为教学用汇总/模拟数据：80 条问卷（诊室、性别、态度）、138 个红细胞值、12 条性别记录、19 条教育等级编码、10 名新生儿身高体重、10 个 MIC 值；**不含姓名、学号、联系方式等身份字段**。
- 公开版把课堂用的 `01-.xlsx`、`02.sav`、`06.xlsx`、`example2_3.xlsx` 导出为 CSV（`data/sex_2025.csv`、`data/education_2025.csv`、`data/newborn_2025.csv`、`data/rbc.csv`），数值逐值比对一致；脚本中相应读取方式同步改为 CSV。
- 课件页码引用已从公开脚本中移除；课件文件、学生数据与个人路径未纳入仓库。
- `output/` 为运行产物，未纳入版本记录（见 `.gitignore`）。

未在 Windows 上实机安装软件；不同系统的编码与字体问题请先阅读 README 的排查说明。
