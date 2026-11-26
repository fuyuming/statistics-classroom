# 当前运行依据

2026-11-26；R 4.6.1、Python 3.12.12、MATLAB R2025a。三语主入口均以场景A的15行数据运行，新会话按当前入口实际分节通过：R 7节、Python 8节、MATLAB 8节。

- 结果目录仅保留6类CSV：`group_summary.csv`、`residuals.csv`、`anova.csv`、`qq.csv`、`variance_tests.csv`、`normality_tests.csv`；三语字段、行序和标签一致，数值最大绝对差约1.03e-12。
- A主例：每组n=5；SSB=998.942813333333、SSE=51.10356、F(2,12)=117.28452734017、P≈1.3287891e-8。Shapiro按组实际运行；Anderson–Darling因n=5明确记录为 `not_run_n_lt_8`，统计量和P为NA。
- `SHOW_PLOTS=true` 分支三语均实际运行三张诊断图；R图例颜色与点颜色一致。
