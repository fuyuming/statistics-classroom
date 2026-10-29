# 统计课堂17：泊松分布拟合优度

这是第11课 CLABSI 计数例题的三语言配套练习。数据用于课堂演示，不标记为可验证的真实研究数据。公开项目目录：[GitHub 项目目录](https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article17_poisson_gof)；[仓库 ZIP](https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip)。

## 下载与运行

下载 ZIP 后完整解压；本项目根目录应为 `statistics-classroom-main/stats_classroom/article17_poisson_gof/`，其中直接包含 `代码/`、`数据/` 和 `运行结果/`。

- R：双击 `统计课堂17-泊松分布拟合优度.Rproj`，打开 `代码/lesson17_poisson.R`，Source 整个脚本或从第一部分顺序运行。只需 base 与 stats。
- Python/Spyder：将本目录设为工作目录，打开 `代码/lesson17_poisson.py`，按 F5 或依次运行四个 `# %%` 单元。脚本不依赖 `__file__`，每段使用前段建立的命名对象。依赖见 `requirements.txt`。
- MATLAB：把 Current Folder 设为本目录，打开 `代码/lesson17_poisson.m`，Run 或依次 Run Section。需要 Statistics and Machine Learning Toolbox 中的 `poisspdf`、`poisscdf`、`chi2cdf`；`fminbnd` 为基础 MATLAB 函数。

若提示找不到数据，请确认当前工作目录中直接有 `数据/` 和 `代码/`。三套脚本默认把结果写入各自的 `运行结果/<语言>/` 子目录。

## 数据与流程

共享 CSV 为 `数据/poisson_weekly_counts.csv`：一周计数 0、1、2、3、4 对应周数 223、142、48、15、4；总周数 432、总事件数 299。原始均数 299/432 只作引入和初查，不与分组 MLE 混用。检验保留 0、1、2、≥3 四组，尾组理论概率包含未观测的 ≥5。

四部分依次为：看数据与假设；选四组并查条件；按同一分组数据估计 λ、计算期望数和 Pearson 贡献；解释结果并完成参数练习。正式拟合为分组 MLE λ≈0.69130516，df=2。练习另行使用预先给定 λ=0.7，重新计算 E/Q，按事先给定参数使用 df=3；比较时不要只替换自由度。

每种语言的 `运行结果/` 保留 4 个实际数值 CSV：主结果、分组期望与贡献、参数比较、λ=0.7 的期望数。
