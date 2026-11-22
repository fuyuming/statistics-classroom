# 本篇计算核验（当前8人版本）

例1逐人输入仅8行：`000` 1人、`001` 6人、`111` 1人；`transition_counts.csv`（例2）保持原样。三语主入口均动态读取各自项目根目录的 `data/`，输出七张表：`key_results`、`q_conditional`、`pairwise`、`categories`、`paths`、`enumeration`、`exercise`。

例1核对值：C=(1,1,7)，T=9，D=12；Cochran Q=12，df=2，渐近P=e^-6=0.002478752176666。固定行和的汇总枚举为28个配置，逐人等概率安排数为3^6=729，Q≥12的条件精确P=1/243=0.00411522633744856。

例1 `pairwise` 仅含 comparison、b、c、rate_diff、p_raw、holm：AB为b=c=0、P=1；AC和BC均为b=0、c=6、P=0.03125，Holm分别为1、0.09375、0.09375；率差依次为0、0.75、0.75，分母为8。例2 `categories` 保留率差及Wald区间；`paths` 保留互转者方向概率及CP区间。

本轮核验要求：R、Python、MATLAB从新复制项目、无旧results、不同工作目录新进程运行，退出码均为0；七张表跨语言逐字段比较，数值绝对误差≤1e-10。例2四张结果表（`categories`、`paths`、`enumeration`、`exercise`）与 `work/raw_records_run/acceptance_project/results` 逐字段数值未变。源码与输入均记录SHA-256；日志和机器可读证据见 `work/eight_person/project_acceptance.json`。
