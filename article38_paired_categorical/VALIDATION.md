# 本篇计算核验（两例均8人）

例1输入保持不变：`000` 1人、`001` 6人、`111` 1人；C=(1,1,7)、T=9、D=12，Cochran Q=12，渐近P=e^-6，Q条件汇总28配置、逐人安排729种、精确P=1/243。例1 `pairwise` 仍为 comparison、b、c、rate_diff、p_raw、holm。

例2新增 `data/transition_records.csv`，8行按A/B/C顺序动态汇总为 M=[[1,2,0],[0,1,2],[1,0,1]]；`transition_counts.csv`为9格对照。行和=(3,3,2)、列和=(2,3,3)、d=(1,0,-1)，S前2=[[3,-2],[-2,4]]，S⁻¹=[[0.5,0.25],[0.25,0.375]]，Stuart–Maxwell=0.5，df=2，χ²近似P=0.778800783071405；Bowker=5，df=3，χ²近似P=0.171797144296733。两者均仅作小样本渐近示例。

例2 `categories` 的(leave,enter)为(2,1)、(2,2)、(1,2)，率差为−0.125、0、0.125，精确P及Holm均为1；`paths`正反为(2,0)、(0,1)、(2,0)，theta为1、0、1，原始P为0.5、1、0.5，Holm全1，CP区间分别为[0.158113883008,1]、[0,0.975]、[0.158113883008,1]。Bowker条件枚举路径总数2/1/2，共18配置，概率和1，右尾P=0.25。M2为另一组8人矩阵[[1,2,0],[0,1,2],[2,0,0]]，SM=0、Bowker=6。

三语主入口均直接读取例2逐人记录；七张结果表保留。新复制、无旧results、不同cwd新进程运行的退出码、源码/输入SHA、七表跨语言比较及例1不变核对见 `work/small_example2/project_acceptance.json`。

## 本轮R条件枚举变更（2026-11-23）

R脚本将Bowker三层循环替换为自编`bowker_exact3(M)`：仅接受3×3且三对反向格合计均>0，按固定合计的二项条件分布对观测统计量作右尾精确枚举，不使用蒙特卡罗；`enumeration`字段与数值保持不变（18、1、0.25）。README同步说明了可选ANSM5核对及300000上限（默认100000；本例门槛6^7=279936）。

本轮仅在新副本`work/exact_r_addition/worker_run`中执行`ARTICLE38_SAVE=1 /usr/local/bin/Rscript scripts/article38_paired_categorical.R`，退出码0。七张R结果表与当前`outputs/.../results/R`基线按数值误差≤1e-10比较全部通过（各表最大绝对误差0）；运行日志、源码/输入SHA、比较JSON及退出码记录在`work/exact_r_addition/`。ANSM5数值为既有实跑记录，非本轮新增执行。
