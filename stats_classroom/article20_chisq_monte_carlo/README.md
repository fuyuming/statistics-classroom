# 统计课堂20：R怎样模拟卡方检验的P值？

本项目沿用课堂ACE基因型例题，比较DN与无DN两组的基因型构成。数据为249名患者的汇总频数；没有患者身份或个体记录。主线是生成一张模拟表、重复模拟、数极端结果、核对R内置结果。

## 下载后怎样运行

完整解压本项目，保留“代码、数据、运行结果、文章配图”的目录层级，任选一种语言运行：

- **RStudio**：打开本目录的`.Rproj`文件，打开`代码/lesson20_chisq_monte_carlo.R`，点击Source完整运行。脚本只使用基础R。公众号末尾额外示范的`effectsize::cramers_v()`需安装effectsize包；项目脚本用基础R实现同一非中心卡方区间算法。
- **Spyder**：以本目录创建“现有文件夹项目”，打开`代码/lesson20_chisq_monte_carlo.py`，按F5运行。依赖为numpy、pandas、scipy、matplotlib，安装时使用Spyder实际解释器。主脚本自动调用同目录`plot_chisq_monte_carlo.py`，不要只复制入口文件。
- **MATLAB**：在Current Folder进入本项目的`代码`目录，打开`lesson20_chisq_monte_carlo.m`后Run。脚本按自身位置寻找数据，使用基础MATLAB，无需Statistics工具箱。Current Folder是当前目录，和“添加搜索路径”是不同操作。

每份脚本都按读题、单次模拟、批量模拟、穷举核对、效应区间、精度练习分节。初次建议从头完整运行，再按节查看中间结果。

## 输入及六份结果文件

`数据/ace_genotype.csv`的group列有DN和noDN，分别表示有、无糖尿病肾病；DD、ID、II列为三种基因型人数。两行是42/48/21和30/72/36。行合计111/138，列合计72/120/57，总人数249。

结果保存在`运行结果/R`、`运行结果/Python`、`运行结果/MATLAB`，每种语言均有六份CSV：

- `cells.csv`：六格的组别、基因型、观察数、期望数和χ²贡献，按DN三格、noDN三格排列。
- `one_table.csv`：先抽取的一张模拟表。核对两组人数和三种基因型总数均未改变。
- `simulation.csv`：重新设置种子后，整批模拟的逐次Pearson χ²，列名Q。
- `enumeration.csv`：4063张合法条件表的第一行、χ²、出现概率与是否属于右尾。第二行可由列合计减去第一行恢复。
- `summary.csv`：观察χ²、渐近P、B、极端次数b、模拟P与SE、精确条件尾概率、Cramér's V及近似95%区间。builtinP仅在R中填写，供手动过程与chisq.test核对。
- `precision_exercise.csv`：同一模拟序列前2000、10000、100000次的b、P与模拟SE。

## 预期结果与解释

观察χ²=7.912694，自由度2，渐近P=0.01913288。穷举4063张表的概率和为1，按Pearson χ²累加得到条件尾概率约0.01928477。

种子123、B=100000时，本机R得到b=1949，P=0.019499805；重设种子后与chisq.test内置结果一致。Python和MATLAB分别约0.018949811和0.018819812，它们使用不同随机算法，均在该条件尾概率的合理Monte Carlo误差范围内。种子相同不意味着跨语言逐次相同。

V=0.1782635，近似双侧95%区间[0.0245108,0.2936630]；区间用非中心χ²反演，未作偏倚校正。模拟SE反映P值计算的随机波动，V区间反映研究效应估计的不确定性。

## 图与复现

Python入口运行后，三图和对应数据CSV保存在`文章配图/`，图标题注明所读语言。它默认绘制本次Python的结果。交付包中三图沿用文章的R模拟序列，因此与正文的1949次极端结果一致。图1对照观察表与一张模拟表，蓝色行列合计不变、橙色格内人数随抽签改变；图2上方保留完整范围，下方放大χ²约7.91–16的右尾区域；图3显示累计P与穷举条件概率。

如需重新生成与文章一致的图，先运行R入口，再从Python中调用：

```python
from pathlib import Path
import runpy
project = Path.cwd()  # Spyder项目根目录
runpy.run_path(str(project / "代码/plot_chisq_monte_carlo.py"),
              init_globals={"SOURCE_LANGUAGE": "R"})
```

完整数据序列保存在simulation.csv；绘图模块另存直方图坐标与累计P坐标，图中的右尾频数总和有断言核验。

## 自我练习

改成另一个种子，再观察2000次和100000次模拟的P值波动。两处“单次演示/批量模拟”的种子各有用途；R内置核对部分须与批量模拟使用同一种子。也可修改B，精度表会保留不超过B的检查点。B影响计算精度，总样本量仍为249。

本篇完整项目：https://github.com/fuyuming/statistics-classroom/tree/main/stats_classroom/article20_chisq_monte_carlo

完整仓库：https://github.com/fuyuming/statistics-classroom

仓库ZIP：https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

完整解压仓库ZIP，进入`statistics-classroom-main/stats_classroom/article20_chisq_monte_carlo/`，先读本README。
