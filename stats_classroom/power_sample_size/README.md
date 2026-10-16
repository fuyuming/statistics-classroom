# 统计功效与样本量：鱼苗例

完整解压后选择一种语言运行即可。三种语言计算同一套数值；文章的四幅图及公式图片由Python生成，R/MATLAB入口只计算并导出结果表。

- R：打开“统计功效.Rproj”，打开代码/power_design.R，点击Source。仅需基础R。
- Python：在Spyder中打开代码/power_design.py，F5运行。用运行Spyder的同一Python环境安装requirements.txt列出的依赖。脚本会显示图，并保存到文章配图。
- MATLAB：打开代码/power_design.m，点击Run。需要Statistics and Machine Learning Toolbox。脚本自动定位文件；Current Folder表示当前目录，不等于永久加入搜索路径。

结果在运行结果/对应语言目录。设计参数在数据/设计参数.csv；mu0、mu1、sigma的单位为cm，alpha和target_power为概率。无缺失值，无真实个体测量数据。来源：李春喜《生物统计学》第6版例4-1数值，另拟单样本右侧设计。

预期：100尾功效0.812405；最小样本量97；96尾功效0.798322，97尾0.801928。观察结果另设双侧情境，差值固定0.40 cm。

练习：修改共享CSV的mu1为7.45，再分别运行，观察最小样本量变为386。或将target_power改为0.90。只在研究设计阶段改变假设，不以观察到的P值决定追加样本。

Python配图中的标题、图注按文章基准参数撰写；修改参数练习后，勿将新图直接用于原稿。文章配图与基准运行结果已随包保存。

plot_power.py、plot_threshold.py与render_formulas.py由Python主入口自动调用，分别生成主图、α权衡图与数学公式。

## Z与t的条件

情境表使用已知总体标准差σ=1.58的Z检验。如果1.58是样本标准差s、总体标准差未知，应改用t分布，自由度n−1，P值和区间也需重新计算；不能只更改检验名称。25、100、400尾均适用这一区分。
