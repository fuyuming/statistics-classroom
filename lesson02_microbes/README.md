# 统计课堂 · Python与MATLAB快速入门：微生物数据练习

本包对应第二篇公众号。第一篇的专业构成与抛硬币例子保留在仓库根目录；本篇请进入 **lesson02_microbes** 再运行，不要混用两课的数据路径。

## 数据是什么

`data/microbes_od600_demo.csv` 是人工构造的模拟教学数据，不来自真实实验，不含个人信息。

设同一非致病微生物在培养基A、B中各有8个独立培养物，在同一固定终点记录扣除空白后的OD600；每组人为设置1个缺失。每一行代表一个独立培养物，不是重复读数，也不是生长曲线的一个时间点。

| 字段 | 含义 |
|---|---|
| culture_id | 虚构培养物编号A01—A08、B01—B08 |
| group | medium_A、medium_B |
| od600 | 模拟OD600，空白代表缺失 |

OD600是无量纲浊度读数，不能直接解释为活菌数或生长速率；模拟数据也不能证明培养基的实际效果。缺失不填零；课堂计算先报告缺失，再使用有效观测。

## 下载与运行

在仓库首页点击 Code → Download ZIP，完整解压后进入lesson02_microbes。直接下载链接：

https://github.com/fuyuming/statistics-classroom/archive/refs/heads/main.zip

Python：打开 `02_python_quickstart.py`，在Spyder中完整运行。需要pandas、numpy、matplotlib；也可在本文件夹的终端输入 `python 02_python_quickstart.py`。脚本使用自身路径定位数据，完整运行时才有 `__file__`。

MATLAB：打开 `matlab_quickstart.m` 并点击Run。仅使用基础MATLAB，文件名须以字母开头，文件采用UTF-8。已在R2025a执行。

R对照：打开 `course_demo.Rproj`，运行 `source("01_R_quickstart.R", encoding="UTF-8")`。此脚本仅需基础R。

每种语言输出到对应的 `output/Python`、`output/MATLAB`、`output/R`：summary.csv、group_plot.png、histogram.png及版本记录。重新运行会覆盖对应练习结果，不修改输入数据。

## 核对统计结果

| group | 总记录 | 有效数 | 缺失 | mean | sd | median | q1 | q3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| medium_A | 8 | 7 | 1 | 0.414286 | 0.028200 | 0.410 | 0.395 | 0.430 |
| medium_B | 8 | 7 | 1 | 0.514286 | 0.028200 | 0.510 | 0.495 | 0.530 |

SD使用n−1分母。四分位数用位置1+(n−1)p的线性插值，三种语言口径一致。两组教学数值特意保留相同的离散程度，方便区分平均水平与离散程度。散点图横线表示样本均值，没有显著性检验或误差区间含义。

直方图采用相同的分箱边界；小样本与边界归箱会影响形状，不能单凭形状判断分布。结果图使用英文标签配合文章中文图注，避免依赖某一操作系统的字体。

## 编程基础与课后练习

`R_basics.R`、`python_basics.py`、`matlab_basics.m` 是独立的小型语法示例，不是另一份微生物实验数据。

1. 先解释数据类型、结构、索引、条件、循环和函数。
2. 选一种语言运行完整OD600流程，解释每组有效数与缺失数。
3. 比较两组平均水平和离散程度，并说明为什么不能由模拟结果宣称某培养基更好。
4. 复制CSV修改一个值，再改脚本输入文件名，观察均值、中位数与SD的变化。
5. 重启软件，从头运行，确认不依赖旧会话变量。

## 中文乱码排查

脚本注释乱码：先按正确编码重新打开，确认后另存UTF-8，不要覆盖乱码文本。

数据乱码：本CSV使用UTF-8 BOM；Python用utf-8-sig，R用UTF-8-BOM，MATLAB使用readtable的Encoding='UTF-8'。其他来源的CSV应核对其实际编码，不要忽略解码错误。

图上中文方框：检查绘图字体。Python用matplotlib.font_manager查询已安装字体；MATLAB用listfonts。选择本机真实存在的中文字体，再重新绘图、导出。修改字体不能修复已经读错的文字。

## 软件入口

- Spyder：https://www.spyder-ide.org/download
- Spyder安装：https://docs.spyder-ide.org/current/installation.html
- MATLAB试用：https://www.mathworks.com/campaigns/products/trials.html
- R：https://cran.r-project.org/

验证：三种语言的微生物统计表已核对，全部数值列误差小于1e-10；输出图已检查。测试环境为macOS，未逐一测试Windows软件安装。

## 项目与工作目录设置

### 先设置项目和工作目录，再读数据

“工作目录”就是程序解析 `data/文件名.csv` 这类相对路径时的起点。**打开了脚本，不一定就已经切换到脚本所在文件夹。** 本篇请把起点设为同时包含脚本和 `data` 的 `lesson02_microbes` 文件夹，不是其中的 `data` 子文件夹。

在 Spyder 中按下面的顺序操作：

1. 完整解压 GitHub 下载包，找到 `statistics-classroom-main/lesson02_microbes`。
2. 选择菜单 **项目（Projects）→ 新建项目（New Project）**。
3. 在创建窗口左侧选择 **现有文件夹（Existing directory）**，选中刚才的 `lesson02_microbes`，然后创建项目。已有数据和脚本时不用选“新建文件夹”，也不用再套一层空目录。
4. 在项目文件树中打开 `02_python_quickstart.py`。以后通过“打开项目”或“最近的项目”返回这里。

Spyder 会使用项目根文件夹设置工作目录；运行配置也可能影响执行时的位置，所以最后在 **IPython 控制台**核对一次。[Spyder 官方项目说明](https://docs.spyder-ide.org/current/panes/projects.html)

```python
from pathlib import Path
print(Path.cwd())
print([p.name for p in Path.cwd().iterdir()])
print((Path.cwd() / "data" / "microbes_od600_demo.csv").is_file())
```

应能看到 `data` 和 `02_python_quickstart.py`，最后返回 `True`。不创建项目也能使用 Spyder：在顶部工作目录栏旁的文件夹按钮中选择同一文件夹，再用上述命令确认。

如果需要用代码切换，可以在控制台输入下面两行，运行时粘贴**自己电脑上该文件夹的完整路径，不加外层引号**：

```python
import os
os.chdir(input("请粘贴 lesson02_microbes 文件夹的完整路径：").strip())
```

本篇完整脚本通过 `__file__` 定位数据；控制台里逐句试验时则使用已经核对的 `Path.cwd()`。这两个起点可能不同，要知道自己正在用哪一个。

### MATLAB：先切换“当前文件夹”，再理解“设置路径”

顶部 **设置路径（Set Path）** 按钮打开的是函数搜索路径列表。把文件夹加进这个列表，并不等于把它设为当前工作目录；也不会替你把 `data/...` 的相对路径起点改到那里。本篇无需在该窗口添加文件夹或保存路径。

实际操作只要两步：

1. 在工具栏下方的**当前文件夹地址栏**输入或浏览到 `lesson02_microbes`。左侧 Files / Current Folder 中应该同时看到 `data` 和 `matlab_quickstart.m`。
2. 打开 `matlab_quickstart.m`，点击 Run；若提示文件不在当前文件夹，并提供 **Change Folder（更改文件夹）**，本练习选择更改到脚本所在文件夹。

也可以在命令窗口用文件夹选择框完成切换，无需手写长路径：

```matlab
folder = uigetdir(pwd, '选择含data和脚本的lesson02_microbes文件夹');
if ~isequal(folder, 0) % 点击取消时返回0，不执行切换
    cd(folder);
end
pwd
dir
isfile(fullfile(pwd, 'data', 'microbes_od600_demo.csv'))
```

最后一行应返回逻辑值 `1`（true）。`pwd` 查看当前位置，`cd` 改变当前位置，`dir` 列出其中的文件。[MATLAB 官方 cd 说明](https://www.mathworks.com/help/matlab/ref/cd.html)

**当前文件夹、脚本所在文件夹、搜索路径，是三个不同概念。** 本篇完整 `.m` 脚本根据自身位置寻找数据；上面的目录设置则方便在命令窗口逐句练习与检查文件。只有以后需要调用其他目录里的公共函数时，才进一步讨论 `addpath` 与持久保存搜索路径。

本地单独练习ZIP解压后的文件夹可能叫 `microbes_demo`；选择同时含脚本和 `data` 的那一层即可，文件夹名称本身不影响运行。
