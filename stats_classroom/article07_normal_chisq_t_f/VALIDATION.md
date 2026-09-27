# 运行验证

已实际运行 R 4.6.1、Python 3.12.12（NumPy 2.4.3、SciPy 1.17.1、Matplotlib 3.10.8）、MATLAB R2025a（Statistics and Machine Learning Toolbox 25.1）。

三套脚本内部恒等式检查均通过。全部理论计算表、三张图的坐标、同一教学数据的 t/ANOVA/回归结果，按 rtol=1e-9、atol=1e-10 交叉核对通过。

```text
critical_values.csv / Python: max absolute difference 2.27e-13
critical_values.csv / MATLAB: max absolute difference 9.95e-13
t_normal_curves.csv / Python: max absolute difference 1.39e-15
t_normal_curves.csv / MATLAB: max absolute difference 2e-15
chisq_curves.csv / Python: max absolute difference 1.95e-14
chisq_curves.csv / MATLAB: max absolute difference 9.99e-16
f_curves.csv / Python: max absolute difference 8.1e-15
f_curves.csv / MATLAB: max absolute difference 3e-15
t_anova_regression_equivalence.csv / Python: max absolute difference 5.33e-14
t_anova_regression_equivalence.csv / MATLAB: max absolute difference 9.97e-18
simple_regression_equivalence.csv / Python: max absolute difference 2.7e-13
simple_regression_equivalence.csv / MATLAB: max absolute difference 0
```

各语言模拟流不同，模拟数据不逐行对齐；各自检查 T²=F 并输出经验分位点。项目路径由脚本定位。UTF-8；无中文绘图字体依赖。未测试其他操作系统的安装流程。
