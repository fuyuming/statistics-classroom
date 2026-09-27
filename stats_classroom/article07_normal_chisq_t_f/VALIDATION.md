# v3 运行验证

重写后的三语言主脚本均已完整运行：R 4.6.1、Python 3.12.12（NumPy 2.4.3、pandas 3.0.1、SciPy 1.17.1、Matplotlib 3.10.8）、MATLAB R2025a（Statistics and Machine Learning Toolbox 25.1）。

理论分位点、错误率、三张曲线全部坐标以及共用教学数据的检验／回归结果，按 rtol=1e-9、atol=1e-10 交叉核对通过。各脚本的恒等式检查通过。

```text
critical_values.csv / Python: max absolute difference 2.27e-13
critical_values.csv / MATLAB: max absolute difference 9.95e-13
error_rates.csv / Python: max absolute difference 1.11e-16
error_rates.csv / MATLAB: max absolute difference 1.6e-15
t_normal_curves.csv / Python: max absolute difference 1.39e-15
t_normal_curves.csv / MATLAB: max absolute difference 2e-15
chisq_curves.csv / Python: max absolute difference 1.95e-14
chisq_curves.csv / MATLAB: max absolute difference 9.99e-16
f_curves.csv / Python: max absolute difference 8.1e-15
f_curves.csv / MATLAB: max absolute difference 3e-15
t_anova_regression_equivalence.csv / Python: max absolute difference 5.33e-14
t_anova_regression_equivalence.csv / MATLAB: max absolute difference 2.98e-13
simple_regression_equivalence.csv / Python: max absolute difference 4.41e-13
simple_regression_equivalence.csv / MATLAB: max absolute difference 1.2e-11
```

五张图在三种语言中均已生成，正文使用R版本。新增构造图和错误率图已目视检查。模拟使用独立标准正态流，跨语言不要求逐行相同；保存全部抽样值与经验分位点。主脚本采用中文逐题注释，绘图排版单独放在plot_distributions配套文件中。

项目包含Rproj、README、requirements，脚本UTF-8。未测试其他操作系统的软件安装流程。
