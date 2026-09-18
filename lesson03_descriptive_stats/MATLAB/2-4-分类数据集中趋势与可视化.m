% 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
% 原R脚本：2-4-分类数据集中趋势与可视化.R；共用函数 medstats_examples.m 按同一数据逐步实现。
% 用MATLAB打开本文件点击运行，或run('本文件完整路径')。
% 依赖Statistics and Machine Learning Toolbox，在“主页→附加功能”安装并确认许可。
% 共用函数medstats_examples.m按同一原数据逐步实现，含详细统计注释。
% 随机模拟不要求跨语言逐位一致；确定性结果统一分母、样本SD和R type=7分位数。
script_dir=fileparts(mfilename('fullpath'));
addpath(script_dir);
result=medstats_examples('2-4');
