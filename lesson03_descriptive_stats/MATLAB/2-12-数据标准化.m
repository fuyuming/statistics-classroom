% 对应：2026 年秋医学统计学·统计描述（第二次课，2026-09-18）；课件为课程内部资料，未随仓库公开
% 原R脚本：2-12-数据标准化.R；共用函数 medstats_examples.m 按同一数据逐步实现。
% 用MATLAB打开本文件点击运行，或run('本文件完整路径')。
% 依赖Statistics and Machine Learning Toolbox（boxplot及形状指标所需函数）。
% 共用函数medstats_examples.m按同一原数据逐步实现，含详细统计注释。
% 本脚本全部为确定性计算（样本SD分母n-1），可与R、Python逐位对照。
script_dir=fileparts(mfilename('fullpath'));
addpath(script_dir);
result=medstats_examples('2-12');
