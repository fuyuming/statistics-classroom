% 统计课堂35｜关联强度与区间。输入为第11课的假设工程试验汇总计数。
% 300台独立记录器随机分3组各100台：A硅胶、B环氧、C热缩封装。
% 弯折100次后接模拟器记录10min，无持续≥1s数据缺口判合格（课堂阈值）。
% 主表按A/B/C行、合格/不合格列排列为[[90,10],[80,20],[65,35]]。
% 在项目根目录按节顺序运行；默认不写文件，ARTICLE35_SAVE_RESULTS=1才保存。

%% 02.1 准备输入
% 完整run可能暂时进入scripts；仅在找不到输入时，用本脚本位置回到项目根目录。
if ~isfile(fullfile('data','device_pass.csv'))
    task_script_path = mfilename('fullpath');
    task_project_root = fileparts(fileparts(task_script_path));
    cd(task_project_root);
end
addpath('helpers');
raw = readtable('data/device_pass.csv');
disp(raw)
tab = table2array(raw(:, {'pass','fail'}));
assert(all(sum(tab,1)>0) && all(sum(tab,2)>0));

%% 02.2 调用函数计算V及双侧95%区间
% 这是本项目配套自编函数，内部使用MATLAB Statistics Toolbox数值接口。
fit = cramers_v_ci(tab,.95);

%% 02.3 读取结果
% V为整体关联强度，V_L/V_U为区间下限/上限，confidence记录置信水平。
V = fit.V;
V_L = fit.V_L;
V_U = fit.V_U;
ci = table(V, V_L, V_U, .95, ...
    'VariableNames', {'V','CI_low','CI_high','confidence'});
disp(ci)

%% 02.4 效应分级
interpretation = table(["lower";"estimate";"upper"], [V_L;V;V_U], ...
    interpret_v_funder([V_L;V;V_U]), repmat("funder2019",3,1), ...
    'VariableNames', {'component','V','label','rule'});
assert(all(interpretation.V >= 0 & interpretation.V <= 1));
disp(interpretation)

%% 02.5 本例解释与报告
% 本段用fprintf组织报告句。
fprintf('本例Cramér V=%.6f，双侧95%%区间=[%.6f, %.6f]；区间表示估计不确定性。\n', ...
    V, V_L, V_U);
% 第1、2、3行分别是下限、点估计、上限的等级。
labels = interpretation.label;
fprintf('按Funder–Ozer规则，点估计为%s；区间两端等级为%s/%s。\n', ...
    labels(2), labels(1), labels(3));

%% 02.6 可选保存
if strcmp(getenv('ARTICLE35_SAVE_RESULTS'),'1')
    if ~isfolder('results')
        mkdir('results');
    end
    writetable(ci,'results/MATLAB_ci.csv');
    writetable(interpretation,'results/MATLAB_interpretation.csv');
end
