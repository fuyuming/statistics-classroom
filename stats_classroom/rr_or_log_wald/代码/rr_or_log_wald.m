%% 统计学堂32：RR、OR与log-Wald区间（MATLAB，UTF-8）
% 数据为人为设定教学数据。先把项目根目录设为 Current Folder，再按各 %% 运行。

%% 0 读取与检查
% 改为true可保存；环境变量用于教师批量核对。
SAVE_RESULTS = false;
if strcmp(getenv('RR_OR_SAVE_RESULTS'), '1')
    SAVE_RESULTS = true;
end
projectDir = pwd;
dataPath = fullfile(projectDir, '数据', 'rehab_rates.csv');
if ~isfile(dataPath)
    % 命令行从任意目录运行时，用脚本自身位置定位项目；分节运行仍优先使用Current Folder。
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        projectDir = fileparts(fileparts(thisFile));
    end
    dataPath = fullfile(projectDir, '数据', 'rehab_rates.csv');
end
if ~isfile(dataPath), error('请把项目根目录设为 Current Folder，或从项目目录运行脚本。'); end
T = readtable(dataPath, 'TextType', 'string');
disp(T)
if height(T) ~= 2 || numel(unique(T.group)) ~= 2 || ~all(ismember(["new"; "control"], T.group)), error('group必须恰好有唯一的new和control两行。'); end
if ~isequal(T.Properties.VariableNames, {'group', 'success', 'failure'}), error('输入字段必须为group, success, failure。'); end
new = T(T.group == "new", :);
control = T(T.group == "control", :);
a = new.success;
b = new.failure;
c = control.success;
d = control.failure;
counts = [a b c d];
if any(~isfinite(counts)) || any(counts <= 0) || any(mod(counts,1) ~= 0), error('success/failure必须是有限非负整数；本练习同时计算RR和OR，请使用四格正整数；零计数需另选方法。'); end
n1 = a + b;
n0 = c + d;
p1 = a / n1;
p0 = c / n0;

%% 1 RD、RR、odds与OR
% RD是概率差；RR是概率比；OR是两组odds之比。
rd = p1 - p0;
rr = p1 / p0;
odds1 = p1 / (1 - p1);
odds0 = p0 / (1 - p0);
or_value = odds1 / odds0;
fprintf('RD=%.6f（%.1f个百分点），RR=%.6f，OR=%.6f\n', rd, rd*100, rr, or_value);
fprintf('odds(new)=%.6f，odds(control)=%.6f\n', odds1, odds0);

%% 2 从二项方差与变换斜率得到log尺度SE
var_p1 = p1 * (1 - p1) / n1;
var_p0 = p0 * (1 - p0) / n0;
slope_log_p1 = 1 / p1;
slope_log_p0 = 1 / p0;
slope_logit_p1 = 1 / (p1 * (1 - p1));
slope_logit_p0 = 1 / (p0 * (1 - p0));
var_log_p1 = slope_log_p1^2 * var_p1;
var_log_p0 = slope_log_p0^2 * var_p0;
var_logit_p1 = slope_logit_p1^2 * var_p1;
var_logit_p0 = slope_logit_p0^2 * var_p0;
se_log_rr = sqrt(var_log_p1 + var_log_p0);
se_log_or = sqrt(var_logit_p1 + var_logit_p0);
se_log_rr_cells = sqrt(1 / a - 1 / n1 + 1 / c - 1 / n0);
se_log_or_cells = sqrt(1 / a + 1 / b + 1 / c + 1 / d);
fprintf('var(p1)=%.8f, var(p0)=%.8f\n', var_p1, var_p0);
fprintf('log斜率=%.6f, %.6f；logit斜率=%.6f, %.6f\n', slope_log_p1, slope_log_p0, slope_logit_p1, slope_logit_p0);
fprintf('变换后方差：log p1=%.8f, log p0=%.8f；logit p1=%.8f, logit p0=%.8f\n', var_log_p1, var_log_p0, var_logit_p1, var_logit_p0);
fprintf('SE_log_RR=%.15f（四格公式=%.15f）；SE_log_OR=%.15f（四格公式=%.15f）\n', se_log_rr, se_log_rr_cells, se_log_or, se_log_or_cells);

%% 3 log尺度95%区间，再exp回原尺度
% erfcinv是基础MATLAB函数；这里等价于标准正态的0.975分位点。
z = -sqrt(2) * erfcinv(2 * 0.975); % 双侧95%，左右各留0.025
log_rr = log(rr);
log_or = log(or_value);
log_rr_lower = log_rr - z * se_log_rr;
log_rr_upper = log_rr + z * se_log_rr;
log_or_lower = log_or - z * se_log_or;
log_or_upper = log_or + z * se_log_or;
% exp把log尺度端点还原为原比值尺度。
rr_lower = exp(log_rr_lower);
rr_upper = exp(log_rr_upper);
or_lower = exp(log_or_lower);
or_upper = exp(log_or_upper);
fprintf('z=%.12f；RR原尺度=[%.8f, %.8f]；OR=[%.8f, %.8f]\n', z, rr_lower, rr_upper, or_lower, or_upper);

%% 4 交换两组：实际重算SE与log区间，再核对倒数端点
swapped_rr = p0 / p1;
swapped_or = odds0 / odds1;
swapped_se_log_rr = sqrt(1 / c - 1 / n0 + 1 / a - 1 / n1);
swapped_se_log_or = sqrt(1 / c + 1 / d + 1 / a + 1 / b);
swapped_log_rr = log(swapped_rr);
swapped_log_or = log(swapped_or);
swapped_log_rr_lower = swapped_log_rr - z * swapped_se_log_rr;
swapped_log_rr_upper = swapped_log_rr + z * swapped_se_log_rr;
swapped_log_or_lower = swapped_log_or - z * swapped_se_log_or;
swapped_log_or_upper = swapped_log_or + z * swapped_se_log_or;
swapped_rr_lower = exp(swapped_log_rr_lower);
swapped_rr_upper = exp(swapped_log_rr_upper);
swapped_or_lower = exp(swapped_log_or_lower);
swapped_or_upper = exp(swapped_log_or_upper);
rr_swap_error = max(abs(swapped_rr_lower - 1 / rr_upper), abs(swapped_rr_upper - 1 / rr_lower));
or_swap_error = max(abs(swapped_or_lower - 1 / or_upper), abs(swapped_or_upper - 1 / or_lower));
fprintf('交换后RR=%.8f，OR=%.8f；重算SE=%.15f, %.15f\n', swapped_rr, swapped_or, swapped_se_log_rr, swapped_se_log_or);
fprintf('RR区间倒数核对最大误差=%.3e；OR区间倒数核对最大误差=%.3e\n', rr_swap_error, or_swap_error);

%% 5 人数同时乘4：点估计不变，SE减半
scale = 4;
se_rr4 = sqrt(1/(scale*a) - 1/(scale*n1) + 1/(scale*c) - 1/(scale*n0));
se_or4 = sqrt(1/(scale*a) + 1/(scale*b) + 1/(scale*c) + 1/(scale*d));
fprintf('人数乘%g后SE_log_RR=%.15f，SE_log_OR=%.15f\n', scale, se_rr4, se_or4);

%% 可选保存：设置环境变量 RR_OR_SAVE_RESULTS=1 才写文件
metrics = ["p1";"p0";"RD";"RR";"OR";"z";"se_log_RR";"se_log_OR";"log_RR_lower";"log_RR_upper";"log_OR_lower";"log_OR_upper";"RR_lower";"RR_upper";"OR_lower";"OR_upper"];
values = [p1;p0;rd;rr;or_value;z;se_log_rr;se_log_or;log_rr_lower;log_rr_upper;log_or_lower;log_or_upper;rr_lower;rr_upper;or_lower;or_upper];
if SAVE_RESULTS
    outDir = fullfile(projectDir, '运行结果', 'MATLAB');
    if ~isfolder(outDir), mkdir(outDir); end
    writetable(table(metrics, values, 'VariableNames', {'metric','value'}), fullfile(outDir, 'results.csv'));
end
