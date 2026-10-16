%% 问题：先算100尾的功效，再反求最小样本量。UTF-8。
% 教材数值另拟右侧研究；需要Statistics and Machine Learning Toolbox。
% 打开本文件，点击Run。由脚本路径定位项目，无须手工改变搜索路径。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir')
    mkdir(out);
end
parameters = readtable(fullfile(root,'数据','设计参数.csv'));
disp(parameters); % 教学参数，不是原始数据；体长单位cm
mu0 = parameters.mu0;
mu1 = parameters.mu1;
sigma = parameters.sigma;
alpha = parameters.alpha;
target = parameters.target_power;
delta = mu1-mu0;
d = delta/sigma;
%% 1. norminv输入左侧累计概率，返回标准正态分位数。
n = 100;
se = sigma/sqrt(n);
cutoff = mu0+norminv(1-alpha)*se;
% normcdf计算指定真实均数下未达到拒绝门槛的概率。
beta = normcdf((cutoff-mu1)/se);
power = 1-beta;
%% 2. ceil向上取整，再检查相邻两个整数。
n_continuous = ((norminv(1-alpha)+norminv(target))/d)^2;
n_min = ceil(n_continuous);
power_previous = normcdf(norminv(1-alpha)-d*sqrt(n_min-1),'upper');
power_minimum = normcdf(norminv(1-alpha)-d*sqrt(n_min),'upper');
quantity = ["se";"cutoff";"beta";"power";"d";"n_continuous";"n_min";"power_previous";"power_minimum"];
value = [se;cutoff;beta;power;d;n_continuous;n_min;power_previous;power_minimum];
result = table(quantity,value);
disp(result);
writetable(result,fullfile(out,'结果.csv'));
assert(power_previous < target && power_minimum >= target);
disp('功效是指定真实差异下的检出概率，不是H1为真的概率。');
%% 3. 改为双侧观察结果，三种情境的差值都是0.40 cm。
observed_n = [25;100;400];
observed_se = sigma./sqrt(observed_n); % ./表示各元素分别相除
p_two_sided = 2*normcdf(abs(0.40./observed_se),'upper');
ci_low = 0.40-norminv(.975)*observed_se;
ci_high = 0.40+norminv(.975)*observed_se;
observed = table(observed_n,repmat(0.40,3,1),p_two_sided,ci_low,ci_high,'VariableNames',{'n','difference','p_two_sided','ci_low','ci_high'});
disp(observed);
writetable(observed,fullfile(out,'观察结果情境.csv'));
