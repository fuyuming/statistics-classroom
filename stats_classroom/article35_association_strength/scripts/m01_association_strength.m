% 统计课堂35｜关联强度与区间。输入为第11课的假设工程试验汇总计数。
% 300台独立记录器随机分3组各100台：A硅胶、B环氧、C热缩封装。
% 弯折100次后接模拟器记录10min，无持续≥1s数据缺口判合格（课堂阈值）。
% 主表按A/B/C行、合格/不合格列排列为[[90,10],[80,20],[65,35]]。
% 在项目根目录按节顺序运行；默认不写文件，ARTICLE35_SAVE_RESULTS=1才保存。

%% 01.1 准备与输入
% 完整run可能暂时进入scripts；仅在找不到输入时，用本脚本位置回到项目根目录。
if ~isfile(fullfile('data','device_pass.csv'))
    task_script_path = mfilename('fullpath');
    task_project_root = fileparts(fileparts(task_script_path));
    cd(task_project_root);
end
addpath('helpers');                    % 从项目根目录加载本项目自编完整函数
raw = readtable('data/device_pass.csv');
disp(raw)
tab = table2array(raw(:, {'pass','fail'}));
assert(all(sum(tab,1)>0) && all(sum(tab,2)>0));
n = sum(tab, 'all');
r = size(tab,1);
s = size(tab,2);

%% 01.2 手算预期频数、Q、P、C与V
% 列向量乘行向量得到每个行列组合的合计乘积；点运算逐格平方/相除。
expected = sum(tab,2) * sum(tab,1) / n;
Q = sum((tab-expected).^2 ./ expected, 'all');
df = (r-1)*(s-1);
m = min(r-1,s-1);
p_value = chi2cdf(Q, df, 'upper');      % upper读取右尾概率
C = sqrt(Q/(Q+n));
Cmax = sqrt(m/(m+1));
V = sqrt(Q/(n*m));
point = table(n,r,s,df,m,Q,p_value,C,Cmax,V, ...
    'VariableNames', {'n','r','s','df','m','Q','p','C','Cmax','V'});
disp(expected)
disp(point)

%% 01.3 A/C子表：有符号phi与0/1相关
ac = tab([1,3], :);                    % MATLAB从1计数，第1/3行是A/C
n_sub = sum(ac,'all');
a = ac(1,1);
b = ac(1,2);
c = ac(2,1);
d = ac(2,2);
phi = (a*d-b*c)/sqrt(prod(sum(ac,1))*prod(sum(ac,2)));
% 按已有四格计数展开类别标签，核对0/1 Pearson相关。
x = repelem([1;1;0;0], [a;b;c;d]);
y = repelem([1;0;1;0], [a;b;c;d]);
phi_pearson = corr(x,y);
flipped = ac([2,1],:);
phi_flipped = det(flipped)/sqrt(prod(sum(flipped,1))*prod(sum(flipped,2)));
subfit = association_summary(ac);
Q_sub = subfit.Q;
V_sub = subfit.V;
phi_result = table(a,b,c,d,n_sub,Q_sub,phi,phi_pearson,phi_flipped,V_sub, ...
    'VariableNames', {'a','b','c','d','n_sub','Q_sub','phi','phi_pearson','phi_row_flipped','V_sub'});
assert(abs(phi-phi_pearson)<1e-12 && abs(phi_flipped+phi)<1e-12);
assert(abs(phi^2-Q_sub/n_sub)<1e-12);
disp(phi_result)

%% 01.4 相同比例与不同样本量的假设试验
scaling = table();
for k = [1,2,4]
    enlarged = k*tab;
    nn = sum(enlarged,'all');
    ee = sum(enlarged,2)*sum(enlarged,1)/nn;
    qq = sum((enlarged-ee).^2./ee,'all');
    pp = chi2cdf(qq,df,'upper');
    qqn = qq/nn;
    cc = sqrt(qq/(qq+nn));
    vv = sqrt(qq/(nn*m));
    row = table(k,nn,qq,pp,qqn,cc,vv, ...
        'VariableNames', {'multiplier','n','Q','p','Q_over_n','C','V'});
    scaling = [scaling; row];
end
assert(all(abs(scaling.V-V)<1e-12));
disp(scaling)

%% 01.5 完整自编函数核对与可选保存
% 函数从输入表重新计算，用sum(N.^2/E)-n这一等价Q式；不复用手算变量。
fit = association_summary(tab);
function_comparison = table([Q;p_value;C;V], [fit.Q;fit.p;fit.C;fit.V], ...
    'VariableNames', {'manual','full_function'});
assert(max(abs(function_comparison.manual-function_comparison.full_function))<1e-10);
disp(function_comparison)
if strcmp(getenv('ARTICLE35_SAVE_RESULTS'),'1')
    if ~isfolder('results')
        mkdir('results');
    end
    writetable(point,'results/MATLAB_point.csv');
    writetable(phi_result,'results/MATLAB_phi.csv');
    writetable(scaling,'results/MATLAB_scaling.csv');
end
