%% 统计课堂09：期望、方差、协方差
% 只需基础MATLAB，不需统计工具箱。保存为UTF-8，Run完整文件。
% 输入为课堂教学例子；模拟不是实验数据。输出在运行结果/MATLAB。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir')
    mkdir(out);
end
logfile = fullfile(out,'逐题结果.txt');
if exist(logfile,'file')
    delete(logfile); % 只替换本脚本的旧运行记录。
end
diary(logfile);
groups = readtable(fullfile(root,'数据','groups.csv'));
categories = readtable(fullfile(root,'数据','categories.csv'));
disp('共享输入：n为样本量，p为给定总体概率，observed为教学计数。');
disp(groups);
disp(categories);
assert(all(groups.n>0) && all(groups.p>0 & groups.p<1));
assert(height(categories)==4 && all(abs(categories.p-.25)<1e-12));

%% 1. 两独立样本率差：两份方差相加
n1 = groups.n(1);
n2 = groups.n(2);
p1 = groups.p(1);
p2 = groups.p(2);
mean_d = p1-p2;
var1 = p1*(1-p1)/n1;
var2 = p2*(1-p2)/n2;
var_d = var1+var2;
se_d = sqrt(var_d);
fprintf('\n问题1：期望=%.6f，方差=%.6f，标准误=%.6f。\n',mean_d,var_d,se_d);
disp('解释：两组独立，所以协方差为0；SE为9.19个百分点。');

%% 2. 分类计数的负协方差与Q
observed = categories.observed;
p = categories.p;
n = sum(observed);
expected = n*p;
% p*p转置给出各对概率乘积；对角线np(1-p)，非对角线-n*p_j*p_l。
cov_counts = n*(diag(p)-p*p');
disp('问题2：分类计数协方差矩阵，对角线为方差。');
disp(cov_counts);
fprintf('全部元素之和=计数总和的方差=%.6f。\n',sum(cov_counts,'all'));
% .^表示每个元素分别平方，./表示对应元素相除。
contribution = (observed-expected).^2 ./ expected;
Q = sum(contribution);
z1 = (observed(1)-observed(2))/sqrt(n/2);
z2 = (observed(3)-observed(4))/sqrt(n/2);
z3 = (observed(1)+observed(2)-observed(3)-observed(4))/sqrt(n);
z = [z1;z2;z3];
fprintf('问题3：Q=%.6f；三个方向平方和=%.6f。\n',Q,sum(z.^2));
disp(z.^2);
disp('解释：一次数据的Q=4；正确且预先给定的模型下E(Q)=3。');
% 可选矩阵核对三个比较的协方差；正文已经逐项展开一个例子。
a = [1,-1,0,0;0,0,1,-1;1,1,-1,-1];
a = a ./ sqrt([n/2;n/2;n]);
cov_z = a*cov_counts*a';
assert(max(abs(cov_z-eye(3)),[],'all')<1e-12);
assert(abs(Q-sum(z.^2))<1e-12);
assert(abs(sum(cov_counts,'all'))<1e-12);
name = {'mean_d';'var_d';'se_d';'var_count';'cov_count';'Q';'z1_sq';'z2_sq';'z3_sq';'EQ'};
value = [mean_d;var_d;se_d;cov_counts(1,1);cov_counts(1,2);Q;z.^2;length(p)-1];
writetable(table(name,value),fullfile(out,'theory.csv'));

%% 3. 用均匀随机数生成0、1与分类结果，无需统计工具箱
rng(20261011);
repeats = 20000;
% rand小于p时记为成功；每一行是一轮样本。
c1 = sum(rand(repeats,n1)<p1,2);
c2 = sum(rand(repeats,n2)<p2,2);
differences = c1/n1-c2/n2;
u = rand(repeats,n);
counts = zeros(repeats,4);
bounds = [0;cumsum(p)];
for j = 1:4
    counts(:,j) = sum(u>=bounds(j) & u<bounds(j+1),2);
end
q_sim = sum((counts-expected').^2 ./ expected',2);
writetable(table(differences,'VariableNames',{'difference'}),fullfile(out,'rate_simulation.csv'));
sim = array2table(counts,'VariableNames',{'A','B','C','D'});
sim.Q = q_sim;
writetable(sim,fullfile(out,'category_simulation.csv'));
fprintf('模拟：率差均值=%.6f，方差=%.6f，Q均值=%.6f。\n',mean(differences),var(differences),mean(q_sim));
disp('var默认n-1分母；有限模拟有误差，同样种子也不保证跨语言序列相同。');
disp('练习：计数改为10、20、25、25，应得到Q=7.5；概率仍为各0.25。');
% 英文标签用于跨系统显示；中文公众号图由配套Python脚本生成。
figure('Color','w');
histogram(q_sim,0:.5:ceil(max(q_sim)+1),'Normalization','pdf','FaceColor',[.65 .8 .78]);
hold on;
x = linspace(.001,20,501);
% 卡方(3)密度的定义，避免依赖chi2pdf工具箱函数。
density = x.^(.5).*exp(-x/2)/(2^1.5*gamma(1.5));
plot(x,density,'LineWidth',2);
xlim([0,20]);
xlabel('Pearson Q');
ylabel('Density');
title('Multinomial simulation and chi-square(3)');
exportgraphics(gcf,fullfile(out,'Q_check.png'),'Resolution',160);
diary off;
