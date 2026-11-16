%% 准备与参数：统计课堂29，独立样本率差A-B（UTF-8）
% 问题：假设A组9/20阳性，B组13/20阳性，估计总体率差双侧区间。
% Current Folder设为本项目根或代码目录，从本节起逐节运行。
% 无工具箱依赖。本文件逐步实现公式，属于教学代码，不是MATLAB内置区间函数。
% 默认命令窗口显示结果、Figure显示图；保存开关打开后写运行结果/MATLAB。
confidence = 0.95;
SAVE_RESULTS = false;
root = pwd;
if ~isfile(fullfile(root,'数据','independent_counts.csv'))
    root = fileparts(root);
end
data_path = fullfile(root,'数据','independent_counts.csv');

%% 问题1：两批独立的人，A率减B率是多少？
% 思路：先看共享CSV字段和汇总行，再按group取出A/B，不补造个体记录。
data = readtable(data_path,'TextType','string');
disp(data.Properties.VariableNames);
disp(data);
assert(all(ismember({'group','positive','total'},data.Properties.VariableNames)));
assert(height(data)==2 && all(sort(data.group)==["A";"B"]),'A/B须各一行');
counts = [data.positive,data.total];
assert(all(isfinite(counts),'all') && all(counts==floor(counts),'all'),'计数须为有限整数');
assert(all(data.total>0) && all(data.positive>=0) && all(data.positive<=data.total));
assert(isscalar(confidence) && isfinite(confidence) && confidence>0 && confidence<1);
% 表格括号中的逻辑条件按行选择，冒号表示保留全部列。
A = data(data.group=="A",:);
B = data(data.group=="B",:);
xA = A.positive;
nA = A.total;
xB = B.positive;
nB = B.total;
pA = xA/nA;
pB = xB/nB;
d = pA-pB;
fprintf('A阳性/阴性=%d/%d；B阳性/阴性=%d/%d。\n',xA,nA-xA,xB,nB-xB);
fprintf('A=%.2f%%，B=%.2f%%；A-B=%.2f个百分点。\n',100*pA,100*pB,100*d);

%% 问题2：Wald怎样从标准误得到区间？
% 独立两组的率差方差为各自方差相加；用样本率估计未知率。
alpha = 1-confidence;
% 标准正态CDF与erf的关系：双侧95%时累计概率.975，对应约1.96。
% erfinv是基础MATLAB函数，因此本项目不需要norminv的工具箱。
z = sqrt(2)*erfinv(2*(1-alpha/2)-1);
se = sqrt(pA*(1-pA)/nA+pB*(1-pB)/nB);
wald_lower = d-z*se;
wald_upper = d+z*se;
fprintf('Wald标准误=%.8f，正态临界值=%.6f。\n',se,z);
fprintf('Wald区间（百分点）：[%.6f, %.6f]。\n',100*wald_lower,100*wald_upper);

%% 问题3：由每组Wilson端点构造Newcombe区间
% 思路：先求每组率的Wilson中心和半宽；两侧距离再“平方相加开平方”。
% 这是Newcombe1998独立样本方法10，不作连续性校正。
centerA = (pA+z^2/(2*nA))/(1+z^2/nA);
halfA = z*sqrt(pA*(1-pA)/nA+z^2/(4*nA^2))/(1+z^2/nA);
LA = centerA-halfA;
UA = centerA+halfA;
centerB = (pB+z^2/(2*nB))/(1+z^2/nB);
halfB = z*sqrt(pB*(1-pB)/nB+z^2/(4*nB^2))/(1+z^2/nB);
LB = centerB-halfB;
UB = centerB+halfB;
% 左端组合A向下、B向上的距离；右端组合A向上、B向下的距离。
newcombe_lower = d-sqrt((pA-LA)^2+(UB-pB)^2);
newcombe_upper = d+sqrt((UA-pA)^2+(pB-LB)^2);
fprintf('A的Wilson区间：[%.6f, %.6f]；B：[%.6f, %.6f]。\n',LA,UA,LB,UB);

%% 问题4：Agresti-Caffo怎样调整，三种答案怎样解释？
% 每组加1阳性、1阴性：阳性数+1，总人数+2；CSV中的原始资料不变。
pA_adjusted = (xA+1)/(nA+2);
pB_adjusted = (xB+1)/(nB+2);
ac_center = pA_adjusted-pB_adjusted;
ac_se = sqrt(pA_adjusted*(1-pA_adjusted)/(nA+2)+pB_adjusted*(1-pB_adjusted)/(nB+2));
ac_lower = ac_center-z*ac_se;
ac_upper = ac_center+z*ac_se;
method = ["Wald";"Newcombe";"Agresti-Caffo"];
lower = [wald_lower;newcombe_lower;ac_lower];
upper = [wald_upper;newcombe_upper;ac_upper];
estimate = repmat(d,3,1); % 三行都报告同一个原始样本率差，AC中心另外说明。
lower_pp = 100*lower;
upper_pp = 100*upper;
width_pp = upper_pp-lower_pp;
results = table(method,estimate,lower,upper,lower_pp,upper_pp,width_pp);
disp(results);
fprintf('AC计算中心=%.4f个百分点；原始率差仍为%.2f。\n',100*ac_center,100*d);
fprintf('Newcombe的%.0f%%区间为%.2f至%.2f个百分点。\n',100*confidence,100*newcombe_lower,100*newcombe_upper);
if newcombe_lower<=0 && newcombe_upper>=0
    disp('该区间包含0，尚不能确定总体率差方向；不等于证明两总体率相等。');
else
    disp('该区间不含0，方向由区间正负判断；差异大小仍需结合研究背景。');
end
disp('这些是近似区间，不保证所有真实率下都恰好达到所标置信水平。');

%% 问题5：画出区间，点与线分别表示什么？
% 英文图签便于跨平台显示；横轴已乘100，单位为百分点。
figure('Position',[100 100 950 550]);
hold on;
% 每次循环画一种方法；点统一标原始估计，AC另加菱形标计算中心。
for k = 1:3
    plot([lower_pp(k),upper_pp(k)],[k,k],'Color',[0 0 .96],'LineWidth',3);
end
h_raw = scatter(repmat(100*d,3,1),(1:3)',40,'k','filled');
h_ac = scatter(100*ac_center,3,65,[.91 .59 .14],'d','filled');
xline(0,'--','Color',[.5 .5 .5]);
yticks(1:3);
yticklabels(method);
set(gca,'YDir','reverse');
ylim([.5 3.5]);
xlabel('Population difference A-B (percentage points)');
title(sprintf('%.0f%% two-sided CIs | A %d/%d, B %d/%d',100*confidence,xA,nA,xB,nB));
legend([h_raw,h_ac],{'Raw A-B','AC center'},'Location','best');
hold off;
if SAVE_RESULTS
    output = fullfile(root,'运行结果','MATLAB');
    if ~isfolder(output)
        mkdir(output);
    end
    writetable(results,fullfile(output,'intervals.csv'));
    exportgraphics(gcf,fullfile(output,'intervals.png'),'Resolution',160);
end
% 练习：①CSV改45/100、65/100；②恢复原数据，confidence改.90；③两组均0/20。
% 每次只改一项，从准备/读取节起顺序运行，观察点估计、区间与图标题。
