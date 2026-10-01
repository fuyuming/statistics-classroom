% ============================================================
% 精读 05｜混杂的四种基本形状（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 05 讲 Elemental Confounds；教材第 5、6 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
% 前置：先跑 R 版生成 运行结果/00_教材数据_华夫饼与离婚率.csv
% ============================================================
clear; clc;
root = getenv('RET_ROOT');
if isempty(root) || ~isfolder(fullfile(root, '代码'))
  root = fileparts(fileparts(mfilename('fullpath')));
  if isempty(root) || ~isfolder(fullfile(root, '代码')), root = pwd; end
end
cd(root);
OUT = fullfile(root, '文章配图'); RES = fullfile(root, '运行结果');
if ~exist(OUT, 'dir'), mkdir(OUT); end
if ~exist(RES, 'dir'), mkdir(RES); end

function s = numFmt(x, d)
  if nargin < 2, d = 6; end
  s = sprintf('%.*f', d, round(x, d));
  s = regexprep(s, '0+$', ''); s = regexprep(s, '\.$', '');
  if isempty(s) || strcmp(s, '-0'), s = '0'; end
end
function [b, sHat] = olsLS(X, y)
  b = (X' * X) \ (X' * y);
  res = y - X * b;
  sHat = sqrt(sum(res .^ 2) / (numel(y) - size(X, 2)));
end
function r = pearsonR(x, y)
  x = x - mean(x); y = y - mean(y);
  r = sum(x .* y) / sqrt(sum(x .^ 2) * sum(y .^ 2));
end
function r = r4(x), r = round(x, 4); end

%% 0) 教材对照
T = readtable(fullfile(RES, '00_教材数据_华夫饼与离婚率.csv'), 'Encoding', 'UTF-8', 'VariableNamingRule', 'preserve');
div = T{:, 2}; waff = T{:, 3}; south = T{:, 5};
rAll = pearsonR(waff, div); rS = pearsonR(waff(south == 1), div(south == 1)); rN = pearsonR(waff(south == 0), div(south == 0));
fid = fopen(fullfile(RES, '01_教材对照_分层相关_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"范围","华夫饼店与离婚率的相关","州数"\n');
fprintf(fid, '全体 50 个州,%s,%d\n', numFmt(r4(rAll), 4), height(T));
fprintf(fid, '南方各州,%s,%d\n', numFmt(r4(rS), 4), sum(south == 1));
fprintf(fid, '非南方各州,%s,%d\n', numFmt(r4(rN), 4), sum(south == 0));
fclose(fid);

%% 2) 叉（混杂）
treatC = [1 1 1 1 0 1 0 0 0 0]'; batchC = [0 0 0 0 0 1 1 1 1 1]';
residC = [0.6 -0.8 1.1 -0.5 0.3 0.9 -1.2 0.4 -0.6 0.8]';
growthC = 10 + 3.0 * treatC + 2.5 * batchC + residC;
fid = fopen(fullfile(RES, '02_我们的案例_叉_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","批次","长势"\n');
for i = 1:10
  bname = '批次1'; if batchC(i) == 1, bname = '批次2'; end
  fprintf(fid, '%d,%d,%s,%s\n', i, treatC(i), bname, numFmt(growthC(i)));
end
fclose(fid);
fNaive = olsLS([ones(10,1), treatC], growthC);
fAdj = olsLS([ones(10,1), treatC, batchC], growthC);
fid = fopen(fullfile(RES, '03_叉_调整前后_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","处理效应估计","真值"\n');
fprintf(fid, '不调整（长势 ~ 处理）,%s,3\n', numFmt(r4(fNaive(2)), 4));
fprintf(fid, '把批次放进模型（长势 ~ 处理 + 批次）,%s,3\n', numFmt(r4(fAdj(2)), 4));
fclose(fid);
fprintf('  叉：不调整 %s，调整 %s，处理与批次相关 %s\n', numFmt(r4(fNaive(2)),4), numFmt(r4(fAdj(2)),4), numFmt(r4(pearsonR(treatC, batchC)),4));

%% 3) 管（中介）
residD = [0.4 -0.6 0.9 -0.3 0.5 -0.8 0.7 -0.2 0.6 -0.9]';
residG = [0.7 0.5 -0.8 0.9 -0.4 0.3 -0.6 0.8 -0.5 0.2]';
densityD = 5 + 2.0 * treatC + residD;
growthD = 8 + 1.5 * densityD + residG;
fid = fopen(fullfile(RES, '04_我们的案例_管_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","菌群密度","长势"\n');
for i = 1:10
  fprintf(fid, '%d,%d,%s,%s\n', i, treatC(i), numFmt(round(densityD(i),6)), numFmt(round(growthD(i),6)));
end
fclose(fid);
pTot = olsLS([ones(10,1), treatC], growthD);
pDir = olsLS([ones(10,1), treatC, densityD], growthD);
fid = fopen(fullfile(RES, '05_管_总效应与直接效应_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","处理效应估计","说明"\n');
fprintf(fid, '不控制中介（长势 ~ 处理）,%s,总效应\n', numFmt(r4(pTot(2)), 4));
fprintf(fid, '控制中介（长势 ~ 处理 + 菌群密度）,%s,直接效应（总效应的一部分被中介带走）\n', numFmt(r4(pDir(2)), 4));
fclose(fid);

%% 4) 对撞（选择）
baseB = [-1.8 -1.2 -0.7 -0.3 0.1 0.4 0.9 1.3 1.7 2.0]';
treatB = [0 1 0 1 0 1 0 1 0 1]';
residS = [0.5 -0.4 0.2 0.8 -0.9 0.3 -0.2 0.6 -0.5 0.1]';
scoreS = 0.6 * treatB + 1.2 * baseB + residS;
selected = double(scoreS > median(scoreS));
fid = fopen(fullfile(RES, '06_我们的案例_对撞_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","基线值","选中","倾向分"\n');
for i = 1:10
  fprintf(fid, '%d,%d,%s,%d,%s\n', i, treatB(i), numFmt(baseB(i)), selected(i), numFmt(round(scoreS(i),6)));
end
fclose(fid);
rcAll = r4(pearsonR(treatB, baseB));
selIdx = find(selected == 1);
rcSel = r4(pearsonR(treatB(selIdx), baseB(selIdx)));
fid = fopen(fullfile(RES, '07_对撞_条件化前后_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"范围","处理与基线值的相关","单元数"\n');
fprintf(fid, '全部 10 个单元,%s,10\n', numFmt(rcAll, 4));
fprintf(fid, '只看被选中的 5 个单元,%s,%d\n', numFmt(rcSel, 4), numel(selIdx));
fclose(fid);

%% 5) 后代（治疗后变量）
residP = [0.3 -0.5 0.8 -0.2 0.6 -0.7 0.4 -0.1 0.5 -0.8]';
residY = [0.6 0.4 -0.7 0.8 -0.3 0.2 -0.5 0.7 -0.4 0.1]';
postD = 3 + 1.6 * treatC + residP;
growthP = 9 + 2.8 * treatC + 0.9 * postD + residY;
fid = fopen(fullfile(RES, '08_我们的案例_后代_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","处理后指标","长势"\n');
for i = 1:10
  fprintf(fid, '%d,%d,%s,%s\n', i, treatC(i), numFmt(round(postD(i),6)), numFmt(round(growthP(i),6)));
end
fclose(fid);
dTot = olsLS([ones(10,1), treatC], growthP);
dCtl = olsLS([ones(10,1), treatC, postD], growthP);
truth = 2.8 + 0.9 * 1.6;
fid = fopen(fullfile(RES, '09_后代_控制前后_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","处理效应估计","真值"\n');
fprintf(fid, '不控制后代（长势 ~ 处理）,%s,%s\n', numFmt(r4(dTot(2)),4), numFmt(r4(truth),4));
fprintf(fid, '把后代当协变量（长势 ~ 处理 + 处理后指标）,%s,%s\n', numFmt(r4(dCtl(2)),4), numFmt(r4(truth),4));
fclose(fid);

%% 图（两张对照图）
f = figure('Position', [100 100 1050 630], 'Color', 'w');
b = bar([r4(fAdj(2)), r4(fNaive(2))], 0.55);
b.FaceColor = 'flat'; b.CData(1,:) = [0.086 0.490 0.502]; b.CData(2,:) = [0.910 0.706 0.659];
hold on; yline(3.0, '--', 'Color', [0.757 0.275 0.173]);
text(1, r4(fAdj(2)) + 0.1, numFmt(r4(fAdj(2)),2), 'HorizontalAlignment','center');
text(2, r4(fNaive(2)) + 0.1, numFmt(r4(fNaive(2)),2), 'HorizontalAlignment','center');
xticklabels({'把批次放进模型（长势 ~ 处理 + 批次）', '不调整（长势 ~ 处理）'});
title('叉（混杂）：不调整批次会把效应估偏'); ylabel('处理效应估计'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-叉.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
scatter(baseB(selected==0), treatB(selected==0), 60, [0.086 0.490 0.502], 'filled');
scatter(baseB(selIdx), treatB(selIdx), 60, [0.757 0.275 0.173], 'filled');
xs = linspace(-1.9, 2.1, 50)';
k1 = polyfit(baseB(selIdx), treatB(selIdx), 1); plot(xs, polyval(k1, xs), 'Color', [0.757 0.275 0.173], 'LineWidth', 1.4);
k2 = polyfit(baseB, treatB, 1); plot(xs, polyval(k2, xs), '--', 'Color', [0.604 0.659 0.659], 'LineWidth', 1.4);
title('对撞（选择）：只看被选中的样本，就出现了假的负相关'); xlabel('基线值'); ylabel('处理（1=接种，0=对照）');
legend({'未被选中','被选中','被选中子集趋势','全体趋势'}, 'Location','southeast'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-对撞.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 10 份 CSV 已写入\n');
