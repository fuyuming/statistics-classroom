% ============================================================
% 精读 04｜分类变量、中心化与曲线（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 04 讲 Categories and Curves；教材第 4、5 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
% 前置：先跑 R 版生成 运行结果/00_教材数据_Howell1*.csv
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
set(0, 'defaultAxesFontName', 'PingFang SC'); set(0, 'defaultTextFontName', 'PingFang SC');
TEAL = [0.086 0.490 0.502]; RED = [0.757 0.275 0.173]; DARK = [0.086 0.239 0.282]; GREY = [0.604 0.659 0.659];

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

%% 0) 教材数据
hw = readmatrix(fullfile(RES, '00_教材数据_Howell1成年人.csv'));
allHw = readmatrix(fullfile(RES, '00_教材数据_Howell1全部.csv'));
ageAll = allHw(:, 1); wtAll = allHw(:, 2);

%% 1) 三组终点厚度
labels = {'对照', '低剂量', '高剂量'};
vals = {[41.2 43.0 39.8 42.5 40.9], [46.5 48.1 45.2 47.4 46.0], [52.8 54.6 51.9 53.7 52.2]};
fid = fopen(fullfile(RES, '01_我们的案例_三组终点厚度_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"处理","厚度"\n');
for g = 1:3
  for v = vals{g}, fprintf(fid, '%s,%s\n', labels{g}, numFmt(v)); end
end
fclose(fid);

gMean = zeros(3,1); gN = zeros(3,1); ss = 0;
for g = 1:3
  gMean(g) = mean(vals{g}); gN(g) = numel(vals{g});
  ss = ss + sum((vals{g} - gMean(g)) .^ 2);
end
dfWithin = sum(gN) - 3;
sigmaPool = sqrt(ss / dfWithin);
postSd = sigmaPool ./ sqrt(gN);
z89 = norminv(0.945);

fid = fopen(fullfile(RES, '02_组均值与后验区间_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"组","样本均值","后验标准差","区间89下","区间89上"\n');
for g = 1:3
  fprintf(fid, '%s,%s,%s,%s,%s\n', labels{g}, numFmt(round(gMean(g),4),4), numFmt(round(postSd(g),4),4), ...
          numFmt(round(gMean(g) - z89*postSd(g),4),4), numFmt(round(gMean(g) + z89*postSd(g),4),4));
end
fclose(fid);

pairNames = {'低剂量 - 对照', '高剂量 - 对照', '高剂量 - 低剂量'};
pairA = [2 3 3]; pairB = [1 1 2];
fid = fopen(fullfile(RES, '03_组间对比_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"对比","均值差","差的标准误","区间89下","区间89上","P大于0"\n');
for i = 1:3
  d = gMean(pairA(i)) - gMean(pairB(i));
  se = sqrt(postSd(pairA(i))^2 + postSd(pairB(i))^2);
  dR = round(d, 4); seR = round(se, 4);     % 先舍入再算区间，与 R／Python 的算法一致
  fprintf(fid, '%s,%s,%s,%s,%s,%s\n', pairNames{i}, numFmt(round(d,4),4), numFmt(round(se,4),4), ...
          numFmt(round(dR - z89*seR,4),4), numFmt(round(dR + z89*seR,4),4), numFmt(round(normcdf(d/se),4),4));
end
fclose(fid);

%% 2) 饱和曲线
tHours = (1:12)';
residFix = [0.5; -1.1; 0.9; -0.6; 1.2; -0.4; 0.7; -1.0; 0.3; 1.1; -0.8; 0.6];
satu = 58 * (1 - exp(-0.32 * tHours)) + residFix;
fid = fopen(fullfile(RES, '04_饱和曲线数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","厚度"\n');
for i = 1:12, fprintf(fid, '%d,%s\n', i, numFmt(round(satu(i),6))); end
fclose(fid);

dLin  = @(t) [ones(numel(t),1), t];
dPoly = @(t) [ones(numel(t),1), t, t.^2, t.^3];
dSpl  = @(t) [ones(numel(t),1), t, max(t - 4, 0), max(t - 8, 0)];
[bLin, sLin]   = olsLS(dLin(tHours), satu);
[bPoly, sPoly] = olsLS(dPoly(tHours), satu);
[bSpl, sSpl]   = olsLS(dSpl(tHours), satu);

fid = fopen(fullfile(RES, '05_三种拟合_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","观测","直线","三次多项式","线性样条"\n');
for i = 1:12
  fprintf(fid, '%d,%s,%s,%s,%s\n', i, numFmt(round(satu(i),6)), ...
          numFmt(round(dLin(tHours(i))*bLin, 4),4), numFmt(round(dPoly(tHours(i))*bPoly, 4),4), ...
          numFmt(round(dSpl(tHours(i))*bSpl, 4),4));
end
fclose(fid);

fid = fopen(fullfile(RES, '06_拟合优度_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"模型","参数个数","残差标准差"\n');
fprintf(fid, '直线,2,%s\n', numFmt(round(sLin,4),4));
fprintf(fid, '三次多项式,4,%s\n', numFmt(round(sPoly,4),4));
fprintf(fid, '线性样条（结点 4、8）,4,%s\n', numFmt(round(sSpl,4),4));
fclose(fid);

fid = fopen(fullfile(RES, '07_外推四小时_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","直线","三次多项式","线性样条"\n');
for t = 13:16
  fprintf(fid, '%d,%s,%s,%s\n', t, numFmt(round(dLin(t)*bLin,4),4), numFmt(round(dPoly(t)*bPoly,4),4), ...
          numFmt(round(dSpl(t)*bSpl,4),4));
end
fclose(fid);

Xl = dLin(ageAll); Xp = dPoly(ageAll);
Xs = [ones(numel(ageAll),1), ageAll, max(ageAll - 10, 0), max(ageAll - 15, 0), max(ageAll - 25, 0), max(ageAll - 40, 0)];
fid = fopen(fullfile(RES, '08_教材对照_拟合优度_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"模型","参数个数","残差标准差"\n');
[~, sL] = olsLS(Xl, wtAll);
fprintf(fid, '直线,2,%s\n', numFmt(round(sL, 4), 4));
[~, sP] = olsLS(Xp, wtAll); [~, sS] = olsLS(Xs, wtAll);
fprintf(fid, '三次多项式,4,%s\n', numFmt(round(sP,4),4));
fprintf(fid, '线性样条（结点 10/15/25/40）,6,%s\n', numFmt(round(sS,4),4));
fclose(fid);

%% 中心化
tbar = mean(tHours); nObs = numel(tHours); Sxx = sum((tHours - tbar) .^ 2);
sdRaw = sLin * sqrt(1/nObs + tbar^2/Sxx); sdCtr = sLin / sqrt(nObs);
fid = fopen(fullfile(RES, '09_中心化对照_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","截距后验均值","截距后验标准差"\n');
fprintf(fid, '未中心化：截距＝时间 0 处的厚度,%s,%s\n', numFmt(round(bLin(1),4),4), numFmt(round(sdRaw,4),4));
fprintf(fid, '中心化：截距＝平均时间处的厚度,%s,%s\n', numFmt(round(bLin(1) + bLin(2)*tbar,4),4), numFmt(round(sdCtr,4),4));
fclose(fid);

%% 图（与 R 版同内容）
f = figure('Position', [100 100 990 600], 'Color', 'w'); hold on;
rng(1);
for g = 1:3
  xj = g + (rand(gN(g),1) - 0.5) * 0.16;
  plot(xj, vals{g}, 'o', 'MarkerSize', 6, 'MarkerFaceColor', RED, 'MarkerEdgeColor', RED);
  plot([g g], [gMean(g) - z89*postSd(g), gMean(g) + z89*postSd(g)], 'Color', DARK, 'LineWidth', 1.4);
  plot(g, gMean(g), 'o', 'MarkerSize', 9, 'MarkerFaceColor', TEAL, 'MarkerEdgeColor', TEAL);
end
xticks(1:3); xticklabels(labels);
title('三个处理组的终点厚度：先把「组均值」的后验画出来'); ylabel('终点厚度（微米）'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-三组终点厚度.png'), 'Resolution', 150); close(f);

xs = linspace(-12, 14, 400);
f = figure('Position', [100 100 1110 630], 'Color', 'w'); hold on;
for i = 1:3
  d = gMean(pairA(i)) - gMean(pairB(i)); se = sqrt(postSd(pairA(i))^2 + postSd(pairB(i))^2);
  plot(xs, normpdf(xs, d, se), 'LineWidth', 1.4);
end
xline(0, '--', 'Color', GREY);
title('报系数不如报对比：三组两两差多少'); xlabel('厚度差（微米）'); ylabel('后验密度');
legend(pairNames, 'Location', 'northeast'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-组间对比后验.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 9 份 CSV 已写入\n');
