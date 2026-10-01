% ============================================================
% 精读 07｜过拟合与正则化（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 07 讲 Overfitting；教材第 7、8 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
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
function [b, e, sHat, sse, r2] = olsLS(X, y)
  b = (X' * X) \ (X' * y);
  e = y - X * b;
  sHat = sqrt(sum(e .^ 2) / (numel(y) - size(X, 2)));
  sse = sum(e .^ 2);
  r2 = 1 - sse / sum((y - mean(y)) .^ 2);
end
function out = looErr(X, y)
  XtXinv = inv(X' * X);
  h = sum((X * XtXinv) .* X, 2);
  [~, e] = olsLS(X, y);
  loo = e ./ (1 - h);
  out = sqrt(mean(loo .^ 2));
end
function v = orthResid(v, X)
  v = v - X * ((X' * X) \ (X' * v));
end
function X = polyXOf(base, k)
  X = [ones(numel(base), 1), base];
  for j = 2:k, X = [X, base .^ j]; end
end
function r = r4(x), r = round(x, 4); end

n = 12;
tX = (1:n)';
ts = (tX - mean(tX)) / std(tX);
tsOf = @(tv) (tv - mean(tX)) / std(tX);
muTrue = 20 + 4 * tX - 0.35 * tX .^ 2 + 0.010 * tX .^ 3;
residG = orthResid([0.8 -1.2 0.6 -0.9 1.1 -0.5 0.7 -1.0 0.4 0.9 -0.6 1.0]', [ones(n,1), tX, tX .^ 2, tX .^ 3]);
yObs = muTrue + residG;

fid = fopen(fullfile(RES, '01_数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","观测","生成均值"\n');
for i = 1:n
  fprintf(fid, '%d,%s,%s\n', i, numFmt(r4(yObs(i)),4), numFmt(r4(muTrue(i)),4));
end
fclose(fid);

%% 阶数与误差
fid = fopen(fullfile(RES, '02_阶数与误差_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"阶数","系数个数","样本内残差标准差","样本内平方误差和","留一交叉验证误差"\n');
bestK = 1; bestLoo = Inf;
Kv = zeros(1, 7); SSEv = zeros(1, 7); LOOk = zeros(1, 7);
for k = 1:7
  X = polyXOf(ts, k);
  [~, ~, sHat, sse] = olsLS(X, yObs);
  lv = looErr(X, yObs);
  fprintf(fid, '%d,%d,%s,%s,%s\n', k, k+1, numFmt(r4(sHat),4), numFmt(r4(sse),4), numFmt(r4(lv),4));
  Kv(k) = k; SSEv(k) = sse; LOOk(k) = lv;
  if lv < bestLoo, bestLoo = lv; bestK = k; end
end
fclose(fid);

tFine = linspace(1, 16, 200)';
truthFine = 20 + 4 * tFine - 0.35 * tFine .^ 2 + 0.010 * tFine .^ 3;
fitCurve = @(k, y) polyXOf(tsOf(tFine), k) * olsLS(polyXOf(ts, k), y);
f3 = fitCurve(3, yObs); f7 = fitCurve(7, yObs);
fid = fopen(fullfile(RES, '03_两条拟合曲线_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","真值形状","三阶","七阶"\n');
for i = 1:200
  fprintf(fid, '%s,%s,%s,%s\n', numFmt(r4(tFine(i)),4), numFmt(r4(truthFine(i)),4), numFmt(r4(f3(i)),4), numFmt(r4(f7(i)),4));
end
fclose(fid);

%% 正则化＝先验宽度
X1 = [ones(n,1), tX];
[~, ~, s1] = olsLS(X1, yObs);
aSeq = 15 + (0:199) * (30 / 199);
bSeq = 0 + (0:199) * (8 / 199);
priorList = {[3.0 0.2], [3.0 1.0], []};
priorNames = {'窄先验 N(3, 0.2)', '中等先验 N(3, 1.0)', '平坦先验'};
fid = fopen(fullfile(RES, '04_先验宽度_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"先验","后验均值斜率","后验标准差","预测第16期"\n');
for p = 1:3
  lp = zeros(200*200, 1); idx = 0;
  for ia = 1:200
    for ib = 1:200
      idx = idx + 1;
      mu = aSeq(ia) + bSeq(ib) * tX;
      v = -0.5 * sum(((yObs - mu) / s1) .^ 2) - n * (log(s1) + 0.5 * log(2 * pi));
      pb = priorList{p};
      if ~isempty(pb), v = v - 0.5 * ((bSeq(ib) - pb(1)) / pb(2)) ^ 2; end
      lp(idx) = v;
    end
  end
  w = exp(lp - max(lp)); w = w / sum(w);
  Bm = kron(ones(200,1), bSeq'); Bm = reshape(Bm, [], 1);   % 与 R 的 expand.grid 一致：a 变得最快
  Am = kron(aSeq(:), ones(200,1));
  bm = sum(Bm .* w); bs = sqrt(sum((Bm - bm) .^ 2 .* w));
  pred16 = sum((Am + Bm * 16) .* w);
  fprintf(fid, '%s,%s,%s,%s\n', priorNames{p}, numFmt(r4(bm),4), numFmt(r4(bs),4), numFmt(r4(pred16),4));
end
fclose(fid);

%% 厚尾 vs 正态
yOut2 = yObs; yOut2(9) = yOut2(9) + 22;
X3 = polyXOf(ts, 3);
[bn, en, sn] = olsLS(X3, yOut2);
z = en / sn;
wT = (4 + 1) ./ (4 + z .^ 2);
fid = fopen(fullfile(RES, '05_厚尾权重_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"观测","标准化残差","厚尾权重"\n');
for i = 1:n
  fprintf(fid, '%d,%s,%s\n', i, numFmt(r4(z(i)),4), numFmt(r4(wT(i)),4));
end
fclose(fid);
[bt, ~, st] = olsLS(X3 .* sqrt(wT), yOut2 .* sqrt(wT));
fid = fopen(fullfile(RES, '06_正态与厚尾_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"模型","三次项系数","残差标准差"\n');
fprintf(fid, '正态最小二乘,%s,%s\n', numFmt(r4(bn(4)),4), numFmt(r4(sn),4));
fprintf(fid, '厚尾（t, df=4，IRLS 一步）,%s,%s\n', numFmt(r4(bt(4)),4), numFmt(r4(st),4));
fclose(fid);

%% 图
f = figure('Position', [100 100 1140 660], 'Color', 'w'); hold on;
plot(Kv, SSEv, '-o', 'Color', [0.086 0.490 0.502], 'LineWidth', 1.4);
plot(Kv, LOOk .^ 2 * 12, '-o', 'Color', [0.757 0.275 0.173], 'LineWidth', 1.4);
xline(bestK, '--', 'Color', [0.604 0.659 0.659]);
title('样本内误差一路变小，样本外误差却先降后升'); xlabel('多项式阶数'); ylabel('误差');
legend({'样本内平方误差和','留一交叉验证误差（放大到同一尺度）'}, 'Location', 'northwest'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-阶数与误差.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1140 660], 'Color', 'w'); hold on;
scatter(tX, yObs, 40, [0.757 0.275 0.173], 'filled');
plot(tFine, f3, 'Color', [0.086 0.490 0.502], 'LineWidth', 1.4);
plot(tFine, f7, 'Color', [0.757 0.275 0.173], 'LineWidth', 1.4);
plot(tFine, truthFine, '--', 'Color', [0.604 0.659 0.659], 'LineWidth', 1.2);
xline(12, '--', 'Color', [0.604 0.659 0.659]);
title('3 阶与 7 阶：样本内都贴得住，出了数据范围就分道扬镳'); xlabel('时间'); ylabel('观测值');
legend({'观测','3 阶','7 阶','真值形状'}, 'Location', 'northwest'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-两阶对照.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 6 份 CSV 已写入\n');
