% ============================================================
% 精读 08｜MCMC（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 08 讲 Markov Chain Monte Carlo；教材第 9 章
%
% 设计：MCMC 天生要用随机数，而本系列要求三个语言给出同一个答案。
% 这里用确定性低差异序列（Halton）代替伪随机数：提议步长由 norminv(Halton) 给出，
% 接受判断用 Halton 的均匀序列。算法结构与真实 Metropolis 一致，但结果可复现。
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
function h = halton(n, base)
  h = zeros(n, 1);
  for i = 1:n
    f = 1; r = 0; k = i;
    while k > 0
      f = f / base; r = r + f * mod(k, base); k = floor(k / base);
    end
    h(i) = r;
  end
end
function q = wQuant(v, w, p)
  [v2, o] = sort(v); w2 = w(o);
  cw = cumsum(w2) / sum(w2);
  k = find(cw >= p, 1);
  if k == 1, q = v2(1); return; end
  q = v2(k-1) + (v2(k) - v2(k-1)) * (p - cw(k-1)) / (cw(k) - cw(k-1));
end
function v = orthResid(v, X)
  v = v - X * ((X' * X) \ (X' * v));
end
function r = r4(x), r = round(x, 4); end

n = 16;
x = (1:n)';
residY = orthResid([0.9 -1.1 0.7 -0.8 1.2 -0.4 0.6 -1.0 0.5 1.1 -0.7 0.8 -0.9 0.4 -0.6 1.0]', [ones(n,1), x]);
y = 12 + 2.5 * x + residY;

fid = fopen(fullfile(RES, '01_数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","时间","观测"\n');
for i = 1:n
  fprintf(fid, '%d,%d,%s\n', i, x(i), numFmt(r4(y(i)),4));
end
fclose(fid);

X0 = [ones(n,1), x];
bOls = (X0' * X0) \ (X0' * y);
sY = sqrt(sum((y - X0 * bOls) .^ 2) / (n - 2));
logLik = @(a, b) -0.5 * sum(((y - (a + b * x)) / sY) .^ 2);

%% 网格后验（对照）
aSeq = 0 + (0:299) * (24 / 299);
bSeq = 0 + (0:299) * (5 / 299);
[Ag, Bg] = meshgrid(aSeq, bSeq);
Ag = Ag'; Bg = Bg';            % 与 R 的 expand.grid 同序：a 变得最快
Ag = Ag(:); Bg = Bg(:);
lp = zeros(numel(Ag), 1);
for i = 1:numel(Ag), lp(i) = logLik(Ag(i), Bg(i)); end
gw = exp(lp - max(lp)); gw = gw / sum(gw);
gridRow = {'网格（300×300 = 9 万格点）', r4(sum(Ag .* gw)), r4(wQuant(Ag, gw, 0.055)), r4(wQuant(Ag, gw, 0.945)), ...
           r4(sum(Bg .* gw)), r4(wQuant(Bg, gw, 0.055)), r4(wQuant(Bg, gw, 0.945))};

%% Metropolis
nStep = 20000; burn = 2000;
function [chain, accept] = metropolisF(start, nStep, stepA, stepB, offset, logLikF)
  chain = zeros(nStep, 2); chain(1, :) = start;
  z1 = norminv(halton(nStep, 2)); z2 = norminv(halton(nStep, 3)); u = halton(nStep, 5);
  lpCur = logLikF(chain(1,1), chain(1,2)); nAcc = 0;
  for i = 2:nStep
    j = mod(i + offset - 2, nStep) + 1;
    prop = chain(i-1, :) + [stepA * z1(j), stepB * z2(j)];
    lpProp = logLikF(prop(1), prop(2));
    if log(u(j)) < (lpProp - lpCur)
      chain(i, :) = prop; lpCur = lpProp; nAcc = nAcc + 1;
    else
      chain(i, :) = chain(i-1, :);
    end
  end
  accept = nAcc / (nStep - 1);
end
[ch1, acc1] = metropolisF([3 0.5], nStep, 0.9, 0.9, 0, logLik);
[ch2, acc2] = metropolisF([3 0.5], nStep, 0.25, 0.025, 0, logLik);
[ch3, acc3] = metropolisF([3 0.5], nStep, 0.05, 0.005, 0, logLik);
[ch2b, ~]  = metropolisF([20 4.0], nStep, 0.25, 0.025, 777, logLik);
keep1 = ch1(burn+1:end, :); keep2 = ch2(burn+1:end, :); keep3 = ch3(burn+1:end, :); keep2b = ch2b(burn+1:end, :);
allKeep = [keep2; keep2b];

acfOf = @(v, lag) sum((v(lag+1:end) - mean(v)) .* (v(1:end-lag) - mean(v))) / sum((v - mean(v)) .^ 2);
acfBig = zeros(20,1); acfGood = zeros(20,1); acfSmall = zeros(20,1);
for l = 1:20
  acfBig(l) = r4(acfOf(keep1(:,2), l)); acfGood(l) = r4(acfOf(keep2(:,2), l)); acfSmall(l) = r4(acfOf(keep3(:,2), l));
end
wVar = mean([var(keep2(:,2)), var(keep2b(:,2))]);
bVar = var([mean(keep2(:,2)), mean(keep2b(:,2))]);
rHat = sqrt((((size(keep2,1)-1)/size(keep2,1)) * wVar + bVar/size(keep2,1)) / wVar);
nEff = size(keep2,1) * (1 - acfOf(keep2(:,2),1)) / (1 + acfOf(keep2(:,2),1));

wEq = ones(numel(allKeep(:,2)),1) / numel(allKeep(:,2));
mcmcRow = {'Metropolis（两条链，预热 2000 步，各留 18000 步）', r4(mean(allKeep(:,1))), r4(wQuant(allKeep(:,1), wEq, 0.055)), ...
           r4(wQuant(allKeep(:,1), wEq, 0.945)), r4(mean(allKeep(:,2))), r4(wQuant(allKeep(:,2), wEq, 0.055)), r4(wQuant(allKeep(:,2), wEq, 0.945))};
fid = fopen(fullfile(RES, '02_网格与MCMC对照_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"方法","截距均值","截距89下","截距89上","斜率均值","斜率89下","斜率89上"\n');
for rowIdx = 1:2
  row = gridRow; if rowIdx == 2, row = mcmcRow; end
  fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n', row{1}, numFmt(row{2},4), numFmt(row{3},4), numFmt(row{4},4), numFmt(row{5},4), numFmt(row{6},4), numFmt(row{7},4));
end
fclose(fid);

fid = fopen(fullfile(RES, '03_诊断_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值"\n');
fprintf(fid, '步长太大：接受率,%s\n', numFmt(r4(acc1),4));
fprintf(fid, '合适步长：接受率,%s\n', numFmt(r4(acc2),4));
fprintf(fid, '步长太小：接受率,%s\n', numFmt(r4(acc3),4));
fprintf(fid, '步长太大：滞后 1 自相关,%s\n', numFmt(acfBig(1),4));
fprintf(fid, '合适步长：滞后 1 自相关,%s\n', numFmt(acfGood(1),4));
fprintf(fid, '步长太小：滞后 1 自相关,%s\n', numFmt(acfSmall(1),4));
fprintf(fid, '合适步长：近似有效样本量（3.6 万步）,%s\n', numFmt(r4(nEff),4));
fprintf(fid, '简化 R-hat（斜率，两条合适步长的链）,%s\n', numFmt(r4(rHat),4));
fprintf(fid, '斜率后验均值（链 1）,%s\n', numFmt(r4(mean(keep2(:,2))),4));
fprintf(fid, '斜率后验均值（链 2）,%s\n', numFmt(r4(mean(keep2b(:,2))),4));
fclose(fid);

fid = fopen(fullfile(RES, '04_维度诅咒_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"参数个数","格点数"\n');
for k = 1:8
  fprintf(fid, '%d,%s\n', k, numFmt(r4(100.0 ^ k),4));
end
fclose(fid);

run1 = cumsum(keep2(:,2)) ./ (1:size(keep2,1))';
run2 = cumsum(keep2b(:,2)) ./ (1:size(keep2b,1))';
fid = fopen(fullfile(RES, '05_运行均值_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"步数","链1运行均值","链2运行均值"\n');
for st = 100:500:size(keep2,1)
  fprintf(fid, '%d,%s,%s\n', st, numFmt(r4(run1(st)),4), numFmt(r4(run2(st)),4));
end
fclose(fid);

%% 图
f = figure('Position', [100 100 1110 690], 'Color', 'w'); hold on;
scatter(Ag(1:37:end), Bg(1:37:end), 2, [0.933 0.961 0.961], 'filled');
scatter(allKeep(1:20:end,1), allKeep(1:20:end,2), 3, [0.757 0.275 0.173], 'filled');
title('MCMC 做的事：在后验上走一圈'); xlabel('截距 a'); ylabel('斜率 b'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-MCMC与网格.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
plot(1:20, acfBig, '-o', 'Color', [0.604 0.659 0.659], 'LineWidth', 1.2, 'MarkerSize', 3);
plot(1:20, acfGood, '-o', 'Color', [0.086 0.490 0.502], 'LineWidth', 1.2, 'MarkerSize', 3);
plot(1:20, acfSmall, '-o', 'Color', [0.757 0.275 0.173], 'LineWidth', 1.2, 'MarkerSize', 3);
title('自相关：步长太大或太小都不好'); xlabel('滞后'); ylabel('自相关');
legend({'步长太大','合适步长','步长太小'}, 'Location', 'northeast'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-自相关.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 5 份 CSV 已写入\n');
