% ============================================================
% 精读 03｜12 个时间点，就敢画一条生长曲线吗（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 03 讲 Geocentric Models；教材第 4 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己 cd 进项目目录）
% 前置：先跑 R 版生成 运行结果/00_教材数据_Howell1成年人.csv
% ============================================================
clear; clc;
% 目录自动定位：优先读环境变量 RET_ROOT（跑批用），否则按脚本位置推（读者直接用这一条）
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
A_LIM = [10 26]; B_LIM = [2.0 6.0]; N_GRID = 200;

function s = numFmt(x, d)
  if nargin < 2, d = 6; end
  s = sprintf('%.*f', d, round(x, d));
  s = regexprep(s, '0+$', ''); s = regexprep(s, '\.$', '');
  if isempty(s) || strcmp(s, '-0'), s = '0'; end
end
function u = halton(n, base)
  u = zeros(n, 1);
  for i = 1:n
    f = 1; r = 0; k = i;
    while k > 0
      f = f / base; r = r + f * mod(k, base); k = floor(k / base);
    end
    u(i) = r;
  end
end
function q = wQuant(v, w, p)              % 加权分位数（线性插值），三语言同法
  v = v(:); w = w(:);
  [vs, o] = sort(v); ws = w(o);
  cw = cumsum(ws) / sum(ws);
  k = find(cw >= p, 1);
  if k == 1, q = vs(1); return; end
  q = vs(k-1) + (vs(k) - vs(k-1)) * (p - cw(k-1)) / (cw(k) - cw(k-1));
end
function [aSeq, bSeq, post] = gridPost(x, y, nGrid, sigma, aLim, bLim, priorA, priorB)
  if nargin < 5 || isempty(aLim), aLim = A_LIM; end
  if nargin < 6 || isempty(bLim), bLim = B_LIM; end
  if nargin < 7, priorA = []; end
  if nargin < 8, priorB = []; end
  stepA = (aLim(2) - aLim(1)) / (nGrid - 1);
  stepB = (bLim(2) - bLim(1)) / (nGrid - 1);
  aSeq = aLim(1) + (0:nGrid-1)' * stepA;
  bSeq = bLim(1) + (0:nGrid-1)' * stepB;
  x = x(:); y = y(:); n = numel(y);
  M = x * bSeq';
  lp = zeros(nGrid * nGrid, 1);
  for i = 1:(nGrid * nGrid)
    ai = mod(i - 1, nGrid) + 1; bi = floor((i - 1) / nGrid) + 1;
    e = (y - (aSeq(ai) + M(:, bi))) / sigma;
    lp(i) = -0.5 * sum(e .^ 2) - n * (log(sigma) + 0.5 * log(2 * pi));
    if ~isempty(priorA)
      lp(i) = lp(i) - 0.5 * ((aSeq(ai) - priorA(1)) / priorA(2)) ^ 2 - log(priorA(2)) - 0.5 * log(2 * pi);
    end
    if ~isempty(priorB)
      lp(i) = lp(i) - 0.5 * ((bSeq(bi) - priorB(1)) / priorB(2)) ^ 2 - log(priorB(2)) - 0.5 * log(2 * pi);
    end
  end
  post = exp(lp - max(lp)); post = post / sum(post);
end
function [a, b, sHat, seB] = olsFit(x, y)
  x = x(:); y = y(:);
  b = sum((x - mean(x)) .* (y - mean(y))) / sum((x - mean(x)) .^ 2);
  a = mean(y) - b * mean(x);
  res = y - (a + b * x);
  s2 = sum(res .^ 2) / (numel(y) - 2);
  sHat = sqrt(s2); seB = sqrt(s2 / sum((x - mean(x)) .^ 2));
end
function [seq_, p] = marg(aSeq, bSeq, post, which)
  nGrid = numel(aSeq);
  if which == 'a', seq_ = aSeq; idx = mod((1:numel(post))' - 1, nGrid) + 1;
  else,            seq_ = bSeq; idx = floor(((1:numel(post))' - 1) / nGrid) + 1;
  end
  p = zeros(nGrid, 1);
  for i = 1:numel(post), p(idx(i)) = p(idx(i)) + post(i); end
  p = p / sum(p);
end
function [lo, hi] = interval89(seq_, p)
  c_ = cumsum(p);
  lo = seq_(find(c_ >= 0.055, 1)); hi = seq_(find(c_ >= 0.945, 1));
end

%% 0) 教材对照数据
T = readmatrix(fullfile(RES, '00_教材数据_Howell1成年人.csv'));
height = T(:, 1); height_c = T(:, 2); weight = T(:, 3);

%% 1) 主线案例
tHours = (1:12)';
residFix = [0.6; -0.9; 1.2; -0.4; 0.8; -1.3; 0.3; 1.1; -0.7; 0.2; -0.5; 1.4];
thickness = 18 + 3.4 * tHours + residFix;
fid = fopen(fullfile(RES, '01_我们的案例_生物膜厚度_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","厚度"\n');
for i = 1:12, fprintf(fid, '%s,%s\n', numFmt(tHours(i)), numFmt(thickness(i))); end
fclose(fid);

%% 2) 网格后验
[~, ~, sB, seB] = olsFit(tHours, thickness);
[aBio, bBio, postBio] = gridPost(tHours, thickness, N_GRID, sB, A_LIM, B_LIM);
[~, ~, sH] = olsFit(height_c, weight);
[aHw, bHw, postHw] = gridPost(height_c, weight, N_GRID, sH, [30 60], [0.4 0.9]);

[~, obB, obS] = olsFit(tHours, thickness);
[~, ohB, ohS] = olsFit(height_c, weight);
[mbSeq, mbP] = marg(aBio, bBio, postBio, 'b');
[mhSeq, mhP] = marg(aHw, bHw, postHw, 'b');
[ivLo1, ivHi1] = interval89(mbSeq, mbP);
[ivLo2, ivHi2] = interval89(mhSeq, mhP);

idxA = mod((1:N_GRID*N_GRID)' - 1, N_GRID) + 1; idxB = floor(((1:N_GRID*N_GRID)' - 1) / N_GRID) + 1;
meanABio = sum(aBio(idxA) .* postBio); meanBBio = sum(bBio(idxB) .* postBio);
meanAHw = sum(aHw(idxA) .* postHw);    meanBHw = sum(bHw(idxB) .* postHw);

fid = fopen(fullfile(RES, '02_后验摘要_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"案例","最小二乘斜率","后验均值斜率","后验均值截距","斜率89下","斜率89上","残差标准差"\n');
fprintf(fid, '我们的案例：厚度~时间,%s,%s,%s,%s,%s,%s\n', numFmt(round(obB,4),4), numFmt(round(meanBBio,4),4), ...
        numFmt(round(meanABio,4),4), numFmt(round(ivLo1,4),4), numFmt(round(ivHi1,4),4), numFmt(round(obS,4),4));
fprintf(fid, '教材对照：体重~身高,%s,%s,%s,%s,%s,%s\n', numFmt(round(ohB,4),4), numFmt(round(meanBHw,4),4), ...
        numFmt(round(meanAHw,4),4), numFmt(round(ivLo2,4),4), numFmt(round(ivHi2,4),4), numFmt(round(ohS,4),4));
fclose(fid);

%% 3) 二项概率 vs 正态的整点区间概率
nBin = 10; ks = (0:nBin)'; sdBin = sqrt(nBin * 0.25);
pk = binopdf(ks, nBin, 0.5);
pn = normcdf(ks + 0.5, nBin / 2, sdBin) - normcdf(ks - 0.5, nBin / 2, sdBin);
fid = fopen(fullfile(RES, '03_二项概率_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"成功数","概率"\n'); for i = 1:numel(ks), fprintf(fid, '%d,%s\n', ks(i), numFmt(pk(i))); end
fclose(fid);
fid = fopen(fullfile(RES, '03_正态的整点区间概率_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"成功数","区间概率"\n'); for i = 1:numel(ks), fprintf(fid, '%d,%s\n', ks(i), numFmt(pn(i))); end
fclose(fid);

%% 4) 先验预测
nPrior = 60;
pa = A_LIM(1) + (A_LIM(2) - A_LIM(1)) * halton(nPrior, 2);
pb = B_LIM(1) + (B_LIM(2) - B_LIM(1)) * halton(nPrior, 3);
fid = fopen(fullfile(RES, '04_先验预测的直线_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"线号","截距","斜率"\n');
for i = 1:nPrior, fprintf(fid, '%d,%s,%s\n', i, numFmt(round(pa(i),6)), numFmt(round(pb(i),6))); end
fclose(fid);

%% 5) 后验预测
nDraw = 2000; u = halton(nDraw, 5); cPost = cumsum(postBio);
rows = zeros(nDraw, 1);
for i = 1:nDraw, rows(i) = find(cPost >= u(i), 1); end
drA = aBio(mod(rows - 1, N_GRID) + 1); drB = bBio(floor((rows - 1) / N_GRID) + 1);
eps = sB * norminv(halton(nDraw, 11));
w = ones(nDraw, 1) / nDraw;
fid = fopen(fullfile(RES, '05_后验预测_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"时间","预测均值","均值89下","均值89上","单次89下","单次89上"\n');
for tt = 1:12
  mu = drA + drB * tt;
  fprintf(fid, '%s,%s,%s,%s,%s,%s\n', numFmt(tt), numFmt(round(mean(mu),3),3), ...
          numFmt(round(wQuant(mu, w, 0.055),3),3), numFmt(round(wQuant(mu, w, 0.945),3),3), ...
          numFmt(round(wQuant(mu + eps, w, 0.055),3),3), numFmt(round(wQuant(mu + eps, w, 0.945),3),3));
end
fclose(fid);

%% 6) 先验敏感性（只动斜率先验）
sensNames = {'平坦先验', '斜率先验 N(3.4, 0.2)', '斜率先验 N(5.0, 1.2)'};
sensPV = {[], [3.4 0.2], [5.0 1.2]};
fid = fopen(fullfile(RES, '06_先验敏感性_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"斜率先验","斜率均值","斜率89下","斜率89上","先验权重"\n');
for s = 1:3
  [aG, bG, pG] = gridPost(tHours, thickness, N_GRID, sB, A_LIM, B_LIM, [], sensPV{s});
  [seqS, pS] = marg(aG, bG, pG, 'b');
  [loS, hiS] = interval89(seqS, pS);
  meanBS = sum(bG(idxB) .* pG);
  if isempty(sensPV{s}), wt = 0; else, wt = seB^2 / (seB^2 + sensPV{s}(2)^2); end
  fprintf(fid, '%s,%s,%s,%s,%s\n', sensNames{s}, numFmt(round(meanBS,4),4), numFmt(round(loS,4),4), ...
          numFmt(round(hiS,4),4), numFmt(round(wt,4),4));
end
fclose(fid);

%% 7) 网格分辨率
resGrid = [20 50 100 200];
fid = fopen(fullfile(RES, '07_网格分辨率_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"网格点数","后验均值斜率","后验密度最高点斜率","最小二乘参照"\n');
for ng = resGrid
  [~, bG2, pG2] = gridPost(tHours, thickness, ng, sB, A_LIM, B_LIM);
  meanB2 = 0; [~, mi2] = max(pG2); modeB2 = 0;
  for i = 1:numel(pG2)
    bi2 = floor((i - 1) / ng) + 1;
    meanB2 = meanB2 + bG2(bi2) * pG2(i);
    if i == mi2, modeB2 = bG2(bi2); end
  end
  fprintf(fid, '%d,%s,%s,%s\n', ng, numFmt(round(meanB2, 4), 4), ...
          numFmt(round(modeB2, 4), 4), numFmt(round(obB, 4), 4));
end
fclose(fid);

%% 图（与 R 版同内容）
f = figure('Position', [100 100 1110 630], 'Color', 'w'); hold on;
for i = 1:nPrior
  plot([0.5 12.5], pa(i) + pb(i) * [0.5 12.5], 'Color', [TEAL 0.20], 'LineWidth', 0.4);
end
plot(tHours, thickness, 'o', 'MarkerSize', 7, 'MarkerFaceColor', RED, 'MarkerEdgeColor', RED);
plot([0.5 12.5], meanABio + meanBBio * [0.5 12.5], 'Color', DARK, 'LineWidth', 1.6);
title('先验给出一个「直线框」，数据把它收成一条窄带'); xlabel('培养时间（小时）'); ylabel('生物膜厚度（微米）');
box on; exportgraphics(f, fullfile(OUT, '08-MATLAB版-先验与后验直线.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 630], 'Color', 'w'); hold on;
bar(ks, pk, 0.55, 'FaceColor', TEAL); plot(ks, pn, '-o', 'Color', RED, 'LineWidth', 1.3, 'MarkerSize', 5);
title('十个小波动加起来，就近似钟形了'); xlabel('成功次数（0 到 10）'); ylabel('概率'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-小波动相加.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 8 份 CSV 已写入\n');
