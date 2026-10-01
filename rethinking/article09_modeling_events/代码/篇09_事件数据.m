% ============================================================
% 精读 09｜事件数据（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 09 讲 Modeling Events；教材第 10、11 章
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
function q = wQuant(v, w, p)
  [v2, o] = sort(v); w2 = w(o);
  cw = cumsum(w2) / sum(w2);
  k = find(cw >= p, 1);
  if k == 1, q = v2(1); return; end
  q = v2(k-1) + (v2(k) - v2(k-1)) * (p - cw(k-1)) / (cw(k) - cw(k-1));
end
function [Ag, Bg, W] = grid2(aSeq, bSeq, loglikF)
  [A, B] = meshgrid(aSeq, bSeq);
  Ag = A'; Bg = B';              % 与 R 的 expand.grid 同序：a 变得最快
  Ag = Ag(:); Bg = Bg(:);
  lp = zeros(numel(Ag), 1);
  for i = 1:numel(Ag), lp(i) = loglikF(Ag(i), Bg(i)); end
  w = exp(lp - max(lp));
  W = w / sum(w);
end
function r = r4(x), r = round(x, 4); end

%% 1) 0/1 结果：定植
dose = (1:8)';
nUnits = repmat(10, 8, 1);
nSucc = [1 2 4 5 7 8 9 10]';
logLikB = @(a, b) sum(nSucc .* log(max(min(1 ./ (1 + exp(-(a + b * dose))), 1 - 1e-12), 1e-12)) + ...
                     (nUnits - nSucc) .* log(max(min(1 - 1 ./ (1 + exp(-(a + b * dose))), 1 - 1e-12), 1e-12)));
fid = fopen(fullfile(RES, '01_定植数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","单元数","成功数","成功比例"\n');
for i = 1:8
  fprintf(fid, '%d,%d,%d,%s\n', dose(i), nUnits(i), nSucc(i), numFmt(r4(nSucc(i)/nUnits(i)),4));
end
fclose(fid);
[Ab, Bb, Wb] = grid2(linspace(-6, 2, 300), linspace(0, 1.2, 300), logLikB);
fid = fopen(fullfile(RES, '02_二项后验_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"参数","后验均值","下界89","上界89"\n');
fprintf(fid, '截距 a,%s,%s,%s\n', numFmt(r4(sum(Ab.*Wb)),4), numFmt(r4(wQuant(Ab, Wb, 0.055)),4), numFmt(r4(wQuant(Ab, Wb, 0.945)),4));
fprintf(fid, '斜率 b,%s,%s,%s\n', numFmt(r4(sum(Bb.*Wb)),4), numFmt(r4(wQuant(Bb, Wb, 0.055)),4), numFmt(r4(wQuant(Bb, Wb, 0.945)),4));
fclose(fid);
doseFine = linspace(0, 9, 60);
fid = fopen(fullfile(RES, '03_定植概率曲线_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","均值","下","上"\n');
for d = doseFine
  pp = 1 ./ (1 + exp(-(Ab + Bb * d)));
  fprintf(fid, '%s,%s,%s,%s\n', numFmt(round(d,3),3), numFmt(round(sum(pp.*Wb),3),3), numFmt(round(wQuant(pp, Wb, 0.055),3),3), numFmt(round(wQuant(pp, Wb, 0.945),3),3));
end
fclose(fid);

%% 2) 计数：泊松
treat = repelem((0:7)', 3);
counts = [1 2 3 1 3 5 2 5 8 2 6 10 4 9 14 4 11 18 5 14 23 6 18 30]';
meanCnt = mean(counts);
within = zeros(8, 2);
for d = 0:7, v = counts(treat == d); within(d+1, :) = [var(v) mean(v)]; end
varCnt = sum(within(:,1)) / sum(within(:,2));
logLikP = @(a, b) sum(counts .* log(exp(a + b * treat)) - exp(a + b * treat) - gammaln(counts + 1));
fid = fopen(fullfile(RES, '04_计数数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"处理水平","菌落数"\n');
for i = 1:numel(treat), fprintf(fid, '%d,%d\n', treat(i), counts(i)); end
fclose(fid);
[Ap, Bp, Wp] = grid2(linspace(-1, 3, 300), linspace(0, 0.6, 300), logLikP);
fid = fopen(fullfile(RES, '05_泊松后验_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"参数","后验均值","下界89","上界89"\n');
fprintf(fid, '截距 a,%s,%s,%s\n', numFmt(r4(sum(Ap.*Wp)),4), numFmt(r4(wQuant(Ap, Wp, 0.055)),4), numFmt(r4(wQuant(Ap, Wp, 0.945)),4));
fprintf(fid, '斜率 b,%s,%s,%s\n', numFmt(r4(sum(Bp.*Wp)),4), numFmt(r4(wQuant(Bp, Wp, 0.055)),4), numFmt(r4(wQuant(Bp, Wp, 0.945)),4));
fclose(fid);
rr = exp(Bp * 7);
fid = fopen(fullfile(RES, '06_率比与过离散_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值"\n');
fprintf(fid, '每皿平均菌落数（观测）,%s\n', numFmt(r4(meanCnt),4));
fprintf(fid, '方差/均值（同一剂量内合并，过离散检查）,%s\n', numFmt(r4(varCnt),4));
fprintf(fid, '率比 exp(7b) 后验均值,%s\n', numFmt(r4(sum(rr.*Wp)),4));
fprintf(fid, '率比 89%% 区间下,%s\n', numFmt(r4(wQuant(rr, Wp, 0.055)),4));
fprintf(fid, '率比 89%% 区间上,%s\n', numFmt(r4(wQuant(rr, Wp, 0.945)),4));
fclose(fid);

%% 3) offset
plate = [10 10 10 10]; obsSmall = [3 5 4 4];
plateBig = [40 40 40 40]; obsBig = [12 18 15 19];
naive = mean(obsBig) / mean(obsSmall);
offCorrect = (sum(obsBig)/sum(plateBig)) / (sum(obsSmall)/sum(plate));
fid = fopen(fullfile(RES, '07_offset对照_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"比较","比值","说明"\n');
fprintf(fid, '直接比计数（忽略体积）,%s,大体积皿当然菌落更多，这是体积的功劳\n', numFmt(r4(naive),4));
fprintf(fid, '按体积折算成率再比,%s,把体积放进模型（对数偏移），比的是每微升的密度\n', numFmt(r4(offCorrect),4));
fclose(fid);

%% 4) 有序类别（切点固定，只估斜率）
ordX = [1 1 2 2 3 3 4 4 5 5 6 6]';
ordY = [0 0 0 1 0 1 1 1 1 2 2 2]';
c1Fix = 1.5; c2Fix = 0;
fid = fopen(fullfile(RES, '08_有序数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","等级"\n');
for i = 1:12, fprintf(fid, '%d,%d\n', ordX(i), ordY(i)); end
fclose(fid);
bGrid = linspace(-0.5, 3, 400)';
lpOrd = zeros(400, 1);
for j = 1:400
  b = bGrid(j);
  p0 = min(max(1 ./ (1 + exp(-(c1Fix - b * ordX))), 1e-9), 1 - 1e-9);
  p1 = min(max(1 ./ (1 + exp(-(c2Fix - b * ordX))), 1e-9), 1 - 1e-9);
  ll = 0;
  for i = 1:12
    if ordY(i) == 0, ll = ll + log(p0(i));
    elseif ordY(i) == 1, ll = ll + log(max(p1(i) - p0(i), 1e-9));
    else, ll = ll + log(1 - p1(i)); end
  end
  lpOrd(j) = ll;
end
wOrd = exp(lpOrd - max(lpOrd)); wOrd = wOrd / sum(wOrd);
fid = fopen(fullfile(RES, '09_有序后验_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"参数","后验均值","下界89","上界89"\n');
fprintf(fid, '斜率 b,%s,%s,%s\n', numFmt(r4(sum(bGrid.*wOrd)),4), numFmt(r4(wQuant(bGrid, wOrd, 0.055)),4), numFmt(r4(wQuant(bGrid, wOrd, 0.945)),4));
fclose(fid);
contCoef = polyfit(ordX, ordY, 1);
fid = fopen(fullfile(RES, '10_连续vs有序_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","轻的累积概率","中的累积概率","当连续变量_预测等级"\n');
for d = 1:6
  pk = sum((1 ./ (1 + exp(-(c1Fix - bGrid * d)))) .* wOrd);
  pz = sum((1 ./ (1 + exp(-(c2Fix - bGrid * d)))) .* wOrd);
  fprintf(fid, '%d,%s,%s,%s\n', d, numFmt(r4(pk),4), numFmt(r4(pz),4), numFmt(r4(polyval(contCoef, d)),4));
end
fclose(fid);

%% 图
f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
scatter(dose, nSucc ./ nUnits, 45, [0.757 0.275 0.173], 'filled');
pp = zeros(60,1); lo = zeros(60,1); hi = zeros(60,1);
for i = 1:60
  pv = 1 ./ (1 + exp(-(Ab + Bb * doseFine(i))));
  pp(i) = sum(pv .* Wb); lo(i) = wQuant(pv, Wb, 0.055); hi(i) = wQuant(pv, Wb, 0.945);
end
plot(doseFine, pp, 'Color', [0.086 0.490 0.502], 'LineWidth', 1.5);
plot(doseFine, lo, '--', 'Color', [0.086 0.490 0.502]);
plot(doseFine, hi, '--', 'Color', [0.086 0.490 0.502]);
title('0/1 结果要的是概率曲线，不是直线'); xlabel('剂量'); ylabel('成功概率'); ylim([0 1]); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-logit.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
vals = zeros(6,3);
for d = 1:6
  vals(d,1) = sum((1 ./ (1 + exp(-(c1Fix - bGrid * d)))) .* wOrd);
  vals(d,2) = sum((1 ./ (1 + exp(-(c2Fix - bGrid * d)))) .* wOrd);
  vals(d,3) = polyval(contCoef, d);
end
bar(1:6, vals(:,3), 0.5, 'FaceColor', [0.910 0.706 0.659]);
plot(1:6, vals(:,1), '-o', 'Color', [0.086 0.490 0.502], 'LineWidth', 1.4);
plot(1:6, vals(:,2), '-o', 'Color', [0.757 0.275 0.173], 'LineWidth', 1.4);
title('两种读法给出不同的东西'); xlabel('剂量'); ylabel('值（概率或折算后的等级）');
legend({'把等级当连续变量','累积 logit：P(等级 ≤ 轻)','累积 logit：P(等级 ≤ 中)'}, 'Location', 'northeast'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-连续vs有序.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 10 份 CSV 已写入\n');
