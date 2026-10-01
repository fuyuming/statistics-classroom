% ============================================================
% 精读 10｜计数里的混杂（MATLAB 版，与 R/Python 逐位一致）
% 对应：McElreath 2023 第 10 讲 Counts and Confounds；教材第 11、12 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
% 主要拟合工具＝显式泊松 IRLS（固定 12 次迭代），偏移项 log(体积) 系数固定为 1。
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
function r = r4(x), r = round(x, 4); end
function [b, info] = poisIRLS(X, y, off, nIter)
  if nargin < 4, nIter = 12; end
  b = zeros(size(X, 2), 1);
  for it = 1:nIter
    eta = X * b + off;
    mu = exp(eta);
    W = mu;
    z = eta + (y - mu) ./ mu - off;
    b = (X' * (W .* X)) \ (X' * (W .* z));
  end
  info = X' * (W .* X);
end

treat = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1]';
batch = [1 1 1 1 1 0 0 0 1 1 1 0 0 0 0 0]';
plate = [10 20 10 20 10 20 10 20 10 20 10 20 10 20 10 20]';
ok    = [1 1 1 0 1 1 1 1 1 1 0 1 1 1 1 1]';
y     = [6 13 5 0 7 6 3 5 8 19 0 9 4 8 5 7]';
offs  = log(plate / 10);
n = 16;

fid = fopen(fullfile(RES, '01_计数数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","批次","体积微升","接种成功","菌落数"\n');
for i = 1:n
  fprintf(fid, '%d,%d,%d,%d,%d,%d\n', i, treat(i), batch(i), plate(i), ok(i), y(i));
end
fclose(fid);

XA = [ones(n,1), treat];
XB = [ones(n,1), treat, batch];
[bA, infoA] = poisIRLS(XA, y, offs);
[bB, infoB] = poisIRLS(XB, y, offs);
fid = fopen(fullfile(RES, '02_率比调整前后_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"模型","率比","下界89","上界89","批次系数"\n');
rrA = exp(bA(2)); sA = seLogF(infoA, 2);
fprintf(fid, '只放处理（含体积偏移）,%s,%s,%s,（未估计）\n', numFmt(r4(rrA),4), ...
        numFmt(r4(exp(log(rrA) - 1.598*sA)),4), numFmt(r4(exp(log(rrA) + 1.598*sA)),4));
rrB = exp(bB(2)); sB = seLogF(infoB, 2);
fprintf(fid, '再加批次（混杂）,%s,%s,%s,%s\n', numFmt(r4(rrB),4), ...
        numFmt(r4(exp(log(rrB) - 1.598*sB)),4), numFmt(r4(exp(log(rrB) + 1.598*sB)),4), numFmt(r4(bB(3)),4));
fclose(fid);

muA = exp(XA * bA + offs);
pZeroPois = mean(exp(-muA));
keep = ok == 1;
[bC, ~] = poisIRLS(XB(keep,:), y(keep), offs(keep));
muC = exp(XB(keep,:) * bC + offs(keep));
fid = fopen(fullfile(RES, '03_零的比例_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值","说明"\n');
fprintf(fid, '观测到的零比例,%s,16 个单元里有 2 个计数为 0（都是没接种上的）\n', numFmt(r4(mean(y == 0)),4));
fprintf(fid, '单一泊松预测的零比例（含全部单元）,%s,泊松把零当成「恰好长了 0 个」，于是严重低估零的出现\n', numFmt(r4(pZeroPois),4));
fprintf(fid, '只看接种成功单元的率比,%s,把没接种上的单元排除后重估\n', numFmt(r4(exp(bC(2))),4));
fprintf(fid, '只看接种成功单元：预测的零比例,%s,排除结构性零之后，泊松对「真零」的预测就合理了\n', numFmt(r4(mean(exp(-muC))),4));
fclose(fid);

muB = exp(XB * bB + offs);
pearson = sum((y - muB).^2 ./ muB);
dfRes = n - size(XB, 2);
phi = pearson / dfRes;
sePois = sB; seQuasi = sePois * sqrt(phi);
fid = fopen(fullfile(RES, '04_过离散_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值"\n');
fprintf(fid, 'Pearson 卡方,%s\n', numFmt(r4(pearson),4));
fprintf(fid, '残差自由度,%s\n', numFmt(r4(dfRes),4));
fprintf(fid, '离散系数 φ（卡方/自由度）,%s\n', numFmt(r4(phi),4));
fprintf(fid, '泊松假设下的 log 率比标准误,%s\n', numFmt(r4(sePois),4));
fprintf(fid, '按 φ 校正后的标准误,%s\n', numFmt(r4(seQuasi),4));
fprintf(fid, '校正倍数,%s\n', numFmt(r4(sqrt(phi)),4));
pearsonFail = sum((y(~keep) - muB(~keep)).^2 ./ muB(~keep));
pearsonShare = pearsonFail / pearson;
phiSucc = sum((y(keep) - muC).^2 ./ muC) / (sum(keep) - size(XB, 2));
fprintf(fid, '两个失败单元贡献的 Pearson 量,%s\n', numFmt(r4(pearsonFail),4));
fprintf(fid, '其占 Pearson 总量的比例,%s\n', numFmt(r4(pearsonShare),4));
fprintf(fid, '只看成功子集的离散系数 φ（14 个观测、3 个参数）,%s\n', numFmt(r4(phiSucc),4));
fclose(fid);

fid = fopen(fullfile(RES, '05_拟合值_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","批次","体积微升","接种成功","菌落数","模型B预测均值"\n');
for i = 1:n
  fprintf(fid, '%d,%d,%d,%d,%d,%d,%s\n', i, treat(i), batch(i), plate(i), ok(i), y(i), numFmt(r4(muB(i)),4));
end
fclose(fid);

bAdj = log(rrB);
fid = fopen(fullfile(RES, '06_校正后的率比_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值"\n');
fprintf(fid, '调整批次后的率比,%s\n', numFmt(r4(exp(bAdj)),4));
fprintf(fid, '泊松假设下的 89%% 区间下,%s\n', numFmt(r4(exp(bAdj - 1.598*sePois)),4));
fprintf(fid, '泊松假设下的 89%% 区间上,%s\n', numFmt(r4(exp(bAdj + 1.598*sePois)),4));
fprintf(fid, '按 φ 校正后的 89%% 区间下,%s\n', numFmt(r4(exp(bAdj - 1.598*seQuasi)),4));
fprintf(fid, '按 φ 校正后的 89%% 区间上,%s\n', numFmt(r4(exp(bAdj + 1.598*seQuasi)),4));
fclose(fid);

%% 图
f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
scatter(plate(batch==0), y(batch==0), 45, [0.086 0.490 0.502], 'filled');
scatter(plate(batch==1), y(batch==1), 45, [0.757 0.275 0.173], 'filled');
title('先看清三件事：体积、批次、以及两个空皿'); xlabel('接种体积（微升）'); ylabel('每皿菌落数');
legend({'批次 0','批次 1'}, 'Location', 'northwest'); box on;
exportgraphics(f, fullfile(OUT, '10-MATLAB版-数据总览.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
vals = [rrA rrB]; lo = [exp(log(rrA)-1.598*sA) exp(log(rrB)-1.598*sB)];
hi = [exp(log(rrA)+1.598*sA) exp(log(rrB)+1.598*sB)];
b1 = bar(1:2, vals, 0.5);
b1.FaceColor = 'flat'; b1.CData = [0.910 0.706 0.659; 0.086 0.490 0.502];
errorbar(1:2, vals, vals - lo, hi - vals, 'LineStyle', 'none', 'Color', [0.086 0.239 0.282], 'CapSize', 8);
yline(1, '--', 'Color', [0.604 0.659 0.659]);
xticks(1:2); xticklabels({'只放处理', '再加批次'});
title('不调整批次，率比就偏了'); ylabel('率比（处理 vs 对照）'); box on;
exportgraphics(f, fullfile(OUT, '11-MATLAB版-率比调整.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 6 份 CSV 已写入\n');
function s = seLogF(info, idx)
  V = inv(info);
  s = sqrt(V(idx, idx));
end
