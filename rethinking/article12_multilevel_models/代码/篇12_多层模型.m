% ============================================================
% 精读 12｜多层模型（MATLAB 版，与 R/Python 逐位一致）
% 对应：McElreath 2023 第 12 讲 Multilevel Models；教材第 13 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
% 模型 y_ij = α_j + b·x_ij + ε，α_j ~ Normal(μ, σ_α)；对 (b, σ_α, σ) 铺三维网格，
% 组截距与 μ 解析积分掉（闭式边际似然）。
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
function q = wQuant(v, w, p)
  [vs, o] = sort(v(:)); ws = w(o);
  cw = cumsum(ws) / sum(ws);
  k = find(cw >= p, 1);
  if isempty(k), q = vs(end); return; end
  if k == 1, q = vs(1); return; end
  q = vs(k-1) + (vs(k) - vs(k-1)) * (p - cw(k-1)) / (cw(k) - cw(k-1));
end
function g = shrinkFcn(sa, s, n_j)
  g = sa.^2 ./ (sa.^2 + s.^2 / n_j);
end
function ll = margLL(b, sa, s, gmean, gxbar, wy, wx, n_j, J, N)
  r = gmean - b * gxbar;
  v = repmat(s^2 / n_j + sa^2, J, 1);
  S = sum(1 ./ v);
  muHat = sum(r ./ v) / S;
  q = sum((r - muHat).^2 ./ v);
  w = wy - b * wx;
  ll = -0.5 * (log(S) + sum(log(v)) + q) - 0.5 * sum(w.^2) / s^2 - (N - J) * log(s);
end

x = repmat([1; 2; 3], 12, 1);
alphaTrue = [4.60; 4.90; 5.15; 5.35; 5.50; 5.60; 5.70; 5.85; 6.05; 6.30; 6.60; 7.00];
resid = [-0.66; 0.11; 0.55;  0.55; -0.66; 0.11; -0.26; 0.78; -0.52;
          0.39; -0.65; 0.26; -0.66; 0.11; 0.55;  0.55; -0.66; 0.11;
         -0.26; 0.78; -0.52; 0.39; -0.65; 0.26; -0.66; 0.11; 0.55;
          0.55; -0.66; 0.11; -0.26; 0.78; -0.52; 0.39; -0.65; 0.26];
gid = repelem((1:12)', 3);
y = alphaTrue(gid) + 1.0 * x + resid;
n_j = 3; J = 12; N = 36;

fid = fopen(fullfile(RES, '01_数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"瓶","剂量","长势"\n');
for i = 1:N, fprintf(fid, '%d,%d,%s\n', gid(i), x(i), numFmt(y(i))); end   % 长势与 R/Python 同为 6 位小数
fclose(fid);

gmean = zeros(J,1);
for j = 1:J, gmean(j) = mean(y(gid == j)); end
gxbar = mean(x); gbar = mean(y);
wy = zeros(N,1); wx = zeros(N,1);
for j = 1:J
  m = gid == j;
  wy(m) = y(m) - gmean(j);
  wx(m) = x(m) - mean(x(m));
end

Xp = [ones(N,1), x];
bp = (Xp' * Xp) \ (Xp' * y);
b_fe = sum(wy .* wx) / sum(wx.^2);
unpooled = gmean - b_fe * gxbar;
pooled = repmat(mean(unpooled), J, 1);

bSeq = linspace(0, 2, 60)'; saSeq = linspace(0.005, 2, 60)'; sSeq = linspace(0.05, 2, 60)';
[Bg, SAg, Sg] = ndgrid(bSeq, saSeq, sSeq);
Bv = reshape(Bg, [], 1); SAv = reshape(SAg, [], 1); Sv = reshape(Sg, [], 1);
lp = zeros(numel(Bv), 1);
for i = 1:numel(Bv), lp(i) = margLL(Bv(i), SAv(i), Sv(i), gmean, gxbar, wy, wx, n_j, J, N); end
post = exp(lp - max(lp)); post = post / sum(post);

muDraw = mean(gmean) - Bv * gxbar;
partial = zeros(J,1);
for j = 1:J
  a = muDraw + (gmean(j) - Bv * gxbar - muDraw) .* shrinkFcn(SAv, Sv, n_j);
  partial(j) = sum(a .* post);
end
bPost = sum(Bv.*post); muPost = sum(muDraw.*post);
saPost = sum(SAv.*post); sPost = sum(Sv.*post);
shrinkVec = shrinkFcn(SAv, Sv, n_j);
shrinkPost = sum(shrinkVec .* post);

fid = fopen(fullfile(RES, '02_三种算法_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"瓶","瓶内均值","完全合并","完全不合并","部分合并","收缩量"\n');
for j = 1:J
  fprintf(fid, '%d,%s,%s,%s,%s,%s\n', j, numFmt(r4(unpooled(j)),4), numFmt(r4(pooled(j)),4), ...
          numFmt(r4(unpooled(j)),4), numFmt(r4(partial(j)),4), numFmt(r4(unpooled(j)-partial(j)),4));
end
fclose(fid);

fid = fopen(fullfile(RES, '03_超参数后验_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"参数","后验均值","下界89","上界89"\n');
fprintf(fid, '共同斜率 b,%s,%s,%s\n', numFmt(r4(bPost),4), numFmt(r4(wQuant(Bv,post,0.055)),4), numFmt(r4(wQuant(Bv,post,0.945)),4));
fprintf(fid, '总体水平 μ,%s,NA,NA\n', numFmt(r4(muPost),4));
fprintf(fid, '组间标准差 σ_α,%s,%s,%s\n', numFmt(r4(saPost),4), numFmt(r4(wQuant(SAv,post,0.055)),4), numFmt(r4(wQuant(SAv,post,0.945)),4));
fprintf(fid, '残差标准差 σ,%s,%s,%s\n', numFmt(r4(sPost),4), numFmt(r4(wQuant(Sv,post,0.055)),4), numFmt(r4(wQuant(Sv,post,0.945)),4));
fprintf(fid, '收缩因子（每瓶都是 3 个孔）,%s,NA,NA\n', numFmt(r4(shrinkPost),4));
fclose(fid);

seExisting = sqrt(sPost^2 + sPost^2/n_j);
seNew = sqrt(sPost^2 + saPost^2);
fid = fopen(fullfile(RES, '04_两种不确定性_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"情形","预测标准差","区间半宽89","说明"\n');
fprintf(fid, '同一个瓶再测一个孔,%s,%s,只多了观测误差与瓶均值自身的不确定性\n', numFmt(r4(seExisting),4), numFmt(r4(1.598*seExisting),4));
fprintf(fid, '换一个全新的瓶再测一个孔,%s,%s,还要加上瓶与瓶之间的真实差异 σ_α\n', numFmt(r4(seNew),4), numFmt(r4(1.598*seNew),4));
fclose(fid);

rr = corrcoef(unpooled - mean(unpooled), unpooled - partial);
fid = fopen(fullfile(RES, '05_汇总_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值"\n');
fprintf(fid, '完全合并的水平,%s\n', numFmt(r4(mean(unpooled)),4));
fprintf(fid, '完全不合并的极差,%s\n', numFmt(r4(max(unpooled)-min(unpooled)),4));
fprintf(fid, '部分合并的极差,%s\n', numFmt(r4(max(partial)-min(partial)),4));
fprintf(fid, '收缩量与偏离的相关,%s\n', numFmt(r4(rr(1,2)),4));
fprintf(fid, '瓶间真实差异（生成时设定）,%s\n', numFmt(r4(std(alphaTrue)),4));
fprintf(fid, '瓶内噪声（生成时设定）,%s\n', numFmt(r4(sqrt(mean(resid.^2))),4));
fclose(fid);

%% 两张对照图
f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
for j = 1:J, scatter(repmat(j,3,1), y(gid==j), 25, [0.086 0.490 0.502], 'filled', 'MarkerFaceAlpha', 0.85); end
scatter(1:J, gmean, 55, [0.757 0.275 0.173], 'filled');
yline(gbar, '--', 'Color', [0.604 0.659 0.659]);
title('12 个培养瓶，每瓶 3 个孔（MATLAB 版对照图）'); xlabel('培养瓶'); ylabel('长势'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-数据.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1050 750], 'Color', 'w'); hold on;
plot([min(unpooled) max(unpooled)], [min(unpooled) max(unpooled)], '--', 'Color', [0.604 0.659 0.659]);
yline(mean(unpooled), ':', 'Color', [0.086 0.239 0.282]);
scatter(unpooled, partial, 45, [0.086 0.490 0.502], 'filled');
for j = 1:J
  quiver(unpooled(j), unpooled(j), 0, partial(j)-unpooled(j), 0, 'Color', [0.757 0.275 0.173], 'LineWidth', 1.1);
end
title('部分合并＝把每个瓶的估计往整体拉一点（MATLAB 版）');
xlabel('完全不合并（瓶内均值）'); ylabel('部分合并'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-收缩.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 5 份 CSV 已写入\n');
