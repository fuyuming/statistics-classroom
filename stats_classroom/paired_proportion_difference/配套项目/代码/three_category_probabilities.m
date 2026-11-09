% 三类人数的完整概率表：20人、+1/0/-1三类（UTF-8）。
%% 读取路径与教学参数
here = fileparts(mfilename('fullpath'));
project = fullfile(here, '..');
dataDir = fullfile(project, '数据');
outDir = fullfile(project, '运行结果', 'MATLAB');
if ~exist(outDir, 'dir'), mkdir(outDir); end
n = 20;
theta = -0.45;
p_t = 0.44;
q_plus = (1 - p_t + theta) / 2;
q_zero = p_t;
q_minus = (1 - p_t - theta) / 2;

%% 观察组合：一种安排概率与安排种数
rObs = 1;
tObs = 14;
sObs = 5;
observedArrangements = nchoosek(n, rObs) * nchoosek(n-rObs, tObs);
singleObservation = q_plus^rObs * q_zero^tObs * q_minus^sObs;
observedProbability = observedArrangements * singleObservation;
observed = table(rObs, tObs, sObs, observedArrangements, singleObservation, observedProbability, ...
    'VariableNames', {'n_plus','n_zero','n_minus','arrangements','single_arrangement_probability','probability'});
disp(observed);

%% 双循环列出全部231个三类人数；用nchoosek手算多项式系数
tail = readtable(fullfile(dataDir, 'wang_lower_tail_n20.csv'));
tailKeys = string(tail.n_plus) + "," + string(tail.n_zero) + "," + string(tail.n_minus);
rowCount = (n+1)*(n+2)/2;
values = zeros(rowCount, 7);
k = 0;
for r = 0:n
    for t = 0:(n-r)
        s = n-r-t;
        k = k+1;
        arrangements = nchoosek(n,r) * nchoosek(n-r,t);
        single = q_plus^r * q_zero^t * q_minus^s;
        probability = arrangements * single;
        key = string(r) + "," + string(t) + "," + string(s);
        inTail = any(tailKeys == key);
        values(k,:) = [r,t,s,arrangements,single,probability,inTail];
    end
end
names = {'n_plus','n_zero','n_minus','arrangements','single_arrangement_probability','probability','in_lower_tail'};
tableOut = array2table(values, 'VariableNames', names);
writetable(tableOut, fullfile(outDir, '三类人数的231组概率.csv'), 'Encoding', 'UTF-8');
tailProbability = sum(tableOut.probability(logical(tableOut.in_lower_tail)));
allProbability = sum(tableOut.probability);
summary = table(n, theta, p_t, q_plus, q_zero, q_minus, height(tableOut), ...
    sum(logical(tableOut.in_lower_tail)), sum(~logical(tableOut.in_lower_tail)), ...
    allProbability, tailProbability, allProbability-tailProbability, observedProbability, observedArrangements, ...
    'VariableNames', {'n','theta','p_t','q_plus','q_zero','q_minus','combinations','tail_combinations','other_combinations', ...
    'all_probability','tail_probability','other_probability','observed_probability','observed_arrangements'});
writetable(summary, fullfile(outDir, '三类人数概率汇总.csv'), 'Encoding', 'UTF-8');
disp(summary);

%% 固定Wang下限尾部的候选theta核对：4001点后在峰值邻域局部优化
candidateTheta = [-0.45, -0.44923, -0.44];
cutoff = 0.025;
counts = [tail.n_plus, tail.n_zero, tail.n_minus];
if height(tail) ~= 141 || any(~isfinite(counts), 'all') || any(counts < 0, 'all') || ...
        any(mod(counts, 1) ~= 0, 'all') || any(sum(counts, 2) ~= n)
    error('wang_lower_tail_n20.csv必须恰好有141行，且每行是非负整数并且总和为20');
end
candidateValues = zeros(numel(candidateTheta), 4);
for j = 1:numel(candidateTheta)
    thetaValue = candidateTheta(j);
    grid = linspace(0, 1-abs(thetaValue), 4001);
    values = zeros(size(grid));
    for i = 1:numel(grid)
        values(i) = wangTailProbability(thetaValue, grid(i), tail, n);
    end
    [~, peak] = max(values);
    left = grid(max(1, peak-1));
    right = grid(min(numel(grid), peak+1));
    [pBest, negProbability] = fminbnd(@(p) -wangTailProbability(thetaValue, p, tail, n), left, right, optimset('TolX', 1e-12));
    candidateValues(j,:) = [thetaValue, pBest, -negProbability, cutoff];
end
candidateOut = array2table(candidateValues, 'VariableNames', {'theta','p_t_at_max','max_tail_probability','cutoff'});
writetable(candidateOut, fullfile(outDir, '候选率差的尾部概率.csv'), 'Encoding', 'UTF-8');
disp(candidateOut);

function probability = wangTailProbability(thetaValue, pValue, tail, n)
qPlusValue = (1 - pValue + thetaValue) / 2;
qMinusValue = (1 - pValue - thetaValue) / 2;
if min([qPlusValue, pValue, qMinusValue]) < -1e-12
    probability = NaN;
    return;
end
qPlusValue = max(0, qPlusValue);
pValue = max(0, pValue);
qMinusValue = max(0, qMinusValue);
coefficient = factorial(n) ./ (factorial(tail.n_plus) .* factorial(tail.n_zero) .* factorial(tail.n_minus));
probability = sum(coefficient .* qPlusValue.^tail.n_plus .* pValue.^tail.n_zero .* qMinusValue.^tail.n_minus);
end
