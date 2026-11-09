% Wang排序前三步：仅演示第三位候选的选择（MATLAB基础版）
%% 相对路径与固定教学参数
here = fileparts(mfilename('fullpath'));
out = fullfile(fileparts(here), '运行结果', 'MATLAB');
if ~exist(out, 'dir')
    mkdir(out);
end
n = 20;
alpha_tail = 0.025;
candidates = {'A','B'};
counts = [19 0 1; 18 2 0];

%% 对给定theta、pT计算临时尾概率
tail_probability = @(theta,p_t,candidate) local_tail(theta,p_t,candidate);

%% 全域网格后，在网格峰附近作fminbnd，并比较两个端点
peak_over_pt = @(theta,candidate) local_peak(theta,candidate,tail_probability);

%% 对每个候选单独求maxQ(theta)=.025的最小根
summary = zeros(2,9);
for j = 1:2
    candidate = candidates{j};
    lower = 0;
    upper = 1;
    for i = 1:70
        middle = (lower + upper)/2;
        peak = peak_over_pt(middle,candidate);
        if peak(2) >= alpha_tail
            upper = middle;
        else
            lower = middle;
        end
    end
    peak = peak_over_pt(upper,candidate);
    summary(j,:) = [j counts(j,:) (counts(j,1)-counts(j,3))/n upper peak(1) peak(2) 0];
end

%% 比较计算出的临时下限，相同最大值一起选入
summary(:,9) = abs(summary(:,6) - max(summary(:,6))) <= 1e-10;

%% A候选在pT=0时的解析式复核
a_at_a = (1 + summary(1,6))/2;
analytic_A = a_at_a^20 + 20*a_at_a^19*(1-a_at_a);
disp(analytic_A)

%% 保存两类教学结果
summary_table = table(categorical(candidates'),counts(:,1),counts(:,2),counts(:,3),summary(:,5),summary(:,6),summary(:,7),summary(:,8),logical(summary(:,9)), ...
    'VariableNames',{'candidate','n_plus','n_zero','n_minus','D','trial_lower','peak_p_t','tail_probability','selected_for_third'});
writetable(summary_table,fullfile(out,'Wang排序第三位试算.csv'));
tail_counts = [20 0 0; 19 1 0; counts(1,:); 20 0 0; 19 1 0; counts(2,:)];
tail_candidates = repelem(categorical(candidates'),3);
tail_table = table(tail_candidates,tail_counts(:,1),tail_counts(:,2),tail_counts(:,3), ...
    'VariableNames',{'candidate','n_plus','n_zero','n_minus'});
writetable(tail_table,fullfile(out,'Wang排序候选尾部.csv'));
disp(summary_table);
assert(all(sum(tail_counts,2) == n));
assert(height(tail_table) == 6);
assert(isequal(tail_counts(1:2,:),tail_counts(4:5,:)));
assert(summary(2,9) == 1);
assert(all(abs(summary(:,8)-alpha_tail) < 1e-10));

% 局部函数：三项尾概率
function value = local_tail(theta,p_t,candidate)
    a = (1-p_t+theta)/2;
    b = (1-p_t-theta)/2;
    value = a^20 + 20*a^19*p_t;
    if candidate == 'A'
        value = value + 20*a^19*b;
    else
        value = value + 190*a^18*p_t^2;
    end
end

% 局部函数：网格、局部优化与端点比较
function result = local_peak(theta,candidate,fun)
    right = 1-abs(theta);
    grid = linspace(0,right,401);
    values = arrayfun(@(p) fun(theta,p,candidate),grid);
    [best_q,k] = max(values);
    best_p = grid(k);
    left = grid(max(1,k-1));
    right_local = grid(min(401,k+1));
    opts = optimset('TolX',1e-14);
    [p_opt,~] = fminbnd(@(p) -fun(theta,p,candidate),left,right_local,opts);
    options_p = [0 right p_opt];
    options_q = arrayfun(@(p) fun(theta,p,candidate),options_p);
    [q_local,idx] = max(options_q);
    if q_local > best_q
        best_q = q_local;
        best_p = options_p(idx);
    end
    result = [best_p best_q];
end
