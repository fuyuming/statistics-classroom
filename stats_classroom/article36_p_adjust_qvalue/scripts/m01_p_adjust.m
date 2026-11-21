% 统计课堂36：P值的四种调整，读取同一份1000项模拟P值。
% 在项目根目录，先addpath('scripts')，再运行m01_p_adjust或在Editor逐节运行。
% 共享输入为1000个教学P值：800项无差异、200项有差异；标签只用于核对。

%% 01.1 读取完整比较族
assert(isfile('data/screening_pvalues.csv'), '请将Current Folder设为项目根目录。');
raw = readtable('data/screening_pvalues.csv', 'TextType', 'string');
raw(1:6,:)
p = raw.p;
m = numel(p);
alpha = .05;
addpath(fullfile(pwd,'helpers'));

%% 01.2 手算Bonferroni与Holm
[sorted_p, order_p] = sort(p);
rank_p = (1:m)';
bonf_manual = min(1,m*p);
% Holm走到某一项时，前面的检验都要通过，所以取前向累计最大值。
holm_candidate = (m-rank_p+1).*sorted_p;
holm_sorted = min(1,cummax(holm_candidate));
holm_manual = zeros(m,1);
holm_manual(order_p) = holm_sorted;
holm_steps = table(rank_p,sorted_p,holm_candidate,holm_sorted);
holm_steps(1:6,:)

%% 01.3 手算BH与BY
bh_candidate = m*sorted_p./rank_p;
% flipud倒序，cummin取累计最小值，最后再倒回来。
% 后方位置通过时，前方项一起入选，因此需要查看整个后段。
bh_sorted = min(1,flipud(cummin(flipud(bh_candidate))));
bh_manual = zeros(m,1);
bh_manual(order_p) = bh_sorted;
harmonic_m = sum(1./rank_p);
by_manual = min(1,harmonic_m*bh_manual);
bh_steps = table(rank_p,sorted_p,bh_candidate,bh_sorted);
bh_steps(1:6,:)

%% 01.4 完整函数核对
% p_adjust_helpers是随项目提供的完整自编函数，不是MATLAB原生函数。
[bonf_function,holm_function,bh_function,by_function] = p_adjust_helpers(p);
errors = [max(abs(bonf_function-bonf_manual)); max(abs(holm_function-holm_manual)); ...
          max(abs(bh_function-bh_manual)); max(abs(by_function-by_manual))];
function_check = table(["bonf";"holm";"BH";"BY"],errors, ...
                       'VariableNames',{'method','max_abs_error'});
function_check
assert(max(errors)<1e-10);
% mafdr属于Bioinformatics Toolbox，本节仅用其BHFDR模式核对BH。
if exist('mafdr','file') == 2
    bh_mafdr = mafdr(p,'BHFDR',true);
    native_bh_error = max(abs(bh_mafdr-bh_manual));
    assert(native_bh_error<1e-10);
end

%% 01.5 两个各自单独构成比较族的三项Holm小例
% 每组三项单独构成m=3的比较族，不改变前面1000项的m。
mini_p_stop = [.020;.023;.040];
mini_p_increasing = [.020;.040;.100];
mini_rank = (1:3)';
mini_stop_candidate = (3-mini_rank+1).*mini_p_stop; % 乘数依次为3、2、1。
mini_stop_holm = cummax(mini_stop_candidate);        % 前向累计最大值。
mini_increasing_candidate = (3-mini_rank+1).*mini_p_increasing;
mini_increasing_holm = cummax(mini_increasing_candidate);
mini_example = table( ...
    ["stop";"stop";"stop";"increasing";"increasing";"increasing"], ...
    [mini_p_stop;mini_p_increasing], ...
    [mini_stop_candidate;mini_increasing_candidate], ...
    [mini_stop_holm;mini_increasing_holm], ...
    'VariableNames',{'example','p','holm_candidate','holm'});
mini_example
disp('stop例第一关未通过，后续步骤停止；调整值仍返回完整三项。');
mini_expected = [.060;.060;.060;.060;.080;.100]; % 仅作核对，不作为计算列。
[~, mini_holm_stop, ~, ~] = p_adjust_helpers(mini_p_stop);
[~, mini_holm_increasing, ~, ~] = p_adjust_helpers(mini_p_increasing);
assert(max(abs(mini_example.holm - mini_expected)) < 1e-10);
assert(max(abs(mini_example.holm - [mini_holm_stop; mini_holm_increasing])) < 1e-10);

%% 01.6 原题三疗法复算：每个2x2比较使用Yates连续性校正
% therapy前缀避免污染前面1000项的p、m和调整结果。
therapy_raw = readtable('data/therapies.csv', 'TextType', 'string');
therapy_table = therapy_raw{:, {'effective','ineffective'}};
therapy_pairs = [1 2; 1 3; 2 3];
therapy_names = ["physical-drug"; "physical-topical"; "drug-topical"];
therapy_raw_p = zeros(3,1);
therapy_expected_min = zeros(3,1);
for therapy_i = 1:3
    % true：明确开启2×2卡方检验的Yates连续性校正。
    [~, therapy_raw_p(therapy_i), ~, therapy_expected] = ...
        therapy_chisq_helpers(therapy_table(therapy_pairs(therapy_i,:),:), true);
    therapy_expected_min(therapy_i) = min(therapy_expected(:));
end
therapy_bonf = min(1, 3 * therapy_raw_p);
[therapy_sorted, therapy_order] = sort(therapy_raw_p);
therapy_holm_sorted = min(1, cummax((3-(1:3)'+1).*therapy_sorted));
therapy_holm = zeros(3,1);
therapy_holm(therapy_order) = therapy_holm_sorted;
therapy_results = table(therapy_names, therapy_raw_p, therapy_bonf, therapy_holm, ...
    therapy_expected_min, 'VariableNames', ...
    {'pair','raw_p_yates','bonferroni_p','holm_p','expected_min'});
therapy_results
assert(max(abs(therapy_raw_p - [.016754306570385496; .000009344615515985061; ...
    .047763774036015214])) < 1e-10);
assert(sum(therapy_bonf <= .05) == 1 && sum(therapy_holm <= .05) == 3);

%% 01.7 查看完整结果与发现数
results = raw;
results.bonf = bonf_manual;
results.holm = holm_manual;
results.BH = bh_manual;
results.BY = by_manual;
counts = [sum(p<=alpha);sum(bonf_manual<=alpha);sum(holm_manual<=alpha); ...
          sum(bh_manual<=alpha);sum(by_manual<=alpha)];
summary_result = table(["raw";"bonf";"holm";"BH";"BY"],counts, ...
                       'VariableNames',{'method','discoveries'});
summary_result

%% 01.8 可选保存
% 需要保存时先setenv('ARTICLE36_SAVE_RESULTS','1')，再顺序运行。
if strcmp(getenv('ARTICLE36_SAVE_RESULTS'),'1')
    if ~isfolder('results')
        mkdir('results');
    end
    writetable(results,'results/MATLAB_01_p_adjust.csv');
    writetable(summary_result,'results/MATLAB_01_summary.csv');
    writetable(mini_example,'results/MATLAB_01_mini_example.csv');
    writetable(therapy_results,'results/MATLAB_01_therapies.csv');
end
