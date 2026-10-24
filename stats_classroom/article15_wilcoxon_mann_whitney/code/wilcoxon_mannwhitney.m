%% 第一部分：看数据→明问题→立假设；支持项目根或 code/ 工作目录。
if isfile(fullfile('data','textbook_example.csv'))
    data_dir = 'data';
elseif isfile(fullfile('..','data','textbook_example.csv'))
    data_dir = fullfile('..','data');
else
    error('找不到共享数据：textbook_example.csv');
end
textbook = readtable(fullfile(data_dir,'textbook_example.csv'));
rd_data = readtable(fullfile(data_dir,'independent_RD.csv'));

% 教材例5-8：甲、乙各5个独立小池，单位万尾；预定方向是甲高于乙。
% H0：两组总体分布相同；前提：两组独立、数据可排序；无并列值。
A = textbook.value(textbook.group == "A");
B = textbook.value(textbook.group == "B");
m = numel(A);
n = numel(B);

%% 第二部分：选方法→查条件；合并编秩，解释原始秩和与U。
% 本例原始秩和R甲=37、R乙=18，扣除组内秩和15后U甲=22、U乙=3。
pooled_rank = tiedrank([A;B]);
RA = sum(pooled_rank(1:m));
RB = sum(pooled_rank(m+1:end));
UA = RA - m*(m+1)/2;
UB = RB - n*(n+1)/2;
win_A = sum(A > B', 'all');
win_B = sum(A < B', 'all');
win_tie = sum(A == B', 'all');
effect_A = UA/(m*n);
textbook_summary = table(m,n,RA,RB,UA,UB,win_A,win_B,win_tie,effect_A, ...
    'VariableNames',{'n_A','n_B','R_A','R_B','U_A','U_B','wins_A','wins_B','ties','effect_A'});
textbook_summary

%% 第三部分：做检验→判显著→看效应；枚举10选5的252个标签分配，U甲≥22有7种。
N = m+n;
assignment = nchoosek(1:N,m);
u_all = sum(assignment,2) - m*(m+1)/2;
count_assignments = numel(u_all);
p_upper = mean(u_all >= UA);
p_lower = mean(u_all <= UA);
p_two_sided = min(1,2*min(p_upper,p_lower));
u_critical_lower = 4;
p_critical_lower = mean(u_all <= u_critical_lower);
lower_count = sum(u_all <= u_critical_lower);

% HL=5；共同形状位置平移模型下，U临界值4对应第(u临界值+1)=5个有序差值2，实际覆盖20/21=.952381。
cross_diff = sort(reshape(A - B',[],1));
HL = median(cross_diff);
lower_index = u_critical_lower + 1;
lower_95_nominal = cross_diff(lower_index);
actual_confidence = 1 - p_critical_lower;
exact_summary = table(count_assignments,sum(u_all >= UA),p_upper,p_lower, ...
    p_two_sided,u_critical_lower,lower_count,p_critical_lower,HL,lower_95_nominal,Inf, ...
    actual_confidence,'VariableNames',{'assignments','upper_count','p_upper', ...
    'p_lower','p_two_sided','u_critical_lower','lower_count','p_critical_lower','HL', ...
    'ci_lower','ci_upper','actual_confidence'});
exact_summary

% MATLAB ranksum返回原始秩和37；A>B的U由秩和换算为22。
[p_exact,~,matlab_stats] = ranksum(A,B,'tail','right','method','exact');
matlab_summary = table(matlab_stats.ranksum,p_exact, ...
    'VariableNames',{'R_A','p_exact'});
matlab_summary

% 交换组回归检查：交换输入后改用B<A的等价方向，单侧P应为0.0277778。
[swap_p_less,~,swap_stats] = ranksum(B,A,'tail','left','method','exact');
swap_u = swap_stats.ranksum - n*(n+1)/2;
swap_p_upper = mean(u_all >= swap_u);
swap_p_lower = mean(u_all <= swap_u);
swap_p_two_sided = min(1,2*min(swap_p_upper,swap_p_lower));
swap_summary = table(swap_u,swap_p_less,swap_p_two_sided, ...
    'VariableNames',{'swapped_U','swapped_p_less','swapped_p_two_sided'});
swap_summary

%% 第四部分：RD课件对照；两组n=10、12，有并列值，秩取平均，近似+连续性校正。
RD1 = rd_data.RD(rd_data.group == 1);
RD2 = rd_data.RD(rd_data.group == 2);
rd_rank = tiedrank([RD1;RD2]);
rd_n1 = numel(RD1);
rd_n2 = numel(RD2);
RA_RD = sum(rd_rank(1:rd_n1));
RB_RD = sum(rd_rank(rd_n1+1:end));
UA_RD = RA_RD - rd_n1*(rd_n1+1)/2;
UB_RD = RB_RD - rd_n2*(rd_n2+1)/2;
rd_win_tie = sum(RD1 == RD2', 'all');
rd_effect = UA_RD/(rd_n1*rd_n2);
[rd_p,~,rd_stats] = ranksum(RD1,RD2,'tail','right','method','approximate');
rd_summary = table(rd_n1,rd_n2,RA_RD,UA_RD,UB_RD,rd_win_tie,rd_effect,rd_p, ...
    'VariableNames',{'n_1','n_2','R_1','U_1','U_2','ties','effect_1','p_approx_continuity'});
rd_summary

%% 第五部分：配对资料；同一问题可写成差值的单样本符号秩或paired调用。
% 来源为颜艳、王彤教材课件的8对大鼠维生素A数据，X与Y按pair对应。
% H0：差值分布关于0对称；d=round(X-Y,2)，零差剔除；本例n=8、无绝对差并列。
paired_data = readtable(fullfile(data_dir,'paired_vitaminA.csv'));
X = paired_data.X;
Y = paired_data.Y;
d = round(X-Y,2);
d_nonzero = d(d ~= 0);
d_rank = tiedrank(abs(d_nonzero));
positive_rank_sum = sum(d_rank(d_nonzero > 0));
negative_rank_sum = sum(d_rank(d_nonzero < 0));
signed_rank_summary = table(numel(d_nonzero),positive_rank_sum,negative_rank_sum, ...
    'VariableNames',{'n','positive_rank_sum','negative_rank_sum'});
signed_rank_summary
[p_one,~,stats_one] = signrank(d,0,'tail','both','method','exact');
[p_pair,~,stats_pair] = signrank(X,Y,'tail','both','method','exact');
paired_summary = table(p_one,p_pair,stats_one.signedrank,stats_pair.signedrank, ...
    'VariableNames',{'one_sample_p','paired_p','one_sample_V','paired_V'});
paired_summary

%% 可选第六部分：保存可分享的关键结果CSV；默认不写文件。
result_export = table("MATLAB",RA,RB,UA,UB,p_upper,p_two_sided,HL, ...
    lower_95_nominal,Inf,actual_confidence,RA_RD,UA_RD,UB_RD,rd_effect,rd_p, ...
    'VariableNames',{'language','R_A','R_B','U_A','U_B','p_upper', ...
    'p_two_sided','HL','ci_lower','ci_upper','actual_confidence','RD_R1', ...
    'RD_U1','RD_U2','RD_effect','RD_p'});
result_export
save_csv = false;
if save_csv
    output_dir = fullfile(data_dir,'..','results');
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    writetable(result_export,fullfile(output_dir,'MATLAB_results.csv'));
end
