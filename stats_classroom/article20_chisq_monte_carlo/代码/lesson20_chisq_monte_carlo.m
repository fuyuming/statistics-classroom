%% 统计课堂20：ACE基因型2×3表的条件Monte Carlo
% 输入为课堂249人汇总表。MATLAB打开本文件后Run完整运行。
% 按脚本自身位置找数据；基础MATLAB即可，无需Statistics工具箱。
% 输出到运行结果/MATLAB，按下面一至六节顺序计算。

%% 一、原题、期望表与Pearson统计量
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir')
    mkdir(out);
end
input = fullfile(root,'数据','ace_genotype.csv');
data = readtable(input,'TextType','string');
disp(data);
O = data{:,{'DD','ID','II'}};
rows = sum(O,2);
cols = sum(O,1);
N = sum(O,'all');
% 列向量rows乘行向量cols，得到每个行合计与列合计的乘积。
E = rows*cols/N;
contribution = (O-E).^2./E;
Q = sum(contribution,'all');
% gammainc的upper选项给χ²右尾；本表不应用Yates校正。
asymptoticP = gammainc(Q/2,2/2,'upper');
threshold = Q*(1-64*eps);
disp(E);
disp(contribution);
disp(Q);
% 转置后展开，按DN三格、noDN三格统一保存。
group = repelem(data.group,3);
genotype = repmat(["DD";"ID";"II"],2,1);
cells = table(group,genotype,reshape(O.',[],1),reshape(E.',[],1), ...
    reshape(contribution.',[],1), ...
    'VariableNames',{'group','genotype','observed','expected','pearson_contribution'});
writetable(cells,fullfile(out,'cells.csv'));

%% 二、用卡片标签无放回抽111人，生成一张表
labels = [ones(1,cols(1)),2*ones(1,cols(2)),3*ones(1,cols(3))];
rng(123,'twister');
picked = randperm(N,rows(1));
% 标签1、2、3分别落入以0.5、1.5、2.5、3.5为边界的三个箱。
first_row = histcounts(labels(picked),0.5:1:3.5);
first = [first_row;cols-first_row];
disp(first);
fprintf('单张模拟表的χ²：%.9f\n',sum((first-E).^2./E,'all'));
one = table(data.group,first(:,1),first(:,2),first(:,3), ...
    'VariableNames',{'group','DD','ID','II'});
writetable(one,fullfile(out,'one_table.csv'));

%% 三、重复B次，边际保持不变，计算右尾次数与P
B = 100000;
rng(123,'twister');
qs = zeros(B,1);
first_rows = zeros(B,3);
for i = 1:B
    picked = randperm(N,rows(1));
    current_first = histcounts(labels(picked),0.5:1:3.5);
    first_rows(i,:) = current_first;
    current = [current_first;cols-current_first];
    qs(i) = sum((current-E).^2./E,'all');
end
second_rows = cols-first_rows;
assert(all(sum(first_rows,2)==rows(1)));
assert(all(sum(second_rows,2)==rows(2)));
assert(all(first_rows+second_rows==cols,'all'));
assert(all(first_rows>=0,'all') && all(second_rows>=0,'all'));
assert(all(first_rows==floor(first_rows),'all'));
assert(all(qs>=0));
b = sum(qs>=threshold);
P = (b+1)/(B+1);
SE = sqrt(P*(1-P)/B);
writetable(table(qs,'VariableNames',{'Q'}),fullfile(out,'simulation.csv'));
fprintf('极端次数=%d，P=%.10f，模拟SE=%.10f\n',b,P,SE);

%% 四、穷举第一行，按条件概率累加Pearson右尾
% 一旦第一行DD、ID确定，II=111-DD-ID；第二行由列合计确定。
% gammaln(n+1)是log(n!)；组合数先取log再求概率，避免巨大阶乘。
records = zeros((cols(1)+1)*(cols(2)+1),5);
index = 0;
log_denominator = gammaln(N+1)-gammaln(rows(1)+1)-gammaln(rows(2)+1);
for dd = 0:cols(1)
    for id_count = 0:cols(2)
        ii = rows(1)-dd-id_count;
        if ii>=0 && ii<=cols(3)
            a = [dd,id_count,ii];
            current = [a;cols-a];
            current_q = sum((current-E).^2./E,'all');
            log_numerator = sum(gammaln(cols+1)-gammaln(a+1)-gammaln(cols-a+1));
            probability = exp(log_numerator-log_denominator);
            index = index+1;
            records(index,:) = [a,current_q,probability];
        end
    end
end
records = records(1:index,:);
en = array2table(records,'VariableNames',{'DD','ID','II','Q','probability'});
en.pearson_tail = en.Q>=threshold;
exact = sum(en.probability(en.pearson_tail));
assert(height(en)==4063 && abs(sum(en.probability)-1)<1e-10);
assert(abs(P-exact)<5*SE);
writetable(en,fullfile(out,'enumeration.csv'));
fprintf('Pearson精确条件尾概率：%.14f\n',exact);

%% 五、Cramer's V与近似95%区间
% 用非中心χ²累计概率反演；本例用0至100的非中心参数区间求根。
effect_v = sqrt(Q/N);
options = optimset('TolX',1e-12);
low_ncp = fzero(@(value) noncentral_cdf(Q,2,value)-.975,[0,100],options);
high_ncp = fzero(@(value) noncentral_cdf(Q,2,value)-.025,[0,100],options);
ci_low = sqrt(low_ncp/N);
ci_high = sqrt(high_ncp/N);
% builtinP仅R列有数值，用于R手工过程与chisq.test核对。
summary = table(N,Q,2,asymptoticP,B,b,P,SE,exact,effect_v,ci_low,ci_high,NaN, ...
    'VariableNames',{'N','Qobs','df','asymptoticP','B','b','P','SE', ...
    'exact_pearson_P','effect_v','ci_low','ci_high','builtinP'});
writetable(summary,fullfile(out,'summary.csv'));
disp(summary);

%% 六、同一模拟序列的不同前缀：精度练习
sizes = [2000;10000;B];
sizes = unique(sizes(sizes<=B));
exercise = zeros(length(sizes),4);
for i = 1:length(sizes)
    count = sum(qs(1:sizes(i))>=threshold);
    estimate = (count+1)/(sizes(i)+1);
    error = sqrt(estimate*(1-estimate)/sizes(i));
    exercise(i,:) = [sizes(i),count,estimate,error];
end
exercise_table = array2table(exercise,'VariableNames',{'B','b','P','SE'});
writetable(exercise_table,fullfile(out,'precision_exercise.csv'));
disp(exercise_table);
fprintf('MATLAB_CHECKS_PASSED：边际、整数、概率总和及Monte Carlo误差。\n');

function probability = noncentral_cdf(x,df,lambda)
% 非中心χ²的CDF是中心χ² CDF的Poisson加权和。
% 本例求根范围0<=lambda<=100，j=0:200足以使遗漏权重远低于数值容差。
    if lambda==0
        probability = gammainc(x/2,df/2);
        return;
    end
    weight = exp(-lambda/2);
    probability = 0;
    for j = 0:200
        probability = probability+weight*gammainc(x/2,df/2+j);
        weight = weight*(lambda/2)/(j+1);
    end
end
