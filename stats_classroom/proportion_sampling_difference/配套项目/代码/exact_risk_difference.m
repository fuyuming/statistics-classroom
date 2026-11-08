function result = exact_risk_difference(a,n1,b,n2,conf_level,grid_points)
% 两个独立率之差的无条件精确推断（UTF-8，基础MATLAB，不调用R）。
% 无参数运行本文件，演示第一组2/10减第二组0/10的教学例题。
% 也可调用 result=exact_risk_difference(2,10,0,10,0.95,501)。
% 按样本率差排序，二项概率相乘再求和，对未知总体率作数值最大化。
% central双侧P为两倍较小单侧P，区间反演两侧各2.5%的检验。
% 对应R扩展的整数排序共同核对项，详见README及exact2x2官方说明。
demo = nargin == 0;
if demo
    a = 2;
    n1 = 10;
    b = 0;
    n2 = 10;
end
if nargin < 5
    conf_level = 0.95;
end
if nargin < 6
    grid_points = 501;
end
validateattributes(a,{'numeric'},{'scalar','integer','finite','nonnegative'});
validateattributes(b,{'numeric'},{'scalar','integer','finite','nonnegative'});
validateattributes(n1,{'numeric'},{'scalar','integer','finite','positive'});
validateattributes(n2,{'numeric'},{'scalar','integer','finite','positive'});
assert(a <= n1 && b <= n2,'阳性人数不能超过样本量。');
validateattributes(conf_level,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(grid_points,{'numeric'},{'scalar','integer','finite','>=',3});
observed = a/n1 - b/n2;
tail_alpha = (1-conf_level)/2;
% 交叉乘积是整数，相同率差保留为并列，避免浮点相减产生误差。
first_counts = (0:n1)';
second_counts = 0:n2;
count_difference = first_counts*n2 - second_counts*n1;
observed_integer = a*n2 - b*n1;
upper_mask = double(count_difference >= observed_integer);
lower_mask = double(count_difference <= observed_integer);
upper_p = maximum_tail(0,n1,n2,upper_mask,grid_points);
lower_p = maximum_tail(0,n1,n2,lower_mask,grid_points);
result.method = "unconditional_difference_central";
result.estimate = observed;
result.p_value = min(1,2*min(upper_p,lower_p));
if observed == -1
    result.lower = -1;
else
    result.lower = endpoint(-1,observed,true,n1,n2,upper_mask,grid_points,tail_alpha);
end
if observed == 1
    result.upper = 1;
else
    result.upper = endpoint(observed,1,false,n1,n2,lower_mask,grid_points,tail_alpha);
end
result.grid_points = grid_points;
if demo
    fprintf('第一组2/10，第二组0/10，计算第一组减第二组。\n');
    fprintf('率差=%.6f，双侧P=%.8f\n',result.estimate,result.p_value);
    fprintf('95%%区间=[%.8f, %.8f]\n',result.lower,result.upper);
    fprintf('区间包含0，在双侧alpha=0.05下，差异证据不足。\n');
    fprintf('本方法按率差排序，对应R扩展的整数排序核对项，详见README。\n');
    fprintf('可把grid_points从501改为1001，比较数值搜索末位变化。\n');
    code_folder = fileparts(mfilename('fullpath'));
    output = fullfile(fileparts(code_folder),'运行结果','MATLAB');
    if ~exist(output,'dir')
        mkdir(output);
    end
    writetable(struct2table(result),fullfile(output,'率差精确_三语言共同方法.csv'));
end
end

function best = maximum_tail(delta,n1,n2,mask,grid_points)
% 给定delta=p1-p2，令q=p2，两个率均在[0,1]确定q的合法区间。
left = max(0,-delta);
right = min(1,1-delta);
if right <= left
    best = tail_probability(left,delta,n1,n2,mask);
    return
end
grid = linspace(left,right,grid_points)';
values = tail_probability(grid,delta,n1,n2,mask);
[best,highest] = max(values);
if best >= 1-1e-14
    best = 1;
    return
end
peaks = find(values(2:end-1) >= values(1:end-2) ...
    & values(2:end-1) >= values(3:end) ...
    & (values(2:end-1) > values(1:end-2) | values(2:end-1) > values(3:end)))+1;
peaks = unique([peaks;highest]);
options = optimset('TolX',1e-12,'Display','off');
for index = peaks'
    if index > 1 && index < grid_points
        objective = @(q) -tail_probability(q,delta,n1,n2,mask);
        [~,negative_max,exitflag] = fminbnd(objective,grid(index-1),grid(index+1),options);
        assert(exitflag > 0,'概率最大化未收敛，请检查输入或增加网格。');
        best = max(best,-negative_max);
    end
end
end

function probability = tail_probability(q,delta,n1,n2,mask)
q = q(:);
p1 = min(1,max(0,q+delta));
p2 = min(1,max(0,q));
first = binomial_probabilities(n1,p1);
second = binomial_probabilities(n2,p2);
% 每行是一个q，mask选择尾部计数组合，联合概率依据独立性相乘。
probability = sum((first*mask).*second,2);
probability = min(1,max(0,probability));
end

function probability = binomial_probabilities(n,p)
% gammaln避免阶乘溢出，p=0和p=1单独赋值，避免0乘log(0)。
k = 0:n;
probability = zeros(numel(p),n+1);
inside = p > 0 & p < 1;
log_choose = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1);
probability(inside,:) = exp(log_choose + log(p(inside))*k + log1p(-p(inside))*(n-k));
probability(p==0,1) = 1;
probability(p==1,end) = 1;
end

function value = endpoint(left,right,upper_tail,n1,n2,mask,grid_points,alpha)
% 下限右尾P随delta递增，上限左尾P随delta递减。
for iteration = 1:60
    middle = (left+right)/2;
    probability = maximum_tail(middle,n1,n2,mask,grid_points);
    if (probability < alpha) == upper_tail
        left = middle;
    else
        right = middle;
    end
    if right-left <= 1e-7
        break
    end
end
% 取外侧括号，避免末位舍入把区间缩窄。
if upper_tail
    value = left;
else
    value = right;
end
end
