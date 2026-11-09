function result = paired_exact(nplus,nzero,nminus,conf_level,grid_points)
% 三项配对率差的原生精确推断，方向为后测减前测（基础MATLAB）。
% 运行paired_exact，默认演示1个+1、14个0、5个-1。
% 按D排序的central构造，与R PairedCI的Wang归纳构造分别报告。
demo = nargin == 0;
if demo
    nplus = 1;
    nzero = 14;
    nminus = 5;
end
if nargin < 4
    conf_level = 0.95;
end
if nargin < 5
    grid_points = 201;
end
validateattributes(nplus,{'numeric'},{'scalar','integer','finite','nonnegative'});
validateattributes(nzero,{'numeric'},{'scalar','integer','finite','nonnegative'});
validateattributes(nminus,{'numeric'},{'scalar','integer','finite','nonnegative'});
n = nplus+nzero+nminus;
assert(n > 0,'总人数须大于0。');
validateattributes(conf_level,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(grid_points,{'numeric'},{'scalar','integer','finite','>=',3});
observed = (nplus-nminus)/n;
tail_alpha = (1-conf_level)/2;
% 枚举合法的三类人数，直接比较整数差，保留相同D的所有组合。
cells = zeros((n+1)*(n+2)/2,3);
index = 0;
for r = 0:n
    for m = 0:(n-r)
        index = index+1;
        cells(index,:) = [r,n-r-m,m];
    end
end
count_difference = cells(:,1)-cells(:,3);
observed_integer = nplus-nminus;
upper_mask = double(count_difference >= observed_integer);
lower_mask = double(count_difference <= observed_integer);
log_coefficient = gammaln(n+1)-sum(gammaln(cells+1),2);
upper_p = maximum_tail(0,cells,log_coefficient,upper_mask,grid_points);
lower_p = maximum_tail(0,cells,log_coefficient,lower_mask,grid_points);
result.method = "paired_fixed_D_central";
result.estimate = observed;
result.p_value = min(1,2*min(upper_p,lower_p));
if observed == -1
    result.lower = -1;
else
    result.lower = endpoint(-1,observed,true,cells,log_coefficient,upper_mask,grid_points,tail_alpha);
end
if observed == 1
    result.upper = 1;
else
    result.upper = endpoint(observed,1,false,cells,log_coefficient,lower_mask,grid_points,tail_alpha);
end
result.grid_points = grid_points;
if demo
    fprintf('20人：1个+1、14个0、5个-1，方向为后测减前测。\n');
    fprintf('率差=%.6f，双侧P=%.8f\n',result.estimate,result.p_value);
    fprintf('95%%区间=[%.8f, %.8f]\n',result.lower,result.upper);
    fprintf('区间包含0，在双侧alpha=0.05下，差异证据不足。\n');
    fprintf('本方法按D排序，R PairedCI的归纳构造区间另见README。\n');
    fprintf('可把grid_points从201改为401，比较数值搜索末位变化。\n');
    code_folder = fileparts(mfilename('fullpath'));
    output = fullfile(fileparts(code_folder),'运行结果','MATLAB');
    if ~exist(output,'dir')
        mkdir(output);
    end
    writetable(struct2table(result),fullfile(output,'共同精确区间.csv'));
end
end

function best = maximum_tail(delta,cells,log_coefficient,mask,grid_points)
% 给定delta=p1-p2，令q=p2，两个率均在[0,1]确定q的合法区间。
left = 0;
right = max(0,1-abs(delta));
if right <= left
    best = tail_probability(left,delta,cells,log_coefficient,mask);
    return
end
grid = linspace(left,right,grid_points)';
values = tail_probability(grid,delta,cells,log_coefficient,mask);
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
        objective = @(q) -tail_probability(q,delta,cells,log_coefficient,mask);
        [~,negative_max,exitflag] = fminbnd(objective,grid(index-1),grid(index+1),options);
        assert(exitflag > 0,'概率最大化未收敛，请检查输入或增加网格。');
        best = max(best,-negative_max);
    end
end
end

function probability = tail_probability(q,delta,cells,log_coefficient,mask)
q = q(:);
qp = min(1,max(0,(1-q+delta)/2));
qm = min(1,max(0,(1-q-delta)/2));
probabilities = [qp,q,qm];
log_mass = repmat(log_coefficient',numel(q),1);
for category = 1:3
    p = probabilities(:,category);
    count = cells(:,category)';
    inside = p > 0;
    log_mass(inside,:) = log_mass(inside,:)+log(p(inside))*count;
    log_mass(~inside,count>0) = -Inf;
end
mass = exp(log_mass);
probability = min(1,max(0,mass*mask));
end

function value = endpoint(left,right,upper_tail,cells,log_coefficient,mask,grid_points,alpha)
% 下限右尾P随delta递增，上限左尾P随delta递减。
for iteration = 1:60
    middle = (left+right)/2;
    probability = maximum_tail(middle,cells,log_coefficient,mask,grid_points);
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
