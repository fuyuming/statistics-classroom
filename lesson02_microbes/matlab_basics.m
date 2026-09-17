%% MATLAB基础练习：类型、结构、索引、条件、循环、函数
% double是常见数值类型；string表示文本；logical表示真假。
weight = 2.3;
count = int32(7);
label = "control";
is_valid = true;
disp(class(weight)); disp(class(count)); disp(class(label)); disp(class(is_valid));
%% 数组和数据容器
x = [1.8 2.0 NaN 2.4];
A = [1 2; 3 4];
C = {"control", [1.8 2.0]}; % cell允许不同内容
S.group = "control"; S.unit = "g"; % struct使用字段
T = table(["C01"; "C02"], [1.8; 2.0], 'VariableNames', {'id','weight'});
disp(x(1)); disp(A(1,2)); disp(C{2}); disp(S.unit); disp(T.weight);
%% 条件判断：numel计数；缺失数值用isnan识别
valid_x = x(~isnan(x));
if numel(valid_x) >= 2
    disp(std(valid_x, 0));
else
    disp("有效观测不足两个");
end
%% 循环：依次取下标
for k = 1:numel(valid_x)
    disp(valid_x(k) * 1000); % 克转成毫克
end
%% 函数：脚本末尾定义局部函数；完整运行此文件
result = mean_valid([1.8 2.0 NaN]);
disp(result); % 1.9
assert(abs(result - 1.9) < 1e-12);
function result = mean_valid(x)
    x = x(~isnan(x));
    if isempty(x)
        result = NaN;
    else
        result = mean(x);
    end
end
