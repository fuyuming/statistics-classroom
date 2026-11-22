function result = enumerate_tables(rows, remain, i, acc)
% 自编递归：一行分完后扣除列名额，最后核对所有列名额为0。
if nargin < 3
    i = 1;
    acc = [];
end
if i > numel(rows)
    if all(remain == 0)
        result = {acc};
    else
        result = {};
    end
    return
end
result = {};
choices = compositions(rows(i), remain);
for j = 1:numel(choices)
    row = choices{j};
    next = enumerate_tables(rows, remain - row, i + 1, [acc; row]);
    result = [result next]; %#ok<AGROW>
end
end

function result = compositions(total, caps)
% 只允许当前行计数不超过各列剩余名额。
if numel(caps) == 1
    if total <= caps(1)
        result = {total};
    else
        result = {};
    end
    return
end
result = {};
for first = 0:min(total, caps(1))
    tails = compositions(total - first, caps(2:end));
    for j = 1:numel(tails)
        result{end + 1} = [first tails{j}]; %#ok<AGROW>
    end
end
end
