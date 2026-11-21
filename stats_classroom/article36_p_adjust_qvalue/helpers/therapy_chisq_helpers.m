function [statistic,p_value,df,expected] = therapy_chisq_helpers(A,correct)
% Base MATLAB Pearson chi-square helper with optional Yates correction.
row_total = sum(A,2);
column_total = sum(A,1);
expected = row_total * column_total / sum(row_total);
delta = abs(A - expected);
if correct && all(size(A)==[2 2])
    delta = max(delta - .5, 0);
end
statistic = sum((delta.^2) ./ expected, 'all');
df = (size(A,1)-1) * (size(A,2)-1);
p_value = gammainc(statistic/2, df/2, 'upper');
end
