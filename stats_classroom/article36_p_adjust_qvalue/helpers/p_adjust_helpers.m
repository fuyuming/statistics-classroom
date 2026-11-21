function [bonf,holm,bh,by] = p_adjust_helpers(p)
% 四个返回值依次为Bonferroni、Holm、BH、BY调整P值。
% 完整自编函数：同一族P值的四种调整，输出保持原输入顺序。
validateattributes(p,{'numeric'},{'vector','nonempty','real','finite','>=',0,'<=',1});
p = p(:);
m = numel(p);
[sorted_p,order_p] = sort(p);
rank_p = (1:m)';
bonf = min(1,m*p);
% Holm为前向累计最大值；BH为后向累计最小值。
holm_sorted = min(1,cummax((m-rank_p+1).*sorted_p));
bh_sorted = min(1,flipud(cummin(flipud(m*sorted_p./rank_p))));
holm = zeros(m,1);
bh = zeros(m,1);
holm(order_p) = holm_sorted;
bh(order_p) = bh_sorted;
by = min(1,bh*sum(1./rank_p));
end
