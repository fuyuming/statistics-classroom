function value = kappa_bootstrap_m(v)
% 第34篇 helper：四格计数整数式κ，避免同值浮点比较。
v = double(v(:));
if numel(v) ~= 4 || any(~isfinite(v)) || any(v < 0) || sum(v)==0
    value = NaN; return
end
a=v(1); b=v(2); c=v(3); d=v(4);
den=(a+b)*(b+d)+(c+d)*(a+c);
if den==0, value=NaN; else, value=2*(a*d-b*c)/den; end
end
