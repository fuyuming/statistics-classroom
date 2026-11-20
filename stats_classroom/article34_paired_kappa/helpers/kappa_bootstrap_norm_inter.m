function out = kappa_bootstrap_norm_inter(t, alpha)
% boot::norm.inter同款：在正态分位尺度上插值，非R quantile(type=7)。
t=sort(t(isfinite(t))); R=numel(t); out=zeros(size(alpha));
for j=1:numel(alpha)
    rk=(R+1)*alpha(j); k=floor(rk);
    if k<=0, out(j)=t(1); elseif k>=R, out(j)=t(end);
    elseif k==rk, out(j)=t(k);
    else
        q0=norminv(k/(R+1)); q1=norminv((k+1)/(R+1));
        out(j)=t(k)+(norminv(alpha(j))-q0)/(q1-q0)*(t(k+1)-t(k));
    end
end
end
