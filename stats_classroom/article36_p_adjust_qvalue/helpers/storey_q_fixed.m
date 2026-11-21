function [pi0,q] = storey_q_fixed(p,lambda)
% 本篇固定lambda、常规FDR估计模式的Storey q自编函数。
if nargin < 2
    lambda = .5;
end
validateattributes(lambda,{'numeric'},{'scalar','real','finite','>=',0,'<',1});
validateattributes(p,{'numeric'},{'vector','nonempty','real','finite','>=',0,'<=',1});
p = p(:);
pi0 = min(1,mean(p>=lambda)/(1-lambda));
assert(pi0>0,'所选lambda右侧无P值，无法使用本篇pi0估计，请核查设置。');
[~,~,bh,~] = p_adjust_helpers(p);
q = pi0*bh;
end
