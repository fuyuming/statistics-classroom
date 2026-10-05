%% 精读03 v5：同一个模型的网格与二次近似
% 在 MATLAB 打开本文件并 Run。仅使用基础 MATLAB，不需要统计工具箱。
% MATLAB 实现二次近似原理；真正的 rethinking::quap 见同目录 R 脚本。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','v5');
if ~exist(out,'dir'), mkdir(out); end
D = readtable(fullfile(root,'运行结果','00_教材数据_Howell1成年人.csv'));
disp(D(1:5,:));
x = D.height-mean(D.height);
y = D.weight;
n = length(y);
sy = sum(y); sx = sum(x); sxx = sum(x.^2); sxy = sum(x.*y); syy = sum(y.^2);
%% 三个参数都未知：先验与作者相同
% a~Normal(60,10), b~LogNormal(0,1), sigma~Uniform(0,10)。
% 下列范围只是积分窗口，扩大窗口检查见 Python 主脚本。
[A,B,S] = ndgrid(linspace(43,47,101),linspace(.45,.8,101),linspace(3.3,5.5,101));
SSE = syy-2*A*sy-2*B*sxy+n*A.^2+2*A.*B*sx+B.^2*sxx;
logw = -n*log(S)-SSE./(2*S.^2)-(A-60).^2/200-log(B)-log(B).^2/2;
w = exp(logw-max(logw(:)));
w = w/sum(w(:)); % 每组权重除以所有权重之和。
gridmean = [sum(w(:).*A(:));sum(w(:).*B(:));sum(w(:).*S(:))];
%% 二次近似：先找峰值，再计算峰附近的弯曲程度
objective = @(t) neglog(t,n,sy,sx,sxx,sxy,syy);
options = optimset('TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',10000);
t = fminsearch(objective,[45,.63,4.2],options);
a=t(1); b=t(2); s=t(3);
sse = syy-2*a*sy-2*b*sxy+n*a*a+2*a*b*sx+b*b*sxx;
H = [n/s^2+.01,sx/s^2,-2*(n*a+b*sx-sy)/s^3;...
 sx/s^2,sxx/s^2-log(b)/b^2,-2*(a*sx+b*sxx-sxy)/s^3;...
 -2*(n*a+b*sx-sy)/s^3,-2*(a*sx+b*sxx-sxy)/s^3,-n/s^2+3*sse/s^4];
V = H\eye(3);
parameter = ["alpha_centered";"beta";"sigma"];
result = table(parameter,gridmean,t',sqrt(diag(V)),'VariableNames',{'parameter','grid_mean','quadratic_mean','quadratic_sd'});
writetable(result,fullfile(out,'MATLAB_comparison.csv'));
disp(result);
disp('斜率描述群体关联；sigma 也有后验分布，没有被固定。');
% 练习：把每维101改成151，检查结果是否稳定，注意组合数随三次方增加。
function v=neglog(t,n,sy,sx,sxx,sxy,syy)
 a=t(1); b=t(2); s=t(3);
 if b<=0 || s<=0 || s>=10, v=1e100; return; end
 sse=syy-2*a*sy-2*b*sxy+n*a*a+2*a*b*sx+b*b*sxx;
 v=n*log(s)+sse/(2*s*s)+(a-60)^2/200+log(b)+log(b)^2/2;
end
