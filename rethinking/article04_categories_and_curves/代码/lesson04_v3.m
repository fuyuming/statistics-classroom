%% 精读04：分类与样条，MATLAB二次近似复算；基础MATLAB，无额外工具箱。
% Run本文件：自动生成并保留六张中文图，也保存PNG。数据随项目提供。
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'运行结果','v3');
if ~exist(out,'dir'),mkdir(out);end
D=readtable(fullfile(root,'数据','v3','Howell1_adults.csv'));h=D.height;y=D.weight;hc=mean(h);sex=D.male;
G=double([sex==0,sex==1]);X2=[G,G.*(h-hc)];
[m1,v1]=category_fit(G,y);[m2,v3]=category_fit(X2,y);
ch=readtable(fullfile(root,'数据','v3','age_height.csv'));
B=readmatrix(fullfile(root,'数据','v3','age_basis.csv'),'NumHeaderLines',1);pred=readmatrix(fullfile(root,'数据','v3','age_basis_prediction.csv'),'NumHeaderLines',1);
year=pred(:,1);BP=pred(:,2:end);X=[ones(size(B,1),1),B];Y=ch.height;k=size(X,2);
% 样条设计矩阵由作者R的bs定义，三个语言共享，避免结点约定不同。
pmean=[120;zeros(k-1,1)];prec=[1;(1/625)*ones(k-1,1)];
l=fminbnd(@(l) profile(l,X,Y,pmean,prec),-1,4,optimset('TolX',1e-10));
b=beta_at(l,X,Y,pmean,prec);r=X*b-Y;v=exp(-2*l);
H=zeros(k+1);H(1:k,1:k)=v*(X'*X)+diag(prec);
H(1:k,end)=-2*v*(X'*r);H(end,1:k)=H(1:k,end)';H(end,end)=2*v*(r'*r)+4;
ms=[b;l];vs=H\eye(k+1);
writematrix([m1,sqrt(diag(v1))],fullfile(out,'MATLAB_m1.csv'));
writematrix([m2,sqrt(diag(v3))],fullfile(out,'MATLAB_m2.csv'));
writematrix([ms,sqrt(diag(vs))],fullfile(out,'MATLAB_ms.csv'));
disp('各组均值和同身高模型的参数：');disp(m1);disp(m2);
addpath(fileparts(mfilename('fullpath')));
lesson04_plots_v3(root,D,m1,v1,m2,v3,ch,year,BP,ms,vs);
function [t,V]=category_fit(X,y)
k=size(X,2);n=length(y);
if k==2, init=[42,49,5.5];else,init=[45,45,.65,.61,4.2];end
options=optimset('TolX',1e-9,'TolFun',1e-9,'MaxFunEvals',30000,'MaxIter',15000);
t=fminsearch(@(t) objective(t,X,y),init,options)';
b=t(1:k);s=t(end);r=X*b-y;
if k==2, p=[.01;.01];else,p=[.01;.01;0;0];end
H=zeros(k+1);H(1:k,1:k)=X'*X/s^2+diag(p);
H(1:k,end)=-2*(X'*r)/s^3;H(end,1:k)=H(1:k,end)';H(end,end)=-n/s^2+3*(r'*r)/s^4;
V=H\eye(k+1);
end
function v=objective(t,X,y)
k=size(X,2);b=t(1:k)';s=t(end);
if s<=0 || s>=10 || (k==4 && any(b(3:4)<=0 | b(3:4)>=1)),v=1e100;return;end
r=X*b-y;v=length(y)*log(s)+(r'*r)/(2*s^2)+sum((b(1:2)-60).^2)/200;
% Uniform(0,1)先验：边界内为常数。
end
function b=beta_at(l,X,y,p,prec)
v=exp(-2*l);b=(v*(X'*X)+diag(prec))\(v*(X'*y)+prec.*p);
end
function v=profile(l,X,y,p,prec)
b=beta_at(l,X,y,p,prec);r=X*b-y;
v=length(y)*l+(r'*r)*exp(-2*l)/2+sum(prec.*(b-p).^2)/2+2*l*l;
end
