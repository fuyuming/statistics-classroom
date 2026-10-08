%% 精读05 v2：MATLAB打开本文件并Run，无额外工具箱。
% 与PPT31页相同模型；用fminsearch寻找后验峰值，解析Hessian作二次近似。
% 数据一行一个地区。两张中文Figure保留，并保存PNG。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','v2','MATLAB');
if ~exist(out,'dir'), mkdir(out); end
raw = readmatrix(fullfile(root,'数据','v2','WaffleDivorce.csv'));
disp(raw(1:5,:));
s = (raw-mean(raw))./std(raw,0,1); % 样本标准差，分母n-1
Y = s(:,1);
X = [ones(size(Y)),s(:,2:3)];
ps = [.2;.5;.5];
% 参数顺序：截距alpha、结婚率斜率beta_M、结婚年龄斜率beta_A、sigma。
% ps为前三个参数的正态先验标准差；sigma另使用Exponential(1)先验。
% 本节对应正文婚姻模型：在年龄相同的比较中估计结婚率斜率。
objective = @(t) loss(t,X,Y,ps);
options = optimset('MaxFunEvals',20000,'MaxIter',10000,'TolX',1e-10,'TolFun',1e-10);
[t,~,flag] = fminsearch(objective,[0;0;-.5;.8],options);
assert(flag>0,'后验优化未收敛');
r = Y-X*t(1:3); sig=t(4);
H=zeros(4);
H(1:3,1:3)=X'*X/sig^2+diag(1./ps.^2);
H(1:3,4)=2*X'*r/sig^3;
H(4,1:3)=H(1:3,4)';
H(4,4)=-length(Y)/sig^2+3*sum(r.^2)/sig^4;
% 逆Hessian保存联合后验的近似协方差，峰值t用作近似均值。
% 这是贝叶斯二次近似，未使用MCMC，也没有把sigma当已知常数。
cv=H\eye(4); sd=sqrt(diag(cv));
q=sqrt(2)*erfinv(2*.945-1);
result=[t,sd,t-q*sd,t+q*sd];
writetable(array2table(result,'VariableNames',{'mean','sd','lower89','upper89'}),fullfile(out,'posterior.csv'));
disp(result);
%% 精确枚举作者的二元规则
% m=1叉，2管，3对撞，4中间变量后代；前三种A为独立占位变量。
% 每行是一个可能组合，概率沿箭头相乘，单个模型的总概率应为1。
rows=[];
for m=1:4
 for x=0:1
  for z=0:1
   for y=0:1
    for a=0:1
     if m==1
      w=.5*bern(x,.1+.8*z)*bern(y,.1+.8*z)*.5;
     elseif m==2
      w=.5*bern(z,.1+.8*x)*bern(y,.1+.8*z)*.5;
     elseif m==3
      p=.2+.7*(x+y>0); w=.25*bern(z,p)*.5;
     else
      w=.5*bern(z,.1+.8*x)*bern(y,.1+.8*z)*bern(a,.1+.8*z);
     end
     rows=[rows;m,x,z,y,a,w];
    end
   end
  end
 end
end
answers=[];
for m=1:4
 for g=0:4
  sub=rows(rows(:,1)==m,:);
  if g>0
   col=3; if g>=3, col=5; end
   sub=sub(sub(:,col)==mod(g-1,2),:);
  end
  % 筛选后重新归一化，得到组内条件概率，再算相关。
  w=sub(:,6)/sum(sub(:,6));
  ex=sum(w.*sub(:,2)); ey=sum(w.*sub(:,4));
  corr=(sum(w.*sub(:,2).*sub(:,4))-ex*ey)/sqrt(ex*(1-ex)*ey*(1-ey));
  answers=[answers;m,g,corr];
 end
end
writetable(array2table(rows,'VariableNames',{'model','X','Z','Y','A','probability'}),fullfile(out,'joint.csv'));
writetable(array2table(answers,'VariableNames',{'model','group','correlation'}),fullfile(out,'associations.csv'));
%% 原书171页：植物理论平均增长；无需随机抽样。
plant=[0,.5,(1-.5)*5+.5*2;1,.1,(1-.1)*5+.1*2];
writetable(array2table(plant,'VariableNames',{'treatment','fungus_probability','mean_growth'}),fullfile(out,'plant_expectation.csv'));
disp(plant);
%% 中文绘图
available=listfonts;
candidates={'PingFang SC','Microsoft YaHei','Noto Sans CJK SC','Heiti SC'};
font='';
for k=1:numel(candidates)
 if any(strcmp(available,candidates{k})), font=candidates{k}; break; end
end
assert(~isempty(font),'请安装Noto Sans CJK SC中文字体');
f=figure('Color','w');
errorbar(t(2:3),[2;1],q*sd(2:3),'horizontal','o');
xline(0,'--'); ylim([.5,2.5]); yticks([1,2]); yticklabels({'结婚年龄','结婚率'});
xlabel('标准化斜率及89%后验区间'); title('年龄相同以后，结婚率的斜率接近零');
set(findall(f,'-property','FontName'),'FontName',font);
exportgraphics(f,fullfile(out,'model.png'),'Resolution',180);
f=figure('Color','w');
vals=answers(answers(:,1)==4 & ismember(answers(:,2),[0,1,3]),3);
bar(vals); xticks(1:3); xticklabels({'全体','固定中间变量Z','固定后代A'}); ylim([0,.75]);
ylabel('理论相关'); title('后代携带中间变量的信息');
set(findall(f,'-property','FontName'),'FontName',font);
exportgraphics(f,fullfile(out,'descendant.png'),'Resolution',180);
fprintf('两张Figure已生成并保留。\n');
function value=loss(t,X,y,ps)
 if t(4)<=0, value=Inf; return; end
 r=y-X*t(1:3);
 % 依次为正态似然两项、正态先验项、速率1的指数先验项；常数已略。
 value=length(y)*log(t(4))+sum(r.^2)/(2*t(4)^2)+sum((t(1:3)./ps).^2)/2+t(4);
end
function value=bern(x,p)
 if x==1, value=p; else, value=1-p; end
end
