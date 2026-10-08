%% 精读05 v2：MATLAB打开本文件并Run，无额外工具箱。
% 与PPT31页相同模型；用fminsearch寻找后验峰值，解析Hessian作二次近似。
% 数据一行一个地区。两张中文Figure保留，并保存PNG。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','v2','MATLAB');
if ~exist(out,'dir'), mkdir(out); end
raw = readmatrix(fullfile(root,'数据','v2','WaffleDivorce.csv'));
disp(raw(1:5,:));
% raw三列顺序为离婚率、结婚率（每千名成年人）、结婚年龄中位数（岁）。
% 每行一个州或特区。换数据时核对顺序、单位与缺失值，再按列标准化。
s = (raw-mean(raw))./std(raw,0,1); % 样本标准差，分母n-1
Y = s(:,1);
% X是设计矩阵：常数1、标准化结婚率M、标准化结婚年龄A；Y是标准化离婚率。
% 此处矩阵X与下文因果图的二元x含义不同。
X = [ones(size(Y)),s(:,2:3)];
% 先验怎么填：ps依次为alpha、beta_M、beta_A的正态先验标准差，不是方差。
% 三个均值均为0；sigma的指数先验速率为1，写在末尾loss函数中。
ps = [.2;.5;.5];
% 参数顺序：截距alpha、结婚率斜率beta_M、结婚年龄斜率beta_A、sigma。
% ps为前三个参数的正态先验标准差；sigma另使用Exponential(1)先验。
% 本节对应正文婚姻模型：在年龄相同的比较中估计结婚率斜率。
objective = @(t) loss(t,X,Y,ps);
options = optimset('MaxFunEvals',20000,'MaxIter',10000,'TolX',1e-10,'TolFun',1e-10);
% 四个优化初值顺序为alpha、beta_M、beta_A、sigma，最后一个必须为正。
% 初值与先验不同：调整先验标准差改ps，不改这四个初值。
% MaxIter/MaxFunEvals是优化上限，不是后验抽样次数；这里没有MCMC或网格遍历。
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
% 中央89%区间用0.945分位点；erfinv在基础MATLAB中换算标准正态分位数。
% 改95%时使用0.975，并同步列名、图注。
q=sqrt(2)*erfinv(2*.945-1);
result=[t,sd,t-q*sd,t+q*sd];
writetable(array2table(result,'VariableNames',{'mean','sd','lower89','upper89'}),fullfile(out,'posterior.csv'));
% 结果四行是alpha/beta_M/beta_A/sigma，四列是均值、后验标准差、区间下限/上限。
% sigma的均值约0.785是剩余波动大小；该行sd是对sigma的不确定性。
% beta_M约-0.065，区间约[-0.306,0.176]，表示方向仍不确定，不是证明效应为0。
% beta_A约-0.614：结婚率相同时，年龄增加1个标准差，平均离婚率降低0.614个标准差。
disp(result);
%% 精确枚举作者的二元规则
% m=1叉，2管，3对撞，4中间变量后代；前三种A为独立占位变量。
% 每行是一个可能组合，概率沿箭头相乘，单个模型的总概率应为1。
% 0.1+0.8*z把z=0/1映射到概率0.1/0.9；bern返回当前状态的概率。
% 这里a为二元后代，不是结婚年龄；枚举状态不等于在参数网格上计算后验。
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
% g=0全体，1/2分别固定Z=0/1，3/4分别固定A=0/1。
% 练习只改m=4最后的bern(a,.1+.8*z)为bern(a,.5)，其余机制不变。
% A不再包含Z的信息，固定A的相关应从0.390回到全体的0.640。
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
% 三列是处理状态、真菌概率、平均增长；0是不处理，1是处理。
% 无真菌长5，有真菌长2，按概率加权得3.5/4.7，总效应差1.2。
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
% 第一张图：点为斜率近似后验均值，横线为89%区间，竖虚线为零。
% 第二张图：三柱是全体/固定Z/固定A的理论相关，不是后验区间。
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
% 若改指数先验速率lambda，最后+t(4)应改为+lambda*t(4)。
% 正态先验和指数先验之外的分布，需要同时重推上面的Hessian。
 value=length(y)*log(t(4))+sum(r.^2)/(2*t(4)^2)+sum((t(1:3)./ps).^2)/2+t(4);
end
function value=bern(x,p)
 if x==1, value=p; else, value=1-p; end
end
