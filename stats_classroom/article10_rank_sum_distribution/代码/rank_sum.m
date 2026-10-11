%% 统计课堂10：秩和精确分布。只需基础MATLAB，不需统计工具箱。
% UTF-8保存；Run完整脚本。两列人工评分分别演示无并列与并列情形。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir')
    mkdir(out);
end
logfile = fullfile(out,'逐题结果.txt');
if exist(logfile,'file')
    delete(logfile);
end
diary(logfile);
data = readtable(fullfile(root,'数据','scores.csv'));
disp('输入为人工教学数据：group组别；value无并列评分；value_tied并列评分。');
disp(data);
is_a = strcmp(data.group,'A');
m = sum(is_a);
n = sum(~is_a);
N = m+n;
assert(m==3 && n==3 && numel(unique(data.value))==N);
allocations = nchoosek(1:N,m); % 每行选三个记录下标，共20行，不删除重复值记录。
mean_w = m*(N+1)/2;
results = zeros(2,4);
columns = {'value','value_tied'};

%% 1. 平均秩、记录分配、精确概率
for scenario = 1:2
    values = data.(columns{scenario});
    ranks = zeros(N,1);
    for i = 1:N
        % 平均秩 = 比它小的个数 + (与它相同的个数+1)/2。
        ranks(i) = sum(values<values(i))+(sum(values==values(i))+1)/2;
    end
    w_obs = sum(ranks(is_a));
    u_obs = w_obs-m*(m+1)/2;
    all_w = sum(ranks(allocations),2);
    var_exact = mean((all_w-mean(all_w)).^2); % 完整分布，分母20，不是19。
    p_two = mean(abs(all_w-mean_w)>=abs(w_obs-mean_w)-1e-12);
    [~,~,tie_group] = unique(values);
    tied_sizes = accumarray(tie_group,1);
    correction = sum(tied_sizes.^3-tied_sizes);
    var_formula = m*n/12*(N+1-correction/(N*(N-1)));
    fprintf('\n问题：%s；合并秩如下：\n',columns{scenario});
    disp(ranks');
    fprintf('W=%.6f；U=%.6f；均值=%.6f；方差=%.6f；双侧P=%.6f。\n',w_obs,u_obs,mean(all_w),var_exact,p_two);
    disp('解释：以离零假设中心至少同样远定义双侧P；20种记录分配各占1/20。');
    alloc_table = array2table([allocations,all_w],'VariableNames',{'record1','record2','record3','W'});
    writetable(alloc_table,fullfile(out,[columns{scenario},'_allocations.csv']));
    [W,~,index] = unique(all_w);
    count = accumarray(index,1);
    probability = count/numel(all_w);
    pmf = table(W,count,probability);
    disp(pmf);
    writetable(pmf,fullfile(out,[columns{scenario},'_pmf.csv']));
    assert(abs(var_exact-var_formula)<1e-12);
    assert(abs(mean(all_w)-mean_w)<1e-12);
    results(scenario,:) = [w_obs,u_obs,var_exact,p_two];
end

%% 2. 连续性校正近似：标准正态CDF用erfc计算，无需统计工具箱
z_left = (results(1,1)+.5-mean_w)/sqrt(results(1,3));
p_normal = erfc(-z_left/sqrt(2)); % 2*Phi(z_left)，本例左侧观察值。
fprintf('本例连续性校正双侧正态近似=%.6f；精确值=0.2。\n',p_normal);
name = {'W';'U';'mean_W';'var_W';'p_exact';'var_tied';'p_tied';'p_normal_cc'};
value = [results(1,1);results(1,2);mean_w;results(1,3);results(1,4);results(2,3);results(2,4);p_normal];
writetable(table(name,value),fullfile(out,'theory.csv'));

%% 3. 两组各20，反复不放回抽名次
rng(20261011);
m_large = 20;
n_large = 20;
N_large = m_large+n_large;
repeats = 20000;
w_sim = zeros(repeats,1);
for b = 1:repeats
    selected = randperm(N_large,m_large);
    w_sim(b) = sum(selected);
end
center = m_large*(N_large+1)/2;
sd_w = sqrt(m_large*n_large*(N_large+1)/12);
z_sim = (w_sim-center)/sd_w;
writetable(table(w_sim,z_sim,'VariableNames',{'W','Z'}),fullfile(out,'large_simulation.csv'));
fprintf('模拟：W中心=%.6f；SD=%.6f；Z模拟均值=%.6f、方差=%.6f。\n',center,sd_w,mean(z_sim),var(z_sim));
disp('解释：模拟帮助观察正态近似，不是证明；随机结果不要求跨语言逐行一致。');
figure('Color','w');
histogram(z_sim,-8:.2:8,'Normalization','pdf','FaceColor',[.65,.8,.77]);
hold on;
x = linspace(-4.5,4.5,501);
plot(x,exp(-x.^2/2)/sqrt(2*pi),'LineWidth',2);
xlim([-4.5,4.5]);
xlabel('Standardized W');
ylabel('Density');
title('Rank-sum null simulation: m=n=20');
exportgraphics(gcf,fullfile(out,'normal_check.png'),'Resolution',160);
diary off;
