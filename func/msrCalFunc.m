function MSR = msrCalFunc(sig,N)
% 计算信号平均谱半径
% sig：时间序列信号
% N：特征值数量

%% 信号预处理
sig_len = length(sig);
T = floor(sig_len/N);
X_NT = reshape(sig(1:N*T),N,T);

%% X_NT按行归一化
X_NT_mean = repmat(mean(X_NT,2),1,T); %按行计算均值
X_NT_std = repmat(std(X_NT,0,2),1,T); %按行计算标准差
X_NT_norm = (X_NT-X_NT_mean)./X_NT_std;

%% 生成Haar矩阵U
Z = randn(N);
U = Z*(Z'*Z)^(-0.5);

%% 计算矩阵X_NT_norm的奇异值矩阵等价形式X_u
X_u = ((X_NT_norm*(X_NT_norm'))^0.5)*U; 

%% X_u按行归一化
X_u_std = repmat(std(X_u,0,2),1,N); 
X_u_norm = X_u./(N.^0.5)./X_u_std;

%% 求解特征值
X_u_norm_eig=eig(X_u_norm);

% % 特征值在复平面的分布图
% figure;scatter(real(X_u_norm_eig), imag(X_u_norm_eig), 'filled');
% xlabel('Real');
% ylabel('Imaginary');
% xticks(-1:0.5:1);
% xlim([-1,1]);
% yticks(-1:0.5:1);
% ylim([-1,1]);


%% 求解平均谱半径
MSR = sum(abs(X_u_norm_eig))/N;

















