function sig_ECA = ecaFunc(sig_surv,sig_ref,delayDot)
% ECA
%%
sig_surv_len = length(sig_surv);
sig_surv = reshape(sig_surv,sig_surv_len,1);
sig_ref_len = length(sig_ref);
sig_ref = reshape(sig_ref,sig_ref_len,1);


%% 信号预处理 
% 使参考信号与监测信号等长
if sig_surv_len >= sig_ref_len
    sig_ref = [sig_ref;zeros(sig_surv_len-sig_ref_len,1)];
else
    sig_ref = sig_ref(1:sig_surv_len);
end

% 增加时延单元
sig_ref = [zeros(delayDot,1);sig_ref]; % 信号前补0

%% 构造杂波投影子空间
% 初始化
X = zeros(length(sig_ref),delayDot);

% 零多普勒杂波子空间（可另加多普勒频率信息）
for i = 1:delayDot
    X(:,i) =  circshift(sig_ref,i-1);

end
X = X(end-sig_surv_len+1:end,:); % 保留后sig_surv_len行

%% 计算时域滤波器权向量
alpha = inv(X'*X)*(X')*sig_surv;

%% 去除监测信号中的杂波分量
sig_ECA = sig_surv-X*alpha;











