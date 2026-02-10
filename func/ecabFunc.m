function sig_ECAB = ecabFunc(sig_surv,sig_ref,delayDot)
% ECAB
N = length(sig_surv);
num = 245760;%每段采样点
L = ceil(N/num);%子带分段数

%信号分段
sig_surv = [sig_surv;zeros(L*num-length(sig_surv),1)];
sig_ref = [sig_ref;zeros(L*num-length(sig_ref),1)];
sig_surv = reshape(sig_surv,[num,L]);
sig_ref = reshape(sig_ref,[num,L]);

sig_ECAB = zeros(num,L);
for i = 1:L
    if any(sig_ref(:,i))
        sig_ECAB(:,i) = ecaFunc(sig_surv(:,i),sig_ref(:,i),delayDot);
    else
        sig_ECAB(:,i) = sig_surv(:,i);
    end
end
sig_ECAB = reshape(sig_ECAB,[num*L,1]);
sig_ECAB = sig_ECAB(1:N);
        
end