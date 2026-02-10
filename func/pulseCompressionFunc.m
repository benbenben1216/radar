function y = pulseCompressionFunc(sig_surv,sig_ref)
% 匹配滤波
N=length(sig_surv);
Xs = fft(sig_ref,N);        % 本地副本的FFT
Xecho = fft(sig_surv,N);  % 输入信号的FFT
Y = conj(Xs).*Xecho;   % 乘法器
y = ifft(Y,N);         % IFFT
end