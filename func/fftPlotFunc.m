function [f,amp] = fftPlotFunc(sig,fs)

sig_len = length(sig);
s_fft = fftshift(fft(sig));
amp = abs(s_fft/sig_len);
f = (-sig_len/2+0.5:1:sig_len/2-0.5)*(fs/sig_len);
figure;plot(f,amp)

end