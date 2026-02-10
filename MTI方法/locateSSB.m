clear,clc
tStart=tic;
%% 读取数据
path='X:\数据\外辐射源雷达\数据\20231102\';

% 监测信道
file_surv='horn_move_30m';
filename_surv=strcat(path,file_surv);
sig_surv=read_complex_binary(filename_surv);
sig_surv=upsample(sig_surv,2);

% 参考信道
file_ref='linear_move_30m';
filename_ref=strcat(path,file_ref);
sig_ref=read_complex_binary(filename_ref);
sig_ref=upsample(sig_ref,2);

fprintf('====================\n');tic;
fprintf('SSB重构\n')
[ssb,pss]= ssbReconstructionFunc(sig_ref,122.88e6);
tic

down_factor = 2;
sig_surv = downsample(sig_surv,down_factor);
sig_ref = downsample(sig_ref,down_factor);
ssb = downsample(ssb,down_factor);
%% 相关参数
c=299792458;%光速
fs=122.88e6/down_factor; %采样率
f0=2.52495e9; %载频
ssb_bw=7.2e6;%带宽
T_pulse = 20e-3; %脉冲周期
Num_fast = fs*T_pulse; %快时间维点数
Num_slow = floor(length(sig_ref) / Num_fast); %慢时间维点数
Num_slow = ceil(Num_slow/2)*2-1; %取奇数
t_fast=(0:1/fs:(Num_fast-1)/fs).'; %快时间向量
t_slow=(0:Num_fast/fs:(Num_slow-1)*Num_fast/fs); %慢时间向量
dpl=(-Num_slow/2+0.5:1:Num_slow/2-0.5)*fs/Num_fast/Num_slow;%多普勒向量

%% 信号按周期排列
sig_surv = sig_surv(1:(Num_fast*Num_slow));
sig_surv = reshape(sig_surv,Num_fast,Num_slow); %监测信号矩阵

sig_ref = sig_ref(1:(Num_fast*Num_slow));
sig_ref = reshape(sig_ref,Num_fast,Num_slow); %参考信号矩阵

ssb = [ssb;zeros(Num_fast-length(ssb),1)];

%% 信号匹配滤波
% 参考信号匹配滤波
sig_ref_pc = zeros(Num_fast,Num_slow); %匹配滤波后的参考信号矩阵
maxindex_sig_ref_pc = zeros(1,Num_slow); %每个周期内最大值位置 

fprintf('====================\n');tic;
fprintf('参考信号匹配滤波\n')
backNum = 0;
for i = 1:Num_slow
    sig_ref_pc(:,i) = pulseCompressionFunc(sig_ref(:,i),ssb);
    [~,maxindex_sig_ref_pc(i)] = max(sig_ref_pc(:,i));
    
    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 监测信号匹配滤波
sig_surv_pc = zeros(Num_fast,Num_slow); %匹配滤波后的监测信号矩阵
maxindex_sig_surv_pc = zeros(1,Num_slow); %每个周期内最大值位置 

fprintf('====================\n');tic;
fprintf('监测信号匹配滤波\n')
backNum = 0;
for i = 1:Num_slow
    sig_surv_pc(:,i) = pulseCompressionFunc(sig_surv(:,i),ssb);
    [~,maxindex_sig_surv_pc(i)] = max(sig_surv_pc(:,i));

    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 作图
figure;

subplot(1,2,1)
x = 1:Num_slow;
y = min(maxindex_sig_ref_pc)-500:max(maxindex_sig_ref_pc)+500;
z = abs(sig_ref_pc(y,:));
imagesc(x,y,z)
xlabel('slow time')
ylabel('fast time')
title('参考信号匹配滤波')

subplot(1,2,2)
x = 1:Num_slow;
y = min(maxindex_sig_surv_pc)-500:max(maxindex_sig_surv_pc)+500;
z = abs(sig_surv_pc(y,:));
imagesc(x,y,z)
xlabel('slow time')
ylabel('fast time')
title('监测信号匹配滤波')

%% 时间同步
zeroindex_ref = maxindex_sig_ref_pc(1);
zeroindex_surv = maxindex_sig_surv_pc(1);

% 参考信号时间同步
sig_ref_aline = zeros(Num_fast,Num_slow); %同步后的参考信号矩阵
sig_ref_pc_aline = zeros(Num_fast,Num_slow); %匹配滤波后同步的参考信号矩阵

fprintf('====================\n');tic;
fprintf('参考信号时间同步\n')
backNum = 0;
for i = 1:Num_slow
    sig_ref_aline(:,i) = circshift(sig_ref(:,i),zeroindex_ref-maxindex_sig_ref_pc(i));
    sig_ref_pc_aline(:,i) = circshift(sig_ref_pc(:,i),zeroindex_ref-maxindex_sig_ref_pc(i));

    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 监测信号时间同步(与参考信号做同样变换)
sig_surv_aline = zeros(Num_fast,Num_slow); %同步后的参考信号矩阵
sig_surv_pc_aline = zeros(Num_fast,Num_slow); %匹配滤波后同步的参考信号矩阵

fprintf('====================\n');tic;
fprintf('监测信号时间同步\n')
backNum = 0;
for i = 1:Num_slow
    sig_surv_aline(:,i) = circshift(sig_surv(:,i),zeroindex_surv-maxindex_sig_surv_pc(i));
    sig_surv_pc_aline(:,i) = circshift(sig_surv_pc(:,i),zeroindex_surv-maxindex_sig_surv_pc(i));
    
    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 作图
figure;

subplot(1,2,1)
x = 1:Num_slow;
y = min(maxindex_sig_ref_pc)-500:max(maxindex_sig_ref_pc)+500;
z = abs(sig_ref_pc_aline(y,:));
imagesc(x,y,z)
xlabel('slow time')
ylabel('fast time')
title('参考信号匹配滤波(同步)')

subplot(1,2,2)
x = 1:Num_slow;
y = min(maxindex_sig_surv_pc)-500:max(maxindex_sig_surv_pc)+500;
z = abs(sig_surv_pc_aline(y,:));
imagesc(x,y,z)
xlabel('slow time')
ylabel('fast time')
title('监测信号匹配滤波(同步)')

%% 相干积累
% 参考信号相干积累
fprintf('====================\n');tic;
fprintf('参考信号相干积累\n')
sig_ref_pc_aline_fft = fftshift(fft(sig_ref_pc_aline,Num_slow,2),2);
toc

% 监测信号相干积累
fprintf('====================\n');tic;
fprintf('监测信号相干积累\n')
sig_surv_pc_aline_fft = fftshift(fft(sig_surv_pc_aline,Num_slow,2),2);
toc

% 作图
figure;

subplot(1,2,1)
x = dpl;
y = min(maxindex_sig_ref_pc)-500:max(maxindex_sig_ref_pc)+500;
z = abs(sig_ref_pc_aline_fft(y,:));
imagesc(x,y,z)
xlabel('dpl')
ylabel('fast time')
title('参考信号相干积累')

subplot(1,2,2)
x = dpl;
y = min(maxindex_sig_surv_pc)-500:max(maxindex_sig_surv_pc)+500;
z = abs(sig_surv_pc_aline_fft(y,:));
imagesc(x,y,z)
xlabel('dpl')
ylabel('fast time')
title('监测信号相干积累')

%% 相位补偿
[~,phaseindex]=find(sig_ref_pc_aline_fft == max(max(sig_ref_pc_aline_fft)));
phase = -dpl(phaseindex);

% 参考信号相位补偿
sig_ref_aline_phase = zeros(Num_fast,Num_slow); %同步后进行相位补偿的参考信号矩阵

fprintf('====================\n');tic;
fprintf('参考信号相位补偿\n')
backNum = 0;
for i = 1:Num_slow
    sig_ref_aline_phase(:,i) = sig_ref_aline(:,i).*exp(1i*2*pi*phase.*t_slow(i));
    
    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 监测信号相位补偿
sig_surv_aline_phase = zeros(Num_fast,Num_slow); %同步后进行相位补偿的监测信号矩阵

fprintf('====================\n');tic;
fprintf('监测信号相位补偿\n')
backNum = 0;
for i = 1:Num_slow
    sig_surv_aline_phase(:,i) = sig_surv_aline(:,i).*exp(1i*2*pi*phase.*t_slow(i));
    
    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 参考信号匹配滤波
sig_ref_aline_phase_pc = zeros(Num_fast,Num_slow); %匹配滤波后的参考信号矩阵

fprintf('====================\n');tic;
fprintf('补偿后参考信号匹配滤波\n')
backNum = 0;
for i = 1:Num_slow
    sig_ref_aline_phase_pc(:,i) = pulseCompressionFunc(sig_ref_aline_phase(:,i),ssb);

    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 监测信号匹配滤波
sig_surv_aline_phase_pc = zeros(Num_fast,Num_slow); %匹配滤波后的监测信号矩阵

fprintf('====================\n');tic;
fprintf('补偿后监测信号匹配滤波\n')
backNum = 0;
for i = 1:Num_slow
    sig_surv_aline_phase_pc(:,i) = pulseCompressionFunc(sig_surv_aline_phase(:,i),ssb);

    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 参考信号相干积累
fprintf('====================\n');tic;
fprintf('补偿后参考信号相干积累\n')
sig_ref_aline_phase_pc_fft = fftshift(fft(sig_ref_aline_phase_pc,Num_slow,2),2);
toc

% 监测信号相干积累
fprintf('====================\n');tic;
fprintf('补偿后监测信号相干积累\n')
sig_surv_aline_phase_pc_fft = fftshift(fft(sig_surv_aline_phase_pc,Num_slow,2),2);
toc

% 作图
figure;

subplot(1,2,1)
x = dpl;
y = min(maxindex_sig_ref_pc)-500:max(maxindex_sig_ref_pc)+500;
z = abs(sig_ref_aline_phase_pc_fft(y,:));
imagesc(x,y,z)
xlabel('dpl')
ylabel('fast time')
title('补偿后参考信号相干积累')

subplot(1,2,2)
x = dpl;
y = min(maxindex_sig_surv_pc)-500:max(maxindex_sig_surv_pc)+500;
z = abs(sig_surv_aline_phase_pc_fft(y,:));
imagesc(x,y,z)
xlabel('dpl')
ylabel('fast time')
title('补偿后监测信号相干积累')

%% MTI
fprintf('====================\n');tic;
fprintf('监测信号MTI\n')
tmp = circshift(sig_surv_aline_phase,1,2);
sig_surv_aline_phase_mti = sig_surv_aline_phase - tmp;
toc

% 监测信号匹配滤波
sig_surv_aline_phase_mti_pc = zeros(Num_fast,Num_slow); %匹配滤波后的监测信号矩阵

fprintf('====================\n');tic;
fprintf('去杂波后监测信号匹配滤波\n')
backNum = 0;
for i = 1:Num_slow
    sig_surv_aline_phase_mti_pc(:,i) = pulseCompressionFunc(sig_surv_aline_phase_mti(:,i),ssb);

    % 显示进度
    fprintf(repmat('\b',1,backNum)); %输出退格          
    backNum = fprintf('%d/%d',i,Num_slow); %输出进度条
end
fprintf('\n');
toc

% 监测信号相干积累
fprintf('====================\n');tic;
fprintf('去杂波后监测信号相干积累\n')
sig_surv_aline_phase_mti_pc_fft = fftshift(fft(sig_surv_aline_phase_mti_pc,Num_slow,2),2);
toc

% 作图
figure;

x = dpl;
y = min(maxindex_sig_surv_pc)-500:max(maxindex_sig_surv_pc)+500;
z = abs(sig_surv_aline_phase_mti_pc_fft(y,:));
imagesc(x,y,z)
xlabel('dpl')
ylabel('fast time')
title('去杂波后监测信号相干积累')

%% 画结果图
% 初始数据
dot_range = ceil(50/down_factor);
% [~,zeroindex_surv] = max(sig_surv_aline_phase_pc_fft(:,Num_slow/2+0.5));
x = dpl;
y = linspace(-dot_range,dot_range,2*dot_range+1);
z = abs(sig_surv_aline_phase_mti_pc_fft(y+zeroindex_surv,:));


% 单位转换
v = x*c/(2*f0);
r = y*c/fs;

% 均值滤波
ave_filter=fspecial('average',[1*dot_range 1*dot_range]);
mag=imfilter(imresize(z,[20*dot_range+1,20*dot_range+1]),ave_filter);
v = imresize(v,[1,20*dot_range+1]);
r = imresize(r,[1,20*dot_range+1]);
mag = mag/max(max(mag));

% 最大值
[r_idx,v_idx] = find(mag == max(mag(:)));

% 作图
figure;
imagesc(v,r,mag);
colorbar;xlabel('速度(m/s)');ylabel('距离(m)');title('RD图');
text(v(v_idx),r(r_idx),['\leftarrow','(',num2str(r(r_idx)),',',num2str(v(v_idx)),',',num2str(mag(r_idx,v_idx)),')']);



%%
fprintf('====================\n运行完成,');
toc(tStart)