%hSSBurstFrequencyCorrect Frequency offset correction of SS burst waveform

%   Copyright 2020 The MathWorks, Inc.


function [rxWaveform, freqOffset, NID2] = hSSBurstFrequencyCorrect_func(rxWaveform,ssbBlockPattern,rxSampleRate,searchBW)
N=2082:2239;
M=[1,3,5];
i=0;j=0;
for n=1:length(N)
    for m=1:length(M)
        Temp=3*N(n)+(M(m)-3)/2;
        if(Temp>=6246&&Temp<=7017)
            i=i+1;
            frq(i)=(N(n)*1200e3+M(m)*50e3)/1e9;%n41频段的同步栅格
            if(frq(i)>=2.515&&frq(i)<=2.8)
                j=j+1;
                frq_save(j)=frq(i);%在2.515GHz--2.615GHz的同步栅格
            end
        end
    end
end
frq_final=(frq_save-2.8)*1e3;%下变频以后可能的同步栅格位置/MHz


scs = hSSBurstSubcarrierSpacing(ssbBlockPattern);
%     scs=30e3;
if nargin < 3
    searchBW = 6*scs;
end

syncNfft = 256; % minimum FFT size to cover SS burst
syncSR = syncNfft*scs*1e3;
nrbSSB = 20;
syncOfdmInfo = nrOFDMInfo(nrbSSB, scs,'SampleRate',syncSR,'Nfft',syncNfft);

% PSS subcarrier locations in a 240 subcarrier (SSB) grid
subsPSS = nrPSSIndices('IndexStyle','subscript');
kPSS = subsPSS(:,1);

% Create reference grid for timing estimation containing one PSS. The
% PSS is placed in the second OFDM symbol of the reference grid to
% avoid the special CP length of the first OFDM symbol
refGrid = zeros([nrbSSB*12 2]);

% Frequency offset and PSS search
fshifts = (-searchBW:0.015e6:searchBW) ; % Half subcarrier step
peak_value = zeros(numel(fshifts),3);
peak_index = zeros(numel(fshifts),3);
t = (0:size(rxWaveform,1)-1).' / rxSampleRate;


for syn=26
    frq_final(syn)=0;
    for fIdx = 1:numel(fshifts)

        coarseFrequencyOffset = fshifts(fIdx);
        %         coarseFrequencyOffset = 0e3;

        rxWaveformFreqCorrected = rxWaveform .*exp(-1i*2*pi*frq_final(syn)*1e6*t).* exp(-1i*2*pi*coarseFrequencyOffset*t);

        % Downsample to the minumum sampling rate to cover SSB bandwidth
        rxWaveformDS = resample(rxWaveformFreqCorrected,syncSR,rxSampleRate);

        for NID2 = [0 1 2]
            %             NID2=0;
            refGrid(kPSS,2,NID2+1) = nrPSS(NID2);

            nSlot = 0;
            [~,corr] = nrTimingEstimate(rxWaveformDS,nrbSSB,scs,nSlot,refGrid(:,:,NID2+1),'SampleRate',syncSR,'Nfft',syncNfft);
            corr = sum(abs(corr),2);

            [peak_value(fIdx,NID2+1),peak_index(fIdx,NID2+1)] = max(corr);
            peak_index(fIdx,NID2+1) = peak_index(fIdx,NID2+1) + syncOfdmInfo.SymbolLengths(1);

        end
    end

    % Plot PSS correlations
    %     figure(syn) ;
    %     hold on;
    %     plot(fshifts/1e3,peak_value);
    %     title('PSS Correlations versus Frequency Offset');
    %     ylabel('Magnitude');
    %     xlabel('Frequency Offset (kHz)');

    % Determine NID2 and coarse frequency offset by finding the strongest
    % correlation
    [fIdx,NID2] = find(peak_value==max(peak_value(:)));
    coarseFrequencyOffset = fshifts(fIdx);
    NID2 = NID2 - 1;

    % Apply coarse frequency correction
    rxWaveformFreqCorrected = rxWaveform .* exp(-1i*2*pi*coarseFrequencyOffset*t);
    %      rxWaveformFreqCorrected = rxWaveform .*exp(-1i*2*pi*(-40.15)*1e6*t).* exp(-1i*2*pi*coarseFrequencyOffset*t);
    % Downsample received waveform after coarse frequency correction
    rxWaveformDS = resample(rxWaveformFreqCorrected,syncSR,rxSampleRate);

    % Plot selected NID2
    % % %     plot(coarseFrequencyOffset/1e3,peak_value(fIdx,NID2+1),'kx','LineWidth',2,'MarkerSize',8);
    % % %     lgd = legend;
    % % %     lgd.Interpreter = 'latex';
    % % %     legends = "$N_{ID}^{(2)}$ = " + num2cell(0:2);
    % % %     legend([legends "coarse $\Delta_f$ = " + num2str(coarseFrequencyOffset) + ", $N_{ID}^{(2)}$ = " + num2str(NID2)],'Location','East');
    % % %     end
end
% Determine timing offset
offset = peak_index(fIdx,NID2+1) - 1;

% Perform fine frequency offset estimation using CP correlation across
% the 4 OFDM symbols of the SSB
fineFrequencyOffset = hSSBurstFineFrequencyOffset(rxWaveformDS(1+offset:end,:),syncOfdmInfo);

freqOffset = coarseFrequencyOffset + fineFrequencyOffset;

% Apply coarse and fine frequency correction to the received waveform
rxWaveform = rxWaveform .* exp(-1i*2*pi*freqOffset*t);
%    rxWaveform = rxWaveform .* exp(-1i*2*pi*frq_final(syn)*1e6*t).* exp(-1i*2*pi*freqOffset*t);

end

function frequencyOffset = hSSBurstFineFrequencyOffset(waveform,ofdmInfo)

% Get 'Lsym', the number of time domain samples in an OFDM symbol,
% which is the sum of 'Lcp' and 'Lu', the number of cyclic prefix
% samples and useful samples respectively
Lcp = ofdmInfo.CyclicPrefixLengths(2);
Lu = ofdmInfo.Nfft;
Lsym = Lcp + Lu;

% Multiply the waveform by itself delayed by Lu samples and conjugated
delayed = [zeros(Lu,size(waveform,2)); waveform(1:end-Lu,:)];
cpProduct = waveform .* conj(delayed);

% Apply a moving sum filter with a window size equal to the CP length
cpXCorr = filter(ones([Lcp 1]),1,cpProduct);

% Moving sum over 4 OFDM symbols (i.e. the size of the SS block)
y = cpXCorr;
cpXCorrDelayed = cpXCorr;
for k = 1:3
    cpXCorrDelayed = [zeros(Lsym,size(waveform,2)); cpXCorrDelayed(1:end-Lsym,:)];
    y = y + cpXCorrDelayed;
end

% Extract the correlation peak, average over the receive antennas,
% then compute the phase and corresponding frequency offset
cpCorrIndex = Lu + Lcp + 3*Lsym; % 4*Lsym = Lu (delayed correlation) + Lcp (sum filter CP) + 3*Lsym (moving sum filter of 4 symbols)
scs = ofdmInfo.SampleRate/ofdmInfo.Nfft;
frequencyOffset =  scs * angle(mean(y(cpCorrIndex,:))) / (2*pi);

end

