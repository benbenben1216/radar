function [ssb,pss]=ssbReconstructionFunc(sigt,rxSampleRate)
%% Receiver Configuration
% To synchronize and demodulate the received waveform, this information is
% needed:
%
% * The waveform sample rate to demodulate the received waveform.
% * The carrier center frequency to apply symbol phase compensation to the
% received waveform.
% * The minimum channel bandwidth to determine CORESET 0 frequency
% resources. TS 38.101-1 Table 5.3.5-1 [ <#19 1> ] describes the channel
% bandwidths for each NR band.
% * The SS block pattern (Case A...E) to determine the subcarrier spacing
% of the SS/PBCH blocks. UE searches for SS block patterns based on the NR
% operating band. For more information, see TS 38.104 Tables 5.4.3.3-1 and
% 5.4.3.3-2 [ <#19 2> ].
% * The number of SS/PBCH blocks in a burst ($L_{max}$) to calculate
% parameters of the PBCH DM-RS sequences and PBCH descrambling. These
% parameters depend on the SS/PBCH block index as described in TS 38.211
% Sections 7.3.3.1 and 7.4.1.4.1 [ <#19 3> ]. TS 38.213 Section 4.1 [ <#19
% 5> ] describes the set of SS/PBCH blocks in a burst in each case. UE
% knows the value of $L_{max}$ based on the SS block pattern and the NR
% operating band.


% Load captured waveform
% sigt=read_complex_binary("..\20230619\movecome_97_30_1_1");
pulse_time=20e-3;
pulse_len=rxSampleRate*pulse_time;
if length(sigt)>pulse_len
    rx=sigt(1:pulse_len);
else
    rx=sigt;
end
rxWaveform = rx;

% Configure receiver sample rate (samples/second)
% rxSampleRate = 122.88e6;

% Symbol phase compensation frequency. Specify the carrier center
% frequency or set to 0 to disable symbol phase compensation
fPhaseComp = 0; % Carrier center frequency (Hz)

% Set the minimum channel bandwidth for the NR band required to
% configure CORESET 0 in FR1 (See TS 38.101-1 Table 5.3.5-1)
minChannelBW = 5; % 5, 10, 40 MHz

% Configure necessary burst parameters at the receiver. The SSB pattern
% can be 'Case A','Case B','Case C' for FR1 or 'Case D','Case E' for
% FR2. The maximum number of blocks L_max can be 4 or 8 for FR1 and 64
% for FR2.
refBurst.BlockPattern = 'Case C';
refBurst.L_max = 8;

% Get OFDM information from configured burst and receiver parameters
nrbSSB = 20;
scsSSB = hSSBurstSubcarrierSpacing(refBurst.BlockPattern);
rxOfdmInfo = nrOFDMInfo(nrbSSB,scsSSB,'SampleRate',rxSampleRate);

% Display spectrogram of received waveform
% % % figure;
% % % nfft = rxOfdmInfo.Nfft;
% % % spectrogram(rxWaveform(:,1),ones(nfft,1),0,nfft,'centered',rxSampleRate,'yaxis','MinThreshold',-130);
% % % title('Spectrogram of the Received Waveform')

%% PSS Search and Frequency Offset Correction
% The receiver performs PSS search and coarse frequency offset estimation
% following these steps:
%
% * Frequency shift the received waveform with a candidate frequency
% offset. Candidate offsets are spaced half subcarrier apart. Use
% |searchBW| to control the frequency offset search bandwidth.
% * Correlate the frequency-shifted received waveform with each of the
% three possible PSS sequences (NID2) and extract the strongest correlation
% peak. The reference PSS sequences are centered in frequency. Therefore,
% the strongest correlation peak provides a measure of coarse frequency
% offset with respect to the center frequency of the carrier. The peak also
% indicates which of the three PSS (NID2) has been detected in the received
% waveform and the time instant of the best channel conditions.
% * Estimate frequency offsets below half subcarrier by correlating the
% cyclic prefix of each OFDM symbol in the SSB with the corresponding
% useful parts of the OFDM symbols. The phase of this correlation is
% proportional to the frequency offset in the waveform.

disp(' -- Frequency correction and timing estimation --')

% Specify the frequency offset search bandwidth in kHz
searchBW = 0.18e6;
rxSampleRate=122.88e6;
[rxWaveform,freqOffset,NID2] = hSSBurstFrequencyCorrect_func(rxWaveform,refBurst.BlockPattern,rxSampleRate,searchBW);
disp([' Frequency offset: ' num2str(freqOffset,'%.0f') ' Hz'])

%% Time Synchronization and OFDM Demodulation
% The receiver estimates the timing offset to the strongest SS block by
% using the reference PSS sequence dete cted in the frequency search
% process. After frequency offset correction, the receiver can assume that
% the center frequencies of the reference PSS and received waveform are
% aligned. Finally, the receiver OFDM demodulates the synchronized waveform
% and extracts the SS block.

% Create a reference grid for timing estimation using detected PSS. The PSS
% is placed in the second OFDM symbol of the reference grid to avoid the
% special CP length of the first OFDM symbol.
refGrid = zeros([nrbSSB*12 2]);
refGrid(nrPSSIndices,2) = nrPSS(NID2); % Second OFDM symbol for correct CP length

% Timing estimation. This is the timing offset to the OFDM symbol prior to
% the detected SSB due to the content of the reference grid
nSlot = 0;
timingOffset = nrTimingEstimate(rxWaveform,nrbSSB,scsSSB,nSlot,refGrid,'SampleRate',rxSampleRate);

% Synchronization, OFDM demodulation, and extraction of strongest SS block
rxGrid = nrOFDMDemodulate(rxWaveform(1+timingOffset:end,:),nrbSSB,scsSSB,nSlot,'SampleRate',rxSampleRate);
rxGrid = rxGrid(:,2:5,:);

% Display the timing offset in samples. As the symbol lengths are measured
% in FFT samples, scale the symbol lengths to account for the receiver
% sample rate.
srRatio = rxSampleRate/(scsSSB*1e3*rxOfdmInfo.Nfft);
firstSymbolLength = rxOfdmInfo.SymbolLengths(1)*srRatio;
str = sprintf(' Time offset to synchronization block: %%.0f samples (%%.%.0ff ms) \n',floor(log10(rxSampleRate))-3);
fprintf(str,timingOffset+firstSymbolLength,(timingOffset+firstSymbolLength)/rxSampleRate*1e3);

%% SSS Search
% The receiver extracts the resource elements associated to the SSS from
% the received grid and correlates them with each possible SSS sequence
% generated locally. The indices of the strongest PSS and SSS sequences
% combined give the physical layer cell identity, which is required for
% PBCH DM-RS and PBCH processing.

% Extract the received SSS symbols from the SS/PBCH block
sssIndices = nrSSSIndices;
sssRx = nrExtractResources(sssIndices,rxGrid);

% Correlate received SSS symbols with each possible SSS sequence
sssEst = zeros(1,336);
for NID1 = 0:335

    ncellid = (3*NID1) + NID2;
    sssRef = nrSSS(ncellid);
    sssEst(NID1+1) = sum(abs(mean(sssRx .* conj(sssRef),1)).^2);

end

% Plot SSS correlations
% % % figure;
% % % stem(0:335,sssEst,'o');
% % % title('SSS Correlations (Frequency Domain)');
% % % xlabel('$N_{ID}^{(1)}$','Interpreter','latex');
% % % ylabel('Magnitude');
% % % axis([-1 336 0 max(sssEst)*1.1]);

% Determine NID1 by finding the strongest correlation
NID1 = find(sssEst==max(sssEst)) - 1;

% Plot selected NID1
% % % hold on;
% % % plot(NID1,max(sssEst),'kx','LineWidth',2,'MarkerSize',8);
% % % legend(["correlations" "$N_{ID}^{(1)}$ = " + num2str(NID1)],'Interpreter','latex');

% Form overall cell identity from estimated NID1 and NID2
ncellid = (3*NID1) + NID2;

disp([' Cell identity: ' num2str(ncellid)])

%% PBCH DM-RS search
% In a process similar to SSS search, the receiver constructs each possible
% PBCH DM-RS sequence and performs channel and noise estimation. The index
% of the PBCH DM-RS with the best SNR determines the LSBs of the SS/PBCH
% block index required for PBCH scrambling initialization.

% Calculate PBCH DM-RS indices
dmrsIndices = nrPBCHDMRSIndices(ncellid);

% Perform channel estimation using DM-RS symbols for each possible DM-RS
% sequence and estimate the SNR
dmrsEst = zeros(1,8);
for ibar_SSB = 0:7

    refGrid = zeros([240 4]);
    refGrid(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
    [hest,nest] = nrChannelEstimate(rxGrid,refGrid,'AveragingWindow',[0 1]);
    dmrsEst(ibar_SSB+1) = 10*log10(mean(abs(hest(:).^2)) / nest);

end

% Plot PBCH DM-RS SNRs
% % % figure;
% % % stem(0:7,dmrsEst,'o');
% % % title('PBCH DM-RS SNR Estimates');
% % % xlabel('$\overline{i}_{SSB}$','Interpreter','latex');
% % % xticks(0:7);
% % % ylabel('Estimated SNR (dB)');
% % % axis([-1 8 min(dmrsEst)-1 max(dmrsEst)+1]);

% Record ibar_SSB for the highest SNR
ibar_SSB = find(dmrsEst==max(dmrsEst)) - 1;

% Plot selected ibar_SSB
% % % hold on;
% % % plot(ibar_SSB,max(dmrsEst),'kx','LineWidth',2,'MarkerSize',8);
% % % legend(["SNRs" "$\overline{i}_{SSB}$ = " + num2str(ibar_SSB)],'Interpreter','latex');

%% Channel Estimation using PBCH DM-RS and SSS
% The receiver estimates the channel for the entire SS/PBCH block using the
% SSS and PBCH DM-RS detected in previous steps. An estimate of the
% additive noise on the PBCH DM-RS / SSS is also performed.

refGrid = zeros([nrbSSB*12 4]);
refGrid(dmrsIndices) = nrPBCHDMRS(ncellid,ibar_SSB);
refGrid(sssIndices) = nrSSS(ncellid);
[hest,nest,hestInfo] = nrChannelEstimate(rxGrid,refGrid,'AveragingWindow',[0 1]);

%% PBCH Demodulation
% The receiver uses the cell identity to determine and extract the resource
% elements associated with the PBCH from the received grid. In addition,
% the receiver uses the channel and noise estimates to perform MMSE
% equalization. The equalized PBCH symbols are then demodulated and
% descrambled to give bit estimates for the coded BCH block.

disp(' -- PBCH demodulation and BCH decoding -- ')

% Extract the received PBCH symbols from the SS/PBCH block
[pbchIndices,pbchIndicesInfo] = nrPBCHIndices(ncellid);
pbchRx = nrExtractResources(pbchIndices,rxGrid);

% Configure 'v' for PBCH scrambling according to TS 38.211 Section 7.3.3.1
% 'v' is also the 2 LSBs of the SS/PBCH block index for L_max=4, or the 3
% LSBs for L_max=8 or 64.
if refBurst.L_max == 4
    v = mod(ibar_SSB,4);
else
    v = ibar_SSB;
end
ssbIndex = v;

% PBCH equalization and CSI calculation
pbchHest = nrExtractResources(pbchIndices,hest);
[pbchEq,csi] = nrEqualizeMMSE(pbchRx,pbchHest,nest);
Qm = pbchIndicesInfo.G / pbchIndicesInfo.Gd;
csi = repmat(csi.',Qm,1);
csi = reshape(csi,[],1);

% Plot received PBCH constellation after equalization
% % % figure;
% % % plot(pbchEq,'o');
% % % xlabel('In-Phase'); ylabel('Quadrature')
% % % title('Equalized PBCH Constellation');
% % % m = max(abs([real(pbchEq(:)); imag(pbchEq(:))])) * 1.1;
% % % axis([-m m -m m]);

% PBCH demodulation
pbchBits = nrPBCHDecode(pbchEq,ncellid,v,nest);

% Calculate RMS PBCH EVM
pbchRef = nrPBCH(pbchBits<0,ncellid,v);
evm = comm.EVM;
pbchEVMrms = evm(pbchRef,pbchEq);

% Display calculated EVM
disp([' PBCH RMS EVM: ' num2str(pbchEVMrms,'%0.3f') '%']);

%% BCH Decoding
% The receiver weights BCH bit estimates with channel state information
% (CSI) from the MMSE equalizer and decodes the BCH. BCH decoding consists
% of rate recovery, polar decoding, CRC decoding, descrambling, and
% separating the 24 BCH transport block bits from the 8 additional
% timing-related payload bits.

% Apply CSI
pbchBits = pbchBits .* csi;

% Perform BCH decoding including rate recovery, polar decoding, and CRC
% decoding. PBCH descrambling and separation of the BCH transport block
% bits 'trblk' from 8 additional payload bits A...A+7 is also performed:
%   A ... A+3: 4 LSBs of system frame number
%         A+4: half frame number
% A+5 ... A+7: for L_max=64, 3 MSBs of   the SS/PBCH block index
%              for L_max=4 or 8, A+5 is the MSB of subcarrier offset k_SSB
polarListLength = 8;
[~,crcBCH,trblk,sfn4lsb,nHalfFrame,msbidxoffset] = ...
    nrBCHDecode(pbchBits,polarListLength,refBurst.L_max,ncellid);

% Display the BCH CRC
disp([' BCH CRC: ' num2str(crcBCH)]);

% Stop processing MIB and SIB1 if BCH was received with errors
if crcBCH
    disp(' BCH CRC is not zero.');
    return
end

% Use 'msbidxoffset' value to set bits of 'k_SSB' or 'ssbIndex', depending
% on the number of SS/PBCH blocks in the burst
if (refBurst.L_max==64)
    ssbIndex = ssbIndex + (bit2int(msbidxoffset,3) * 8);
    k_SSB = 0;
else
    k_SSB = msbidxoffset * 16;
end

if 0
%% Generate waveform containing SS burst and SIB1
    % Configure the cell identity
    config = struct();
    config.NCellID = ncellid;

    % Configure an SS burst
    config.BlockPattern = 'Case C';         % FR1: 'Case A','Case B','Case C'. FR2: 'Case D','Case E'
    config.TransmittedBlocks = ones(1,8);   % Bitmap of SS blocks transmitted
    config.SubcarrierSpacingCommon = 30;    % SIB1 subcarrier spacing in kHz (15 or 30 for FR1. 60 or 120 for FR2)
    config.EnableSIB1 = 1;                  % Set to 0 to disable SIB1
    % Set the minimum channel bandwidth for the NR band required to
    % configure CORESET0 in FR1 (See TS 38.101-1 Table 5.3.5-1)
    config.MinChannelBW = 5; % 5, 10, 40 MHz

    % Configure and generate a waveform containing an SS burst and SIB1
    wavegenConfig = hSIB1WaveformConfiguration(config);
    wavegenConfig.SampleRate=122.88e6;
    % wavegenConfig.NumSubframes=0;    %Above30_a2_1 348  %Shift_30m27m_a2_1 16

    wavegenConfig.SSBurst.DMRSTypeAPosition=2;
    wavegenConfig.SSBurst.CellBarred=1;
    wavegenConfig.SSBurst.Power=0;
    wavegenConfig.SSBurst.IntraFreqReselection=0;
    bits24 = trblk;
    wavegenConfig.SSBurst.DataSource=double(bits24);
    % wavegenConfig.SSBurst.NCRBSSB=6;
    % wavegenConfig.SSBurst.PDCCHConfigSIB1
    [txWaveform,waveInfo] = nrWaveformGenerator2(wavegenConfig);
    txOfdmInfo = waveInfo.ResourceGrids(1).Info;

    % Introduce a beamforming gain by boosting the SNR of one SSB and
    % associated SIB1 PDCCH and PDSCH
    ssbIdx = ssbIndex; % Index of the SSB to boost (0-based)
    boost = 6; % SNR boost in dB
    txWaveform = hSIB1Boost(txWaveform,wavegenConfig,waveInfo,ssbIdx,boost);
    % txWaveform=txWaveform(1+8800:1228800+8800);
    % Add white Gaussian noise to the waveform
    rng('default'); % Reset the random number generator
    SNRdB = 20; % SNR for AWGN

    % rxWaveform = awgn(txWaveform,SNRdB-boost,-10*log10(double(txOfdmInfo.Nfft)));
    rxWaveform = txWaveform;
    % Configure receiver
    % Sample rate
    rxSampleRate = txOfdmInfo.SampleRate;

    % Symbol phase compensation frequency (Hz). The function
    % nrWaveformGenerator does not apply symbol phase compensation to the
    % generated waveform.
    fPhaseComp = 0; % Carrier center frequency (Hz)

    % Minimum channel bandwidth (MHz)
    minChannelBW = config.MinChannelBW;

    % Configure necessary burst parameters at the receiver
    refBurst.BlockPattern = config.BlockPattern;
    refBurst.L_max = numel(config.TransmittedBlocks);


    % Get OFDM information from configured burst and receiver parameters
    nrbSSB = 20;
    scsSSB = hSSBurstSubcarrierSpacing(refBurst.BlockPattern);
    rxOfdmInfo = nrOFDMInfo(nrbSSB,scsSSB,'SampleRate',rxSampleRate);

    % Display spectrogram of received waveform
    % % % figure;
    % % % nfft = rxOfdmInfo.Nfft;
    % % % spectrogram(rxWaveform(:,1),ones(nfft,1),0,nfft,'centered',rxSampleRate,'yaxis','MinThreshold',-130);
    % % % title('Spectrogram of the Received Waveform')

    %% MIB and BCH Parsing
    % The example parses the 24 decoded BCH transport block bits into a MIB
    % message and creates the |initialSystemInfo| structure with initial system
    % information. This process includes reconstituting the 10-bit system frame
    % number (SFN) |NFrame| from the 6 MSBs in the MIB and the 4 LSBs in the
    % PBCH payload bits. It also includes incorporating the MSB of the
    % subcarrier offset |k_SSB| from the PBCH payload bits in the case of
    % L_max=4 or 8 SS/PBCH blocks per burst.

    % Parse the last 23 decoded BCH transport block bits into a MIB message.
    % The BCH transport block 'trblk' is the RRC message BCCH-BCH-Message,
    % consisting of a leading 0 bit and 23 bits corresponding to the MIB. The
    % leading bit signals the message type transmitted (MIB or empty sequence).
    mib = fromBits(MIB,trblk(2:end)); % Do not parse leading bit

    % Create set of subcarrier spacings signaled by the 7th bit of the decoded
    % MIB, the set is different for FR1 (L_max=4 or 8) and FR2 (L_max=64)
    if (refBurst.L_max==64)
        commonSCSs = [60 120];
    else
        commonSCSs = [15 30];
    end

    initialSystemInfo = struct();
    initialSystemInfo.NFrame = mib.systemFrameNumber*2^4 + bit2int(sfn4lsb,4);
    initialSystemInfo.SubcarrierSpacingCommon = commonSCSs(mib.subCarrierSpacingCommon + 1);
    initialSystemInfo.k_SSB = k_SSB + mib.ssb_SubcarrierOffset;
    initialSystemInfo.DMRSTypeAPosition = 2 + mib.dmrs_TypeA_Position;
    initialSystemInfo.PDCCHConfigSIB1 = info(mib.pdcch_ConfigSIB1);
    initialSystemInfo.CellBarred = mib.cellBarred;
    initialSystemInfo.IntraFreqReselection = mib.intraFreqReselection;

    % Display the MIB structure
    disp(' BCH/MIB Content:')
    disp(initialSystemInfo);
    %ssb=txWaveform;
    ssb=txWaveform(131713:149249);
    pss=ssb(1:floor(length(ssb)/4));
end

if 1
%     ncellid=1002;
%     ssbIndex=4;
    % PSS SSS
    ssblock = zeros([240 4]);
    pssSymbols = nrPSS(ncellid);
    pssIndices = nrPSSIndices;
    ssblock(pssIndices) = pssSymbols;
    sssSymbols = nrSSS(ncellid);
    sssIndices = nrSSSIndices;
    ssblock(sssIndices) = 2*sssSymbols;

    % PBCH
    pbchSymbols=pbchRef;
    pbchIndices = nrPBCHIndices(ncellid);
    ssblock(pbchIndices) = 3*pbchSymbols;
    % PBCH DM-RS
    ibar_SSB = ssbIndex;
    dmrsSymbols = nrPBCHDMRS(ncellid,ibar_SSB);
    dmrsIndices = nrPBCHDMRSIndices(ncellid);
    ssblock(dmrsIndices) =  4*dmrsSymbols;

    %%时域ssb
    ofdmInfo=importdata("ofdmInfo.mat");
    initialNSlot=0;
    overlapping=0;
    hasSampleRate=1;
    %grid=[zeros(240,1),ssblock(:,1)];
    grid = ssblock;
    ssb=nr5g.internal.OFDMModulate(grid,ofdmInfo,initialNSlot,overlapping,hasSampleRate);
    pss=ssb(1:floor(length(ssb)/4));

    disp([' Reconstruction NcellID: ',num2str(ncellid)]);
    disp([' Reconstruction SSB Index: ',num2str(ssbIndex)]);
end

