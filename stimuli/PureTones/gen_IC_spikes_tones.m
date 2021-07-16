% use BOSSA to create IC spikes, as the input to the AIM network

% path of BOSSA
% bossapath = ('C:\Users\Kenny\Desktop\GitHub\BOSSA');
bossapath = fullfile('..','BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'HRTF'])
%addpath([bossapath filesep 'recon'])
%addpath([bossapath filesep 'ObjectiveMeasure'])

% stimuli
fs = 44100;
numTrials = 64;
f0 = logspace(log10(500),log10(16000),numTrials); %hz
tmax = 150/1000; %milliseconds; time axis
t = 0:1/fs:tmax; %
t = [t zeros(1,length(t))]; % pad equial length of time with zeros
tones = zeros(length(f0),length(t));
for i = 1:length(f0)
    tones(i,:) = cos(2*pi*f0(i)*t);
end
imagesc((1:length(t))/fs,f0,tones)

% HRTFs
% hrtfpath = 'C:\Users\Kenny\Desktop\GitHub\BOSSA\HRTF';
hrtfpath = fullfile(bossapath,'HTRF');
talker_azs = 0;
numtalkers = length(talker_azs);
for i = 1:numtalkers
    hrtf_name = ['kemar_small_horiz_' num2str(talker_azs(i)) '_0.mat'];
    hrir = load([hrtfpath filesep hrtf_name]);
    hrtfL(:,i) = hrir.hrir_left;
    hrtfR(:,i) = hrir.hrir_right;
end

ICrms = 150; %input gain
data = struct();
h = figure('position',[300 400 1500 300],'visible','off');
for i = 1:size(tones,1)
    mixedL = fftfilt(hrtfL,tones(i,:));
    mixedR = fftfilt(hrtfR,tones(i,:));
    mixed = [mixedL; mixedR]';
    
    mixed = mixed./mean(rms(mixed)).*ICrms; %adjust average rms of mixture

    % The next section of code is copied from BOSSA's demo file :)

    % ====================== Peripheral model ======================
    % The first part of BOSS is the peripheral filtering:
    % mixedL and mixedR are the two channels of the sound mixture, with
    % sampling frequency fs

    %peripheral filter parameters
    low_freq = 500; %min freq of the filter
    high_freq = 16000;
    numChannel = 64;
    [nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
    fcoefs=MakeERBFilters(fs,cf,low_freq);
    
    ILDs = gen_ILD(low_freq,high_freq,numChannel,fs,0);
    sig_filt = ERBFilterBank(tones(i,:));
    for j = 1:numChannel
        s_filt.sL(j,:) = sig_filt*ILD(j);
    end
    s_filt.sR = sig_filt;
    
    % inputs before filtering should have rms of ~150;
    s_filt.sL=ERBFilterBank(mixed(:,1),fcoefs); %freq*time
    s_filt.sR=ERBFilterBank(mixed(:,2),fcoefs);
    s_filt.band = 'narrow';
    s_filt.BW = bw;
    s_filt.flist = cf;
    s_filt.nf = nf;
    s_filt.fs = fs;
    s_filt.lowFreq = low_freq;
    s_filt.highFreq = high_freq;

    % ===================== Fischer's IC model ======================
    % The second part of BOSS takes the filtered sound mixture and computes
    % neural responses that correspond to each spatial location and frequency

    % inputs: filtered L and R input channels and other parameters above
    options.rand = 0;
    azList = [-90,-45,0,45,90]; % Neuron preferred directions
    tic
    [spk_IC, firingrate] = ICmodel(s_filt,azList,options);
    toc

    % save spikes
    save(sprintf('IC_spks_tones_%02i.mat',i),'spk_IC','fcoefs','fs','cf');
    
    %visualize
    spikes = spkTime2Train(spk_IC,fs);
    LineFormat.LineWidth = 1;
    for loc = 1:5
        subplot(1,5,loc);
        spikes1 = squeeze(spikes(:,:,loc))';
        plotSpikeRasterFs(logical(spikes1), 'PlotType','vertline', 'Fs',fs,'LineFormat',LineFormat);
        xlim([0 170])
    end
    subplot(1,5,1); 
    ylabel('frequency (kHz)')
    yticks(1:8:64)
    yticklabels(round(cf(1:8:64)/1000,2))
    suptitle(['stimulus frequency: ' num2str(f0(i),'%0.2f') 'hz'])
    saveas(h,sprintf('IC spikes tones_%02i.tiff',i))
    clf
end

