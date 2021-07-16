% use BOSSA to create IC spikes, as the input to the AIM network
% stimuli into the IC model are pure tones
% BOSSA needs to have "best ILDs & best ITDs" = 0 for all frequencies.

% path of BOSSA
% bossapath = ('C:\Users\Kenny\Desktop\GitHub\BOSSA');
bossapath = fullfile('..','BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'HRTF'])
addpath([bossapath filesep 'recon'])
addpath([bossapath filesep 'plotting'])
%addpath([bossapath filesep 'ObjectiveMeasure'])

% cosine ramp
N = 20;
r = linspace(-pi,0,N);
w = (cos(r)+1)/2;

% stimuli
lowFreq = 200;
highFreq = 16000;
numTrials = 256;
rmsGain = 20; %input gain
fs = 44100;
f0 = logspace(log10(lowFreq),log10(highFreq),numTrials); %hz
tmax = 150/1000; %milliseconds; time axis
t = 0:1/fs:tmax; %
tones = zeros(length(f0),length(t));
for i = 1:length(f0)
    tones(i,:) = cos(2*pi*f0(i)*t);
    tones(i,1:N) = tones(i,1:N).*w;
    tones(i,end-N+1:end) = tones(i,end-N+1:end).*fliplr(w);
    tones(i,:) = tones(i,:)./rms(tones(i,:))*rmsGain;
end
tones = [tones zeros(size(tones))];
% imagesc((1:length(t))/fs,f0,tones)

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


data = struct();
h = figure('position',[300 400 1500 300],'visible','off');
for i = 1:size(tones,1)
    % The next section of code is copied from BOSSA's demo file :)

    % ====================== Peripheral model ======================
    % The first part of BOSS is the peripheral filtering:
    % mixedL and mixedR are the two channels of the sound mixture, with
    % sampling frequency fs

    %peripheral filter parameters
    low_freq = lowFreq; %min freq of the filter
    high_freq = highFreq;
    numChannel = numTrials;
    [nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
    fcoefs=MakeERBFilters(fs,cf,low_freq);
    
    ILDs = gen_ILD(low_freq,high_freq,numChannel,fs,[-90 0 90]);
    sig_filt = ERBFilterBank(tones(i,:),fcoefs);
    s_filt.sR = sig_filt;
    s_filt.sL = s_filt.sR;
    
    % inputs before filtering should have rms of ~150;
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
    inc = numChannel/8;
    ylabel('frequency (kHz)')
    yticks(1:inc:numChannel)
    yticklabels(round(cf(1:inc:numChannel)/1000,2))
    suptitle(['stimulus frequency: ' num2str(f0(i),'%0.2f') 'hz'])
    saveas(h,sprintf('IC spikes tones_%02i.tiff',i))
    clf
end

%% plot IC response to all stimuli
% stimuliRoot = 'C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\stimuli\PureTones';
stimuliRoot = fullfile('AIM_network','stimuli','PureTones');
total = [];
for trial = 1:numTrials
    spks = load([stimuliRoot filesep sprintf('IC_spks_tones_%02i.mat',trial)]);
    spk_IC = spkTime2Train(spks.spk_IC,spks.fs,fs*0.3);
    spk_IC = squeeze(spk_IC(:,:,3));
    total(:,trial) = sum(spk_IC);
end
figure;imagesc(total)
numTicks = 10;
inc = floor(numTrials/numTicks);
xlabel('stimulus frequency (kHz)')
xticks(1:inc:numTrials)
xticklabels(round(f0(1:inc:numTrials)/1000,2))
ylabel('frequency channel')
yticks(mod(numTrials,inc):inc:numTrials)
yticklabels(round(cf(mod(numTrials,inc):inc:numTrials)/1000,2))
