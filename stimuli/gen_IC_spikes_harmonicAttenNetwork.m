% use BOSSA to create IC spikes, as the input to the AIM network

% path of BOSSA
bossapath = fullfile('..','BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'recon'])
addpath([bossapath filesep 'HRTF'])
addpath([bossapath filesep 'ObjectiveMeasure'])
addpath([bossapath filesep 'Plotting'])

% load stimuli files
load(fullfile('Stimuli','CRM Stimuli TM 0deg.mat'),'wavs');

ICrms = 150; %input gain
data = struct();
speakerSet = 1:19;
for i = speakerSet
    mixed = wavs(i).mixed;
    mixed = mixed./mean(rms(mixed)).*ICrms; %adjust average rms of mixture

    % The next section of code is copied from BOSSA's demo file :)

    % ====================== Peripheral model ======================
    % The first part of BOSS is the peripheral filtering:
    % mixedL and mixedR are the two channels of the sound mixture, with
    % sampling frequency fs

    %peripheral filter parameters
    low_freq = 50; %min freq of the filter
    high_freq = 8000;
    numChannel = 256;
    fs = 40000;
    [nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
    fcoefs=MakeERBFilters(fs,cf,low_freq);

    % inputs before filtering should have rms of ~150;
    s_filt.sL=ERBFilterBank(mixed(:,1),fcoefs); %freq*time
    s_filt.sR=s_filt.sL;
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
    randomness = 0;
    azList = [0]; % Neuron preferred directions
    tic
    [spk_IC, firingrate] = ICmodel(s_filt,azList,randomness);
    toc

    % save spikes
    save(sprintf('CRM MF 0deg colocated%sIC_spks_0DegMF_Set_%02i.mat',filesep,i),'spk_IC','fcoefs','fs');
end

%% debug/evaluate
% plotSpikeRasterFs(spk_IC,'PlotType','vertline','Fs',40000)
% set(gca,'Ytick',1:nf/16:nf,'YtickLabel',cf(1:nf/16:nf))
% spks = spkTime2Train(spk_IC,fs,0.3*fs);
% winLen = fs*0.01;
% fr = conv2(spks,ones(winLen,1));
% figure;imagesc(fr')
% 
%
% addpath('C:\Users\Kenny\Desktop\GitHub\Plotting')
% load('stimuli\ICrms100\IC_spks_trioSet_01.mat','spk_IC')
% spks = spkTime2Train(spk_IC,fs,length(mixed));
% figure;
% for i = 1:5
%     subplot(1,5,i)
%     icSpikes = logical(squeeze(spks(:,:,i))');
%     plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',40000);
%     xlim([0 2000])
%     if i==1, ylabel('IC spikes'); end
%     set(gca,'Ytick',[1:64],'YtickLabel',cf)
% end