% use BOSSA to create IC spikes, as the input to the AIM network

% path of BOSSA
bossapath = ('C:\Users\Kenny\Desktop\GitHub\BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'recon'])
addpath([bossapath filesep 'HRTF'])
addpath([bossapath filesep 'ObjectiveMeasure'])

% define stimuli
rmsGain = 150;
stimuliLoc = [-80:10:80]; %match loudspeaker locations
stimuliDur = 0.15; %seconds
fs = 40000; %default hrtf fs
WGN = wgn(stimuliDur*fs,1,1,'real');
WGN = WGN/rms(WGN)*rmsGain;

for i = 1:length(stimuliLoc)
    hrtf_name = ['kemar_small_horiz_' num2str(stimuliLoc(i)) '_0.mat'];
    hrir = load([bossapath filesep 'HRTF' filesep hrtf_name]);
    hrtfL(:,i) = hrir.hrir_left;
    hrtfR(:,i) = hrir.hrir_right;
end

% apply impulse responses
WGNL = fftfilt(hrtfL,WGN');
WGNR = fftfilt(hrtfR,WGN');

low_freq = 200; %min freq of the filter
high_freq = 8000;
numChannel = 64;
 
%apply peripheral filters
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

data = struct();
for i = 1:length(stimuliLoc)
    tic
    s_filt = struct();
    s_filt.sL=ERBFilterBank(WGNL(:,i),fcoefs); %freq*time
    s_filt.sR=ERBFilterBank(WGNR(:,i),fcoefs);
%     s_filt.sL = WGNL(:,i)';
%     s_filt.sR = WGNR(:,i)';
    s_filt.band = 'narrow';
    s_filt.BW = bw;
    s_filt.flist = cf;
    s_filt.nf = nf;
    s_filt.fs = fs;
    s_filt.lowFreq = low_freq;
    s_filt.highFreq = high_freq;

    % IC model module.
    randomness = 1;
    tic
    azList = [-90:5:90]; %model neuron locations
    [spks(i).spk_IC, spks(i).FR] = ICmodel(s_filt,azList,randomness);
    toc
    
    saveLoc = 'C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\stimuli\WGNs\';
    if ~exist(saveLoc,'dir')
        mkdir(saveLoc);
    end
    spk_IC = spks(i).spk_IC;
    FR = spks(i).FR;
    saveName = sprintf('Spatial Tuning SingleWGN Position %02i.mat',stimuliLoc(i));
    save([saveLoc filesep saveName],'spk_IC','FR','stimuliLoc','azList','fs','cf','fcoefs');
    toc
end
disp('done')