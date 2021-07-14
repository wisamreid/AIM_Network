% This script creates the stimuli & pre-cortical responses for the
% middlebrooks simulation. It uses BOSSA to create IC spikes from the 
% spatialized WGNs.

% path of BOSSA
bossapath = ('C:\Users\Kenny\Desktop\GitHub\BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'recon'])
addpath([bossapath filesep 'HRTF'])
addpath([bossapath filesep 'ObjectiveMeasure'])

% data storage location
saveLoc = 'C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\stimuli\WGNs\';
if ~exist(saveLoc,'dir')
    mkdir(saveLoc);
end
    
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
    
    spk_IC = spks(i).spk_IC;
    FR = spks(i).FR;
    saveName = sprintf('Spatial Tuning SingleWGN Position %02i.mat',stimuliLoc(i));
    save([saveLoc filesep saveName],'spk_IC','FR','stimuliLoc','azList','fs','cf','fcoefs');
    toc
end
disp('done')