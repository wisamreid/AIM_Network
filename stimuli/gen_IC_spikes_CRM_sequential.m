% use BOSSA to create IC spikes, as the input to the AIM network

% path of BOSSA
bossapath = ('C:\Users\Kenny\Desktop\GitHub\BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'recon'])
addpath([bossapath filesep 'HRTF'])
addpath([bossapath filesep 'ObjectiveMeasure'])

% load stimuli files - may already be spatialized
load('Stimuli/CRM Stimuli TM 90deg.mat','wavs'); 
% talker_azs = [0 90];
% numtalkers = length(talker_azs);
% for i = 1:numtalkers
%     hrtf_name = ['kemar_small_horiz_' num2str(talker_azs(i)) '_0.mat'];
%     hrir = load([hrtfpath filesep hrtf_name]);
%     hrtfL(:,i) = hrir.hrir_left;
%     hrtfR(:,i) = hrir.hrir_right;
% end


ICrms = 150; %input gain
saveFolder = 'CRM MF 90deg separated sequential';
if ~exist(saveFolder,'dir'), mkdir(saveFolder); end
data = struct();
for i = 1:20 
    mixed = [wavs(i).tgt; wavs(i).m1];    
    mixed = mixed./mean(rms(mixed)).*ICrms; %adjust average rms of mixture

    % The next section of code is copied from BOSSA's demo file :)

    % ====================== Peripheral model ======================
    % The first part of BOSS is the peripheral filtering:
    % mixedL and mixedR are the two channels of the sound mixture, with
    % sampling frequency fs

    %peripheral filter parameters
    low_freq = 50; %min freq of the filter
    high_freq = 8000;
    numChannel = 64;
    fs = 40000;
    [nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
    fcoefs=MakeERBFilters(fs,cf,low_freq);

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
    randomness = 0;
    azList = [-90,-45,0,45,90]; % Neuron preferred directions
    tic
    [spk_IC, firingrate] = ICmodel(s_filt,azList,randomness);
    toc

    % save spikes
    save(sprintf('%s/IC_spks_set_%02i.mat',saveFolder,i),'spk_IC','fcoefs','fs');
end
