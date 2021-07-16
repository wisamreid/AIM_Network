% use BOSSA to create IC spikes, as the input to the AIM network

% path of BOSSA
bossapath = fullfile('..','BOSSA');
addpath([bossapath filesep 'peripheral'])
addpath([bossapath filesep 'IC'])
addpath([bossapath filesep 'recon'])
addpath([bossapath filesep 'HRTF'])
addpath([bossapath filesep 'ObjectiveMeasure'])

% load stimuli files - may already be spatialized
% load('Stimuli/CRM Stimuli TM 90deg.mat','wavs'); 
% talker_azs = [0 90];
% numtalkers = length(talker_azs);
% for i = 1:numtalkers
%     hrtf_name = ['kemar_small_horiz_' num2str(talker_azs(i)) '_0.mat'];
%     hrir = load([hrtfpath filesep hrtf_name]);
%     hrtfL(:,i) = hrir.hrir_left;
%     hrtfR(:,i) = hrir.hrir_right;
% end


ICrms = 150; %input gain
saveFolder = 'CRM MF Xdeg separated';
if ~exist(saveFolder,'dir'), mkdir(saveFolder); end
data = struct();
for i = 1:20
    stimFile = sprintf('Stimuli/CRM Stimuli TM Xdeg/set%02i.mat',i);
    load(stimFile,'wavs','talker_azs'); 
    
    for maskerLoc = talker_azs(2:end)
    
        mixed = wavs.(['mixed' num2str(maskerLoc)]);
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
        azList = [0:15:90]; % Neuron preferred directions
        tic
        [spk_IC, firingrate] = ICmodel(s_filt,azList,randomness);
        toc

        % save spikes
        save(sprintf('%s/IC_spks_m%02i_set_%02i.mat',saveFolder,maskerLoc,i),'spk_IC','fcoefs','fs','azList');

        % evaluate
        spks = spkTime2Train(spk_IC,fs,length(mixed));

        % ===================== reconstruction  ===================== 

        % reconstruction parameters
        params = struct();
        params.fcoefs = fcoefs;
        params.cf = cf;
        params.fs = fs;
        params.delay = 0;

        params.maskRatio = 0.5;
        params.tau = 0.02;
        masks = calcSpkMask(spks,fs,params);

        % FR-mask reconstruction
        [rstim, maskedWav] = applyMask(masks.FR(:,:,azList==maskerLoc),s_filt.sL,s_filt.sR,1,'filt');
        st = runStoi(rstim,wavs.(['m' num2str(maskerLoc)]),fs);
        st2 = runStoi(rstim,wavs.tgt,fs);
        data(end+1).maskTarget = 'Masker';
        data(end).recon = rstim;
        data(end).st_onTarget = st;
        data(end).st_offTarget = st2;
        data(end).trial = i;
        data(end).masker = maskerLoc;
        stoisAttendM(i,1) = st;
        stoisAttendM(i,2) = st2;

        [rstim, maskedWav] = applyMask(masks.FR(:,:,1),s_filt.sL,s_filt.sR,1,'filt');
        st = runStoi(rstim,wavs.tgt,fs);
        st2 = runStoi(rstim,wavs.(['m' num2str(maskerLoc)]),fs);
        data(end+1).maskTarget = 'Target';
        data(end).recon = rstim;
        data(end).st_onTarget = st;
        data(end).st_offTarget = st2;
        data(end).trial = i;
        data(end).masker = maskerLoc;
        stoisAttendT(i,1) = st;
        stoisAttendT(i,2) = st2;

%     [rstim, maskedWav] = applyMask(masks.FR(:,:,5),s_filt.sL,s_filt.sR,1,'filt');
%     st = runStoi(rstim,wavs(i).m2,fs);
%     data(end+1).type = 'Right';
%     data(end).recon = rstim;
%     data(end).st = st;
%     data(end).trial = i;
%     stois(i,3) = st;    
    end
end
data(1) = []; %remove empty index
save([saveFolder filesep 'recon_data.mat'],'data','stoisAttendT','stoisAttendM')
%% debug/evaluate


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

%%
subdata = data(contains({data.maskTarget},'Target')); %attend to "target"
scatter([subdata.masker],[subdata.st_onTarget])
xlabel('masker location')
title('attend to target, STOI w/ target is high')

scatter([subdata.masker],[subdata.st_offTarget])
title('attend to target, STOI w/ masker is low')

