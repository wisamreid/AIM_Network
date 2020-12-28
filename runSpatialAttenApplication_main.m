% runSpatialAttenApplication_main.m
% a use-case demo selective spatial attention

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('BOSSA',pathCell)), rmpath(genpath('../BOSSA')); end
addpath('mechs')
addpath('network_params')
addpath('util')
addpath('util\plotting')
addpath('util\eval_scripts')
addpath(genpath('..\dynasim'))

% =============== folders and paths ==============
% stimuliRoot = 'stimuli\CRM MF 90deg separated';
stimuliRoot = 'stimuli\CRM MF 90deg separated sequential';
expPrefix = '051';
expName = sprintf('%s_SpatialAttendDemo',expPrefix);
expRootLoc = pwd; 
dataFolder = [expRootLoc filesep 'data' filesep expName]; %things are saved to here
if ~exist(dataFolder,'dir'), mkdir(dataFolder);end
rasterFileNames = 'undecided';
seqSim = contains(stimuliRoot,'sequential'); %is sequential?

% ================ parameters ==============
attend = 0; 
attendAz = 90; %0 or 90
azList = [-90:45:90];
nLocs = 5;
nFreqs = 64;
[~,tgtTDidx] = min(abs(azList - attendAz));
trials = 8:20;

i = 1;
for trial = trials
varies = struct();
varies(1).conxn = 'IC->IC';
varies(end).param = 'trial'; 
varies(end).range = trial;

varies(end+1).conxn = 'XI';
varies(end).param = 'Itonic'; 
varies(end).range = 4;

varies(end+1).conxn = 'TD';
varies(end).param = 'Itonic'; 
varies(end).range = 3.5;

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic'; 
varies(end).range = 0;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN'; 
varies(end).range = 1.25;

varies(end+1).conxn = 'TD->XI';
varies(end).param = 'gSYN'; 
varies(end).range = 2.25;

varies(end+1).conxn = 'E->XI';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'IC->E';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

if ~seqSim
    varies(end+1).conxn = 'E->C';
    varies(end).param = 'gSYN'; 
    varies(end).range = 2.25;
else
    varies(end+1).conxn = 'E->C';
    varies(end).param = 'gSYN'; 
    varies(end).range = 2;
end


% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run' expName currentTime]);
mkdir(fullfile(study_dir, 'solve'));
    
% prep IC inputs

% if seqSim, simLen = simLen*2; end
% for trial = trials
    % ------------ load spk_IC here -----------------
    spks = load([stimuliRoot filesep sprintf('IC_spks_set_%02i.mat',trial)]);
    spk_IC = spkTime2Train(spks.spk_IC,spks.fs,spks.fs);
    spk_IC(1:400,:) = 0; %remove artifact
    simLen = size(spk_IC,1);
    save(fullfile(study_dir, 'solve', sprintf('IC_spks_t%02i.mat',trial)),'spk_IC')
% end

% other network options
options.expName = expName;
options.rasterFileNames = rasterFileNames;
options.saveLoc = dataFolder;
options.plotRasters = 0;
options.vizNetwork = 0;
options.nLocs = nLocs;
options.nFreqs = nFreqs;
options.fs = spks.fs;
options.fcoefs = spks.fcoefs;
options.saveRasters = 0;
options.simLen = simLen;

options.ECnetcon = ones(nLocs,1);
options.attend = attend;
options.tgtTDchan = tgtTDidx;
options.constantAttendLoc = 1;
TD_Imask = ones(nLocs,1);
if attend, TD_Imask(tgtTDidx) = 0; end
options.TD_Imask = TD_Imask;

% run network
tic;
[data(i),pops] = runSpatialAtten(study_dir,spk_IC,varies,options);
toc;

i = i+1;
end
%% plot C rasters
%peripheral filter parameters - match peripheral model
low_freq = 50; %min freq of the filter
high_freq = 8000;
numChannel = 64;
fs = 40000;
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

h = figure;

for trial = 1:length(trials)
    clf
    Cspikes = data(trial).C_V_spikes;
    plotSpikeRasterFs(logical(Cspikes)', 'PlotType','vertline', 'Fs',options.fs);
    xlim([0 size(spk_IC,1)/options.fs*1000])
    ylabel('frequency(kHz)')
    xlabel('time (ms)')
    majorTickLabels = [0.1,0.2,0.5,1,2,4]; %kHz
    [~,yIdx] = min(abs(round(cf/1000,2) - majorTickLabels));
    yticks(fliplr(yIdx))
    yticklabels(round(cf(fliplr(yIdx))/1000,2))

    suffix = '';
    if seqSim
        suffix = ' seq';
    end
    saveas(h,[dataFolder filesep sprintf('Cspikes trial%d attend%d azimuth%d%s.pdf',trials(trial),attend,attendAz,suffix)])
    save([dataFolder filesep sprintf('Cspikes trial%d attend%d azimuth%d%s.mat',trials(trial),attend,attendAz,suffix)],'Cspikes');
end
return;

%% reconstruct and evaluate STOIs
% load stimuli
load(['stimuli' filesep 'CRM Stimuli TM 90deg.mat'],'wavs');
seqSim = contains(stimuliRoot,'sequential'); %is sequential?

reconData = struct();
for trial = trials
    if seqSim
        mixed = wavs(trial).tgt_m1;
        wavs(trial).tgtLong = [wavs(trial).tgt; zeros(size(wavs(trial).m1))];
        wavs(trial).m1Long = [zeros(size(wavs(trial).tgt)); wavs(trial).m1];
        ref0 = wavs(trial).tgtLong;
        ref90 = wavs(trial).m1Long;
    else
        mixed = wavs(trial).mixed;
        ref0 = wavs(trial).tgt;
        ref90 = wavs(trial).m1;
    end
    
    s_filt.sL=ERBFilterBank(mixed(:,1),fcoefs); %freq*time
    s_filt.sR=ERBFilterBank(mixed(:,2),fcoefs);
    
    % reconstruction parameters
    params = struct();
    params.fcoefs = fcoefs;
    params.cf = cf;
    params.fs = fs;
    params.delay = 0;

    params.maskRatio = 0.5;
    params.tau = 0.02;
    masks = calcSpkMask(data(1).C_V_spikes,fs,params);

    % FR-mask reconstruction
    [rstim, maskedWav] = applyMask(masks.FR,s_filt.sL,s_filt.sR,1,'filt');
    st90 = runStoi(rstim,ref90,fs);
    st0 = runStoi(rstim,ref0,fs);
    reconData(end+1).type = 'Right';
    reconData(end).recon = rstim;
    reconData(end).tf = maskedWav;
    reconData(end).st0 = st0;
    reconData(end).st90 = st90;
    reconData(end).trial = trial;
end

%% evalute results - xcorr
tic

if seqSim
    ref(1).wav = wavs.tgtLong;
    ref(2).wav = wavs.m1Long;
else
    ref(1).wav = wavs.tgt;
    ref(2).wav = wavs.m1;
end
ref(1).label = '0 degree';
ref(2).label = '90 degree';
h1 = figure;
for i = 1:2 % 2 talkers
    s_filt.sL=ERBFilterBank(ref(i).wav(:,1),fcoefs); %freq*time
    s_filt.sR=ERBFilterBank(ref(i).wav(:,2),fcoefs);
    
    % "references"
    figure; plot_db(s_filt.sL+s_filt.sR,80)
    yticks(1:2:nf);
    yticklabels(round(cf(1:2:end)/1000,2));
    ylabel('frequency (kHz)')
    xlabel('time')
    title(ref(i).label)
    
    ref(i).tf = s_filt.sL+s_filt.sR;
    ev(i).cc = movxcorrKC(ref(i).tf,maskedWav.L+maskedWav.R,fs*0.02);
    ev(i).scc = smoothdata(ev(i).cc,'movmean',10000);
    toc
    
    figure(h1);plot(ev(i).scc,'linewidth',2); hold on;
end

%% reconstruced signal
figure; plot_db(maskedWav.L+maskedWav.R,80)
yticks(1:2:nf);
yticklabels(round(cf(1:2:end)/1000,2));
ylabel('frequency (kHz)')
xlabel('time')
title(['reconstruction, attend ' num2str(attendAz)])