% runSAN_sims

% in a consecutive talker scenario: concatenate .wav files and spike files,
% then run the network

% note: it is not recommended to set study_dir on a different network (i.e.
% research drive) from the local machine. Saving run data will take
% forever.

% @Kenny Chou, Boston Univ
% 2019 September - created
% 2020-02-05 - combined Consecutive/simultaneous stimuli cases

% to do:
% - move stimuli prep to a new file

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('BOSSA',pathCell)), rmpath(genpath('../BOSSA')); end
addpath('dependencies')
addpath('mechs')
addpath('plotting')
addpath('..\BOSSA\plotting')
addpath('..\BOSSA\peripheral')
addpath('eval_scripts')
addpath(genpath('..\dynasim'))

stimuliRoot = 'stimuli\PureTones';
expPrefix = '018';
fs = 44100;

% ================ parameters ==============
trials = 1:64;
tgtTDchan = 23;
tgtTDgsyn = 0.5;

expName = sprintf('%s_PureTones_PassiveSTRF',expPrefix);
expVar = 'trial';
% expVar = 'i-itonic';
expRootLoc = pwd; %'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network';
dataFolder = [expRootLoc filesep 'data' filesep expName]; %things are saved to here

% network parameters
varies = struct();
varies(1).conxn = 'E->E';
varies(end).param = 'trial';
varies(end).range = trials;

varies(end+1).conxn = 'E->E';
varies(end).param = 'Iapp';
varies(end).range = -0.5;

% network structure:
nFreqs = 64;
nLocs = 1;
tgtI2Chan = {[tgtTDchan,1]}; % I2 channels to turn off; pairs of [freq,loc] channels
% tgtI2Chan = {[]};

% frequency convergence onto C neurons, for EC netcon
paramH.EC_BW= 0.25; %bandwidth, octaves
paramH.BTM= 0.01 ; %BTM, modulation
paramH.t0= 8.9; % f0, peak f
paramH.phase= 0; % phase

% x-frequency inhibition from I to E neurons
paramH.IE_BW= 1; %bandwidth, octaves
paramH.BTM= 0.01 ; %BTM, modulation
paramH.t0= 2.9; % f0, peak f
paramH.phase= 0; % phase

% name of raster plot files based on parameters
% name = expVar x vary.range
a = cellstr(num2str(varies(1).range','%.03f'));
b = {expVar};
[Bx,Ax] = ndgrid(1:numel(b),1:numel(a));
rasterFileNames = strcat(b(Bx(:)),'_',a(Ax(:)))

% =================== end params ==================
    
% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run' expName currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% prep IC inputs
% ------------ load spk_IC here -----------------
for trial = trials
    spks = load([stimuliRoot filesep sprintf('IC_spks_tones_%02i.mat',trial)]);
    spk_IC = spkTime2Train(spks.spk_IC,spks.fs,0.3*fs);
    spk_IC = squeeze(spk_IC(:,:,3));
    simLen = length(spk_IC);
    fs = spks.fs;
    fcoefs = spks.fcoefs;
    cf  = spks.cf;
    save(fullfile(study_dir, 'solve', sprintf('IC_spks_t%02i.mat',trial)),'spk_IC')
end
% prep I2 input
addpath('mechs');
genI2input([nFreqs,nLocs],tgtI2Chan,simLen,fs,study_dir);

% run network
options.expName = expName;
options.rasterFileNames = rasterFileNames;
options.saveLoc = dataFolder;
options.plotRasters = 1;
options.nLocs = nLocs;
options.nFreqs = nFreqs;
options.fs = fs;
options.paramH = paramH;

options.tgtTDgsyn = tgtTDgsyn;
options.tgtTDchan = tgtTDchan;

options.rcorr = 1;
options.fcoefs = fcoefs;
options.cf = cf;

tic;
[data,pops] = runFreqAttenNetwork(study_dir,spk_IC,varies,options);
save([study_dir filesep 'options.mat'],'options');
toc

% vizNetwork(pops)
% suptitle([num2str(paramH.BW) ' octaves'])
%% plot additional results
% addpath('plotting')
% range = 41:60;
% subdata = data(range);
% time_end = floor(simLen/fs*1000);
% 
% plotTightRastersDS(subdata,round(cf(range)/1000,1),200,fs)

range = trials;
f0 = logspace(log10(500),log10(16000),64);
faxis = round(cf(range)/1000,2);
taxis = round(f0/1000,2);
plotAvgActivitiesDS(data,taxis,faxis);
%% plot STRF
f0 = logspace(log10(500),log10(16000),64); %hz

% get spiking responses of target neuron
targetFreq = 7.5; % kHz
clear combinedSpks
[~, idx]=min(abs(cf-targetFreq*1000)); 
dur = 5; %ms, bin duration
binsize = ceil(5/(1000/fs));
for trial = trials
    tgtNeuron = data(trial).C_V_spikes(:,idx);     
    % bin responses
    i = 1;
    while binsize+binsize*(i-1) < length(tgtNeuron)
        combinedSpks(trial,i) = sum(tgtNeuron(1+binsize*(i-1):binsize+binsize*(i-1)));
        i = i+1;
    end
end
figure
x = (1:31)*dur;
imagesc(x,log2(f0/1000),combinedSpks)
ylabel('pure tone frequency kHz')
xlabel('time (ms)')
title(['neuron with cf ' num2str(cf(idx))])
set(gca,'ydir','normal')
ylim([1 4])
yticks([1:4])
yticklabels(2.^[1:4])
caxis([0,7])
colorbar;
