% This file solves the full I, C, R, I2 model with IC input
% 
% This file produces the exact same output as the initial implementation
% (before combining frequency and location into one dimension). It's a good
% starting point for further investigation.
% Runtime for dsSimulate is ~30s in the test system (down from ~3 min in
% the previous implementation)


%% solver params
time_end = 1.9e3; % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
dt = 0.025; %ms % the IC input is currently at dt=0.025

study_dir = fullfile(pwd, 'run', mfilename);
addpath('dependencies')
addpath('mechs')
addpath(genpath('../dynasim'))
%% copy IC spikes to study dir
if exist(study_dir, 'dir')
  rmdir(study_dir, 's');
end
mkdir(fullfile(study_dir, 'solve'));

spkLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\CISPA2.0\Data\001 IC_spk 64Chan150-8000hz\TM at 0_90 -FreqGainNorm talker4\';
spkList = ls([spkLoc '*IC.mat']);
spkFile = spkList(1,:);
spkICCopy = fullfile(study_dir, 'solve', 'IC_spks.mat');
copyfile( fullfile(spkLoc, spkFile), spkICCopy);

if ~exist('spk_IC','var')
    %spk_IC should have dimensions [time x freq chan x spatial chan]
    spk_IC = load(fullfile(spkICCopy), 'spk_IC','fcoefs','cf');
    spk_IC = spk_IC.spk_IC;
    spk_IC_fs = 40e3;
    spk_IC_t = 1/spk_IC_fs:1/spk_IC_fs:size(spk_IC,1)/spk_IC_fs;
end

%% neurons
% model structure
s=[];

nLocs = 5;
nFreqs = 64;
nCells = nLocs*nFreqs;

noise = 0.01; % low noise

% neuron populations
s.populations(1).name = 'I';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',.51, 'noise',noise}; % 10-20 Hz spiking at rest
%tonic = bias - cells spontaneous firing

s.populations(2).name='R';
s.populations(2).equations = 'chouLIF';
s.populations(2).size = nCells;
s.populations(2).parameters = {'noise',noise};

s.populations(3).name='C';
s.populations(3).equations = 'chouLIF';
s.populations(3).size = nFreqs;
s.populations(3).parameters = {'noise',noise};

s.populations(4).name='st';
s.populations(4).equations = 'chouLIF';
s.populations(4).size = nCells;
s.populations(4).parameters = {'Itonic',0, 'noise',noise};
%% connection mechanisms
% explicitly defined mechs
% synAlpha={
%   'gSYN = 1; ESYN = 0; tauD = 1; delay = 0'
%   'f(x) = (exp(-x/tauD)).*(x>0)'
%   'netcon=ones(N_pre,N_post)'
%   'synAlpha(X,t) = gSYN .* ( f(t - tspike_pre - delay) * netcon ).*(X - ESYN)'
%   '@isyn += synAlpha(V_post,t)'
%   };

synDoubleExp={
  'gSYN = 1; ESYN = 0; tauD = 4; tauR = 1; delay = 0'
  'f(x) = (exp(-x/tauD) - exp(-x/tauR)).*(x>0)'
  'netcon=ones(N_pre,N_post)'
  'synDoubleExp(X,t) = gSYN .* ( f(t - tspike_pre - delay) * netcon ).*(X - ESYN)'
  '@isyn += synDoubleExp(V_post,t)'
  };

s.mechanisms(1).name='synDoubleExp';
s.mechanisms(1).equations=synDoubleExp;

% s.mechanisms(2).name='synAlpha';
% s.mechanisms(2).equations=synAlpha;

% build I->R netcon matrix
% netcons are [N_pre,N_post]

irNetconSmall = ones(nLocs)-diag(ones(1,nLocs));
irNetconCell = repmat({irNetconSmall},1,nFreqs);
irNetconGroupByLoc = blkdiag(irNetconCell{:});
irNetcon = regroup(irNetconGroupByLoc, [nFreqs,nLocs]);

i2iNetconGroupByLoc = diag(ones(1,nCells));
i2iNetcon = regroup(i2iNetconGroupByLoc, [nFreqs,nLocs]);

rcNetconCell = repmat({ones(5,1)},1,nFreqs);
rcNetconGroupByLoc = blkdiag(rcNetconCell{:});
rcNetcon = regroup(rcNetconGroupByLoc, [nFreqs,nLocs]);
%% mechanisms
s.connections(1).direction='I->I';
s.connections(1).mechanism_list='IC';
% s.connections(end).parameters={'g_postIC',0.11}; % 1000 hz spiking
s.connections(end).parameters={'g_postIC',0.01}; % 100 hz spiking

s.connections(end+1).direction='R->R';
s.connections(end).mechanism_list='IC';
% s.connections(end).parameters={'g_postIC',0.07};
s.connections(end).parameters={'g_postIC',0.02};

s.connections(end+1).direction='I->R';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.15, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; 
%reversal potential ESYN = inhibitory

s.connections(end+1).direction='R->C';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.1, 'tauR',0.4, 'tauD',2, 'netcon',rcNetcon}; 

s.connections(end+1).direction='st->st';
s.connections(end).mechanism_list='initI2';
s.connections(end).parameters={'g_preI2',0.02,'nFreq',nFreqs}; % 100 hz spiking


s.connections(end+1).direction='st->I';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.12, 'tauR',0.4, 'tauD',10, 'netcon',i2iNetcon, 'ESYN',-80}; 

%% vary
vary = {
    'I->R', 'gSYN', 0.12:0.01:0.17;  
};
nVary = calcNumVary(vary);
parfor_flag = 0; %double(nVary > 1); % use parfor if multiple sims


%% simulate
tic
compile_flag = 1;
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',1, 'vary',vary, 'verbose_flag',1,...
  'parfor_flag',parfor_flag,'compile_flag',compile_flag);
toc

%% insert spikes
V_spike = 50;
for iData = 1:length(data)
  for pop = {'I','R','C','st'}
    pop = pop{1};
    data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
  end
end


%%
temp = [data(6).I_V_spikes];
length(find(temp)) %number of spikes present

%%
for j = 1:length(data)
    figure;
    IVspikes = logical([data(j).I_V_spikes])';
    RVspikes = logical([data(j).R_V_spikes])';
    I2Vspikes = logical([data(j).st_V_spikes])';
    for i = 1:5
        idx = 1+nFreqs*(i-1):nFreqs*i;
        subplot(4,5,i)
        plotSpikeRasterFs(I2Vspikes(idx,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
        xlim([0 2000])
        if i==1, ylabel('I2 spikes'); end
        subplot(4,5,i+5)
        plotSpikeRasterFs(RVspikes(idx,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
        xlim([0 2000])
        if i==1, ylabel('R spikes'); end
        subplot(4,5,i+10)
        plotSpikeRasterFs(IVspikes(idx,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
        xlim([0 2000])
        if i==1, ylabel('I spikes'); end
        subplot(4,5,i+15)
        icSpikes = logical(squeeze(spk_IC(:,:,i))'); 
        plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',spk_IC_fs);
        xlim([0 2000])
        if i==1, ylabel('IC spikes'); end
    end
    figure;
    icSpikes = logical(data(j).C_V_spikes)'; 
    plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',spk_IC_fs);
    xlim([0 2000]); title('C spikes')
end
%% Evaluate Output Intelligibilty
addpath('eval_scripts')
cSpikes = ([data.C_V_spikes])';
maskC = calcSpkMask(cSpikes',40000,'alpha',.2);
% maskC = calcSpkMask(squeeze(spk_IC(:,3,:)),40000,'alpha',.2);
figure;
imagesc(maskC); title('C mask');



%%
IC_info = load(fullfile(spkICCopy), 'fcoefs','cf');
targetLoc = 'eval_data\TM_00_90_set_01_target.wav';
targetSpatializedLoc = 'eval_data\TM_00_90_set_01_target_conv.wav';
mixedLoc = 'eval_data\TM_00_90_set_01_mixed.wav';

tgt = audioread(targetLoc);
targetFiltmono = ERBFilterBank(tgt,IC_info.fcoefs);
figure;plot_db(targetFiltmono,90);

fs = 40000;
params.fcoefs = IC_info.fcoefs;
params.cf = IC_info.cf;
params.fs = fs;
% params.frgain = 1;
% params.low_freq = 200; %min freq of the filter
% params.high_freq = 8000;
% params.numChannel = length(IC_info.cf);
for j = 1:length(data)
    [out(j,:),rstim1,rstim2,rstim3] = recon_eval(data(j),targetLoc,targetSpatializedLoc,mixedLoc,params);
end

