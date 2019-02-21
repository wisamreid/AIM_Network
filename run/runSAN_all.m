% This file solves the full I, C, R model with IC input



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
copyfile( fullfile('eval_data', 'TM_00_90_set_01_SpkIC.mat'), fullfile(study_dir, 'solve', 'IC_spks.mat') );

if ~exist('spk_IC','var')
    spk_IC = load(fullfile('eval_data', 'TM_00_90_set_01_SpkIC.mat'), 'spk_IC');
    spk_IC = spk_IC.spk_IC;
    spk_IC_fs = 40e3;
    spk_IC_t = 1/spk_IC_fs:1/spk_IC_fs:size(spk_IC,1)/spk_IC_fs;
end
%% neurons
% model structure
s=[];

nFreqs = 36;
nCells = 5*nFreqs;

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
s.populations(3).size = 1;
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

irNetcon = ones(nCells)-diag(ones(1,nCells));
i2iNetcon = diag(ones(1,nCells));
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
s.connections(end).parameters={'gSYN',.1, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; 
%reversal potential ESYN = inhibitory

s.connections(end+1).direction='R->C';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.1, 'tauR',0.4, 'tauD',2, 'netcon','ones(N_pre,N_post)'}; 

s.connections(end+1).direction='st->st';
s.connections(end).mechanism_list='initI2';
s.connections(end).parameters={'g_preI2',0.01}; % 100 hz spiking


s.connections(end+1).direction='st->I';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.09, 'tauR',0.4, 'tauD',10, 'netcon',i2iNetcon, 'ESYN',-80}; 

%% vary
vary = {
    'I->R', 'gSYN', 0.1;  
};
nVary = calcNumVary(vary);
parfor_flag = double(nVary > 1); % use parfor if multiple sims


%% simulate
tic
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1/dt, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',1, 'vary',vary, 'verbose_flag',1,...
  'parfor_flag',parfor_flag);
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
temp = [data.R_V_spikes];
length(find(temp)) %number of spikes present

%%
figure;
IVspikes = logical([data.I_V_spikes])';
RVspikes = logical([data.R_V_spikes])';
I2Vspikes = logical([data.st_V_spikes])';
idx = 1:5:nCells*36;
for i = 1:5
    subplot(4,5,i)
    plotSpikeRasterFs(I2Vspikes(idx+i-1,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
    xlim([0 2000])
    if i==1, ylabel('I2 spikes'); end
    subplot(4,5,i+5)
    plotSpikeRasterFs(RVspikes(idx+i-1,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
    xlim([0 2000])
    if i==1, ylabel('R spikes'); end
    subplot(4,5,i+10)
    plotSpikeRasterFs(IVspikes(idx+i-1,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
    xlim([0 2000])
    if i==1, ylabel('I spikes'); end
    subplot(4,5,i+15)
    icSpikes = logical(squeeze(spk_IC(:,i,:))'); 
    plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',spk_IC_fs);
    xlim([0 2000])
    if i==1, ylabel('IC spikes'); end
end

%% Evaluate Output Intelligibilty
addpath('eval_scripts')
cSpikes = ([data.C_V_spikes])';
maskC = calcSpkMask(cSpikes',40000,'alpha',.2);
% maskC = calcSpkMask(squeeze(spk_IC(:,3,:)),40000,'alpha',.2);
figure;
imagesc(maskC); title('C mask');



%%
IC_info = load(fullfile('eval_data', 'TM_00_90_set_01_SpkIC.mat'), 'fcoefs','cf');
targetLoc = 'eval_data\TM_00_90_set_01_target.wav';
targetSpatializedLoc = 'eval_data\TM_00_90_set_01_target_conv.wav';
mixedLoc = 'eval_data\TM_00_90_set_01_mixed.wav';

fs = 40000;
params.fcoefs = IC_info.fcoefs;
params.cf = IC_info.cf;
params.fs = fs;
% params.frgain = 1;
% params.low_freq = 200; %min freq of the filter
% params.high_freq = 8000;
% params.numChannel = length(IC_info.cf);

[out,rstim1,rstim2,rstim3] = recon_eval(data,targetLoc,targetSpatializedLoc,mixedLoc,params);


