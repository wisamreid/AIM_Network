% This file solves the full I, C, R model with IC input

%% solver params
time_end = 1.9e3; % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
dt = 0.025; %ms % the IC input is currently at dt=0.025

study_dir = fullfile(pwd, 'run', mfilename);

%% copy IC spikes to study dir
if exist(study_dir, 'dir')
  rmdir(study_dir, 's');
end
mkdir(fullfile(study_dir, 'solve'));
copyfile( fullfile('loadData', 'IC_spks.mat'), fullfile(study_dir, 'solve', 'IC_spks.mat') );

spk_IC = load(fullfile('loadData', 'IC_spks.mat'), 'spk_IC');
spk_IC = spk_IC.spk_IC;
spk_IC_fs = 40e3;
spk_IC_t = 1/spk_IC_fs:1/spk_IC_fs:size(spk_IC,1)/spk_IC_fs;

%% neurons
% model structure
s=[];

nCells = 5;

% neuron populations
s.populations(1).name = 'I';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',.51}; % 10-20 Hz spiking at rest

s.populations(2).name='R';
s.populations(2).equations = 'chouLIF';
s.populations(2).size = nCells;

s.populations(3).name='C';
s.populations(3).equations = 'chouLIF';
s.populations(3).size = 1;
% s.populations(3).parameters = {'V_reset',-80};


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
  'synDoubleExp(X,t) = gSYN .* ( f(t - tspike_pre - delay) * netcon ).*(ESYN - X)'
  '@isyn += synDoubleExp(V_post,t)'
  };

s.mechanisms(1).name='synDoubleExp';
s.mechanisms(1).equations=synDoubleExp;

% s.mechanisms(2).name='synAlpha';
% s.mechanisms(2).equations=synAlpha;

% build I->R netcon matrix
irNetcon = zeros(nCells);
irNetcon(3,:) = 1;
irNetcon(3,3) = 0;

%% mechanisms
s.connections(1).direction='I->I';
s.connections(1).mechanism_list='IC';
% s.connections(end).parameters={'g_postIC',0.11}; % 1000 hz spiking
s.connections(end).parameters={'g_postIC',0.01}; % 100 hz spiking

s.connections(end+1).direction='R->R';
s.connections(end).mechanism_list='IC';
% s.connections(end).parameters={'g_postIC',0.07};
s.connections(end).parameters={'g_postIC',0.01};

s.connections(end+1).direction='I->R';
s.connections(end).mechanism_list='synDoubleExp';
% s.connections(end).parameters={'gSYN',.07, 'tauR',3, 'tauD',50, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',1, 'tauD',1000, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; % ds iGABAa
s.connections(end).parameters={'gSYN',.05, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; % ds iGABAa

s.connections(end+1).direction='R->C';
s.connections(end).mechanism_list='synDoubleExp';
% s.connections(end).parameters={'gSYN',.07, 'tauR',1, 'tauD',4, 'netcon','ones(N_pre,N_post)'};
% s.connections(end).parameters={'gSYN',.07, 'tauR',1, 'tauD',3, 'netcon','ones(N_pre,N_post)'};
% s.connections(end).parameters={'gSYN',.07, 'tauR',0.4, 'tauD',2, 'netcon','ones(N_pre,N_post)'}; % ds iAMPA
s.connections(end).parameters={'gSYN',.1, 'tauR',0.4, 'tauD',2, 'netcon','ones(N_pre,N_post)'}; % ds iAMPA

%% vary
vary = {
  '(I->I,R->R)', 'freqInd', 1;
%   '(I->R,R->C)', 'gSYN', logspace(-5,-1, 5);

%   'I->R', 'gSYN', logspace(-5,0, 5);
%   'I->R', 'gSYN', logspace(-2,-1, 10);
  'I->R', 'gSYN', linspace(0.01,0.1, 10);
  
%   'R->C', 'gSYN', logspace(-2, 0, 5);
%   'R->C', 'gSYN', linspace(.005, 0.05, 10);
};
nVary = calcNumVary(vary);
parfor_flag = double(nVary >= 4); % use parfor if multiple sims


%% simulate
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',1, 'vary',vary, 'verbose_flag',0,...
  'parfor_flag',parfor_flag);

%% insert spikes
V_spike = 50;
for iData = 1:length(data)
  for pop = {'I','R','C'}
    pop = pop{1};
    data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
  end
end

%% plot
dsPlot(data);
% dsPlot(data, 'variable','I_V');
% dsPlot(data, 'variable','R_V');
% dsPlot(data, 'variable','C_V');
% dsPlot(data(1), 'plot_type','raster');

%% firing rates
% calc fr
data2 = dsCalcFR(data, 'bin_size', 100, 'bin_shift',20); % 100ms bin, 20% shift

figure

% C firing rate
sp(1) = subplot(211);
if nVary == 1
  plot(data2.time_FR, data2.C_V_spikes_FR);
  xlabel('Time [ms]')
  ylabel('Firing Rate [Hz]')
  title('C Firing Rate')
else
  cFR = cat(2, data2.C_V_spikes_FR);
  
  % plot
  imagesc(data2(1).time_FR, [], cFR');
  colorbar
  xlabel('Time [ms]')
  ylabel('Sims')
  title('C Firing Rate')
end

% IC input
sp(2) = subplot(212);
freqInd = vary{strcmp(vary(:,2),'freqInd'), 3};
icSpikes = logical(spk_IC(:,:,freqInd)'); % take all inputs for this freq
plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'AxHandle',sp(2), 'Fs',spk_IC_fs);

if nVary > 1
  colorbar
end
title('IC Input')

linkaxes(sp, 'x');