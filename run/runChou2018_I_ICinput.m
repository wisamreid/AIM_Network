% this file solves just I cells

%% solver params
time_end = 1e3; % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
dt = 0.05; %ms % the IC input is currently at dt=0.025

study_dir = fullfile(pwd, 'run', mfilename);

%% copy IC spikes to study dir
if exist(study_dir, 'dir')
  rmdir(study_dir, 's');
end
mkdir(fullfile(study_dir, 'solve'));

%% neurons
% model structure
s=[];

nCells = 5;

% neuron populations
s.populations(1).name = 'I';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',0.51, 'noise',0.03};


%% mechanisms
s.connections(1).direction='I->I';
s.connections(1).mechanism_list='IC';
s.connections(1).parameters={'g_postIC',0.11};

%% vary
vary = {
  'I->I', 'g_postIC', [0 logspace(-3,-2,5) 0.07 0.11];
  'I->I', 'freqInd', 1; % only choose 1 freqInd
};
nVary = calcNumVary(vary);
parfor_flag = 0;
% parfor_flag = double(nVary > 4); % use parfor if multiple sims


%% simulate
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',1, 'vary',vary, 'verbose_flag',0,...
  'parfor_flag',parfor_flag);

%% insert spikes
V_spike = 50;
for iData = 1:length(data)
  for pop = {'I'}
    pop = pop{1};
    data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
  end
end

%% plot
dsPlot(data);
% dsPlot(data, 'plot_type','raster');

%% firing rates
% calc fr
data2 = dsCalcFR(data, 'bin_size', 100, 'bin_shift',10); % 100ms bin, 10% shift

figure

% C firing rate
sp(1) = subplot(211);
if nVary == 1
  plot(data2.time_FR, data2.I_V_spikes_FR);
  colorbar
  xlabel('Time [ms]')
  ylabel('Firing Rate [Hz]')
  title('I Firing Rate')
else
  iFR = cat(2, data2.I_V_spikes_FR);
  iFR = iFR(:,3:5:end);
  
  % plot
  imagesc(data2(1).time_FR, [], iFR');
  colorbar
  xlabel('Time [ms]')
  ylabel('g-postIC')
  title('I Firing Rate')
end
yticks(1:nVary);
yticklabels(vary{strcmp(vary(:,2),'g_postIC'), 3});


% IC input
sp(2) = subplot(212);
freqInd = vary{strcmp(vary(:,2),'freqInd'), 3};

spk_IC = load(fullfile('loadData', 'IC_spks.mat'), 'spk_IC');
spk_IC = spk_IC.spk_IC;
spk_IC_fs = 40e3;
spk_IC_t = 1/spk_IC_fs:1/spk_IC_fs:size(spk_IC,1)/spk_IC_fs;

icSpikes = logical(spk_IC(:,:,freqInd)'); % take all inputs for this freq
plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'AxHandle',sp(2), 'Fs',spk_IC_fs);
colorbar
title('IC Input')

linkaxes(sp, 'x');