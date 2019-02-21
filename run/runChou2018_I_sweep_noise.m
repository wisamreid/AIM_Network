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

% neuron populations
s.populations(1).name = 'I';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = 1;
s.populations(1).parameters = {'Itonic',0.51};


%% mechanisms
% s.connections(1).direction='I->I';
% s.connections(1).mechanism_list='IC';
% s.connections(end).parameters={'g_postIC',0.11};

%% vary
vary = {
  'I', 'noise', [0 linspace(0.01, 0.1,10)];
};
nVary = calcNumVary(vary);
% parfor_flag = 0;
parfor_flag = double(nVary > 4); % use parfor if multiple sims


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
data2 = dsCalcFR(data, 'bin_size', 100, 'bin_shift',20); % 100ms bin, 20% shift

figure

% C firing rate
if nVary == 1
  plot(data2.time_FR, data2.I_V_spikes_FR);
  colorbar
  xlabel('Time [ms]')
  ylabel('Firing Rate [Hz]')
  title('I Firing Rate')
else
  iFR = cat(2, data2.I_V_spikes_FR);
  
  % plot
  imagesc(data2(1).time_FR, [], iFR');
  colorbar
  xlabel('Time [ms]')
  ylabel('Sims')
  title('I Firing Rate')
end