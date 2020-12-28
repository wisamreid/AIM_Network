 function [data,s] = runSpatialAtten(study_dir,spk_IC,varies,options)
%% solver params
fs = options.fs;
if isfield(options,'simLen')
    simLen = options.simLen;
else
    simLen = length(spk_IC);
end
time_end = floor(simLen/fs*1000); % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
% dt = 0.025; %ms % the IC input is currently at dt=0.025
dt = 1/fs*1000;
%% neurons
% model structure
s=[];

nLocs = options.nLocs;
nFreqs = options.nFreqs;
nCells = nLocs*nFreqs;

noise = 0.01; % low noise

% neuron populations
s.populations(1).name='IC';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'noise',0};

s.populations(end+1).name='E';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'noise',noise,'Itonic',0};

s.populations(end+1).name='C';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nFreqs;
s.populations(end).parameters = {'noise',noise};

s.populations(end+1).name = 'XI';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',3 'noise',noise,'ton',2,'toff',time_end};
%tonic = bias - cells spontaneous firing

if options.constantAttendLoc
    ImaskSmall = options.TD_Imask;
    Imask = repmat(ImaskSmall,nFreqs,1);
    Imask = regroup(Imask,[nFreqs,nLocs]);
    s.populations(end+1).name='TD';
    s.populations(end).equations = 'LIF_Iapp';
    s.populations(end).size = nCells;
    s.populations(end).parameters = {'Itonic',8, 'noise',0,'Imask',Imask'};
else
    s.populations(end+1).name='TD';
    s.populations(end).equations = 'LIF_variable_Imask';
    s.populations(end).size = nCells;
    s.populations(end).parameters = {'Itonic',8, 'noise',0};
end
%% connection mechanisms - in mech files

%% Connection matrices; netcons are [N_pre,N_post]
% each frequency channel is independent of one another
% small netcons are netcons for each frequency;

if isfield(options,'ICEnetcon')
    ICEsmall = options.ICEnetcon;
else
    ICEsmall = eye(nLocs);
end
ICEnetcon = extendAndRegroup(ICEsmall,nFreqs,nLocs); 

ECsmall = options.ECnetcon;
ECnetcon = extendAndRegroup(ECsmall,nFreqs,nLocs);

EIsmall = eye(nLocs,nLocs);
EInetcon = extendAndRegroup(EIsmall,nFreqs,nLocs);

IEsmall = ones(nLocs,nLocs) - eye(nLocs);
IEnetcon = extendAndRegroup(IEsmall,nFreqs,nLocs);

%% mechanisms
s.connections(1).direction='IC->IC';
s.connections(end).mechanism_list='IC';
s.connections(end).parameters={'g_postIC',10,'tauR',0.5, 'tauD',1.5};

s.connections(end+1).direction='IC->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3,'tauR',0.4, 'tauD',2,'netcon',ICEnetcon}; 

s.connections(end+1).direction='E->XI';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',0.4, 'tauD',2, 'netcon',EInetcon}; 

s.connections(end+1).direction='XI->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',4, 'tauR',1, 'tauD',10, 'netcon',IEnetcon, 'ESYN',-80}; 
%reversal potential ESYN = inhibitory

s.connections(end+1).direction='E->C';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',0.4, 'tauD',2, 'netcon',ECnetcon}; 

s.connections(end+1).direction='TD->XI';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',1, 'tauD',10,'netcon',eye(nCells),'ESYN',-80}; 

s.connections(end+1).direction='TD->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',4, 'tauR',1, 'tauD',10,'netcon',eye(nCells),'ESYN',-80}; 

if options.vizNetwork, vizNetwork(s); end

%% vary
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end
nVary = calcNumVary(vary);
parfor_flag = 0;
% parfor_flag = double(nVary > 1); % use parfor if multiple sims


%% simulate
tic
compile_flag = 0;
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',0, 'save_results_flag',0,...
  'study_dir',study_dir, 'debug_flag',0, 'vary',vary, 'verbose_flag',0,...
  'parfor_flag',parfor_flag,'compile_flag',compile_flag);
toc

%% insert spikes
V_spike = 50;
for iData = 1:length(data)
  for pop = {s.populations.name}
    pop = pop{1};
    data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
  end
end

%% plot
if options.plotRasters
    plotNetworkRasters;
end