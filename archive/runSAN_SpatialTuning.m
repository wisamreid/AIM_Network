% This file solves the full I, C, R, I2 model with IC input
% 
% Load IC spikes and concatenate them to simulate consecutive talkers
% Processes the concatenated spike trains through the network

addpath('dependencies')
addpath('mechs')
addpath(genpath('..\DynaSim'))

h = figure('Position',[50,50,850,690]);
set(0, 'DefaultFigureVisible', 'off')

%% load IC spikes, save to study dir
fileLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\WGNs';
sourceLocs = [-90:5:90];
fs = 44100;
for iii = [19]%1:length(sourceLocs)
    study_dir = fullfile(pwd, 'run', [mfilename 'WGNLoc' num2str(sourceLocs(iii))]);
    mkdir(fullfile(study_dir, 'solve'));
    dataFile = ls([fileLoc filesep sprintf('* %02ideg.mat',sourceLocs(iii))]);
    temp = load([fileLoc filesep dataFile]);
    copyfile([fileLoc filesep dataFile],[study_dir filesep 'IC_spks.mat']);

%% solver params
simLen = fs/1000*20; %20 ms
time_end = floor(simLen/fs*1000); % ms
spk_IC = spkTime2Train(temp.spk_IC,fs,simLen);

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
% dt = 0.025; %ms % the IC input is currently at dt=0.025
dt = 1/fs*1000;
%% neurons
% model structure
s=[];

nLocs = size(temp.spk_IC,2);
nFreqs = 64;
nCells = nLocs*nFreqs;

noise = 0.01; % low noise

% neuron populations
s.populations(1).name = 'I';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',.2, 'noise',noise}; % 10-20 Hz spiking at rest
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

rcNetconCell = repmat({ones(nLocs,1)},1,nFreqs);
rcNetconGroupByLoc = blkdiag(rcNetconCell{:});
rcNetcon = regroup(rcNetconGroupByLoc, [nFreqs,nLocs]);
%% mechanisms
s.connections(1).direction='I->I';
s.connections(1).mechanism_list='IC';
% s.connections(end).parameters={'g_postIC',0.11}; % 1000 hz spiking
s.connections(end).parameters={'g_postIC',0.015};

s.connections(end+1).direction='R->R';
s.connections(end).mechanism_list='IC';
s.connections(end).parameters={'g_postIC',0.015};

s.connections(end+1).direction='I->R';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.4, 'tauR',0.4, 'tauD',20, 'netcon',irNetcon, 'ESYN',-80}; 
%reversal potential ESYN = inhibitory

s.connections(end+1).direction='R->C';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.2, 'tauR',0.4, 'tauD',2, 'netcon',rcNetcon}; 

s.connections(end+1).direction='st->st';
s.connections(end).mechanism_list='initI2';
s.connections(end).parameters={'g_preI2',0.02,'nFreq',nFreqs,'nLocs',nLocs,'simLen',simLen+200,'fs',fs};


s.connections(end+1).direction='st->I';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.2, 'tauR',0.4, 'tauD',10, 'netcon',i2iNetcon, 'ESYN',-80}; 

%% vary
vary = {
%     'R->R', 'g_postIC', 0.02:0.01:0.05;  
    'I->I','g_postIC',0.015;
};
nVary = calcNumVary(vary);
parfor_flag = double(nVary > 1); % use parfor if multiple sims


%% simulate
tic
compile_flag = 0;
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',1, 'vary',vary, 'verbose_flag',1,...
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

%% plot all cells
set(0,'currentfigure',h);
for j = 1:length(data)
    IVspikes = logical([data(j).I_V_spikes])';
    RVspikes = logical([data(j).R_V_spikes])';
    I2Vspikes = logical([data(j).st_V_spikes])';
    for i = 1:nLocs
        idx = 1+nFreqs*(i-1):nFreqs*i;
        subplot(4,nLocs,i)
        plotSpikeRasterFs(I2Vspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        if i==1, ylabel('I2 spikes'); end
        subplot(4,nLocs,i+nLocs)
        plotSpikeRasterFs(RVspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        if i==1, ylabel('R spikes'); end
        subplot(4,nLocs,i+2*nLocs)
        plotSpikeRasterFs(IVspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        if i==1, ylabel('I spikes'); end
        subplot(4,nLocs,i+3*nLocs)
        icSpikes = logical(squeeze(spk_IC(:,:,i))'); 
        plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200]); 
        if i==1, ylabel('IC spikes'); end
    end
end

%% save
saveLoc = [fileLoc filesep 'AIMnetwork NoAttend'];
cSpikes = ([data.C_V_spikes])';
parSave([saveLoc filesep sprintf('WGN%02ideg Cspks.mat',sourceLocs(iii))],cSpikes,fs,sourceLocs,time_end)
saveas(gcf,[saveLoc filesep sprintf('WGN%02ideg Cspks.pdf',sourceLocs(iii))]);
clf;
end
%%
% addpath('PISPA2.0\plotting')
% figure;
% spksPerTime = sum(cSpikes);
% plot(spksPerTime);
% cSpksBinned = ones(1,length(spksPerTime)/(fs*0.02));
% for i = 1:length(cSpksBinned)
%     binLen = fs*0.02;
%     cSpksBinned(i) = sum(spksPerTime(1+binLen*(i-1):binLen*i));
% end
% plot(cSpksBinned)
% xlabel('time (ms)')
% 
% avgRate = sum(spksPerTime(fs*0.5:fs))/0.5/64;
% figure;plot(avgRate)

% figure;
% plotSpikeRasterFs(logical(cSpikes), 'PlotType','vertline', 'Fs',fs);
% xlim([0 time_end+200]); 