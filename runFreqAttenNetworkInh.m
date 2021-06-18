function [data,s] = runFreqAttenNetworkInh(study_dir,spk_IC,varies,options)
% Simulates how attention can make neural responses more negative
% This is done by using  with-in channel inhibition.

%% solver params
fs = options.fs;
simLen = length(spk_IC);
time_end = floor(simLen/fs*1000)-3; % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
% dt = 0.025; %ms % the IC input is currently at dt=0.025
dt = 1/fs*1000;
%% neurons
% model structure
s=[];

% treat Locs like Freqs and Freqs like Locs; hacky, but it's the simplest
% way to achieve our goal
nLocs = options.nLocs;
nFreqs = options.nFreqs;
nCells = nLocs*nFreqs;

noise = 0.0; % low noise

% neuron populations
% tonic = bias - cells spontaneous firing
s.populations(1).name='IC';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'noise',0};

s.populations(end+1).name='E';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'noise',noise};

s.populations(end+1).name='C';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nFreqs;
s.populations(end).parameters = {'noise',noise};

s.populations(end+1).name = 'XI';
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',.03, 'noise',noise}; % 10-20 Hz spiking at rest
% s.populations(1).parameters = {'Itonic',.01, 'noise',0.1}; % 10-20 Hz spiking at rest

s.populations(end+1).name = 'S'; % S for sharpening, which represents w/in channel
s.populations(end).equations = 'chouLIF';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',0, 'noise',noise}; % 10-20 Hz spiking at rest

Imask = ones(1,nCells);
if options.attend, Imask(options.tgtTDchan) = 0; end
s.populations(end+1).name='TD';
s.populations(end).equations = 'LIF_Iapp';
s.populations(end).size = nCells;
s.populations(end).parameters = {'Itonic',9, 'noise',0,'Imask',Imask};


%% connection mechanisms
% explicitly defined mechs
synAlpha={
  'gSYN = 1; ESYN = 0; tau = 1; delay = 0'
  'f(x) = (1/tau)*t.*(exp(1-t/tau)).*(t>0)'
  'netcon=ones(N_pre,N_post)'
  'synAlpha(X,t) = gSYN .* sum( f(t - tspike_pre - delay) * netcon ).*(X - ESYN)'
  '@isyn += synAlpha(V_post,t)'
  };

% synDoubleExp={
%   'gSYN = 1; ESYN = 0; tauD = 4; tauR = 1; delay = 0'
%   'f(x) = (exp(-x/tauD) - exp(-x/tauR)).*(x>0)'
%   'netcon = ones(N_pre,N_post)'
%   'I(X,t) = gSYN .* ( f(t - tspike_pre - delay) * netcon ).*(X - ESYN)'
%   '@isyn += I(V_post,t)'
%   'monitor I'
%   };

% currentPulse={
%   'gCP = 1; delay = 0; ESYN = 0;'
%   'g(x) = 1.*(x < 25)'
%   'netcon = zeros(N-pre,N_post)'
%   'current_pulse(X,t) =  gCP .* ( sum(g(t - tspike_pre - delay)) * netcon ).*(X - ESYN)'
%   '@isyn += current_pulse(V_post,t)'
%   };

s.mechanisms(1).name='synAlpha';
s.mechanisms(1).equations=synAlpha;
% s.mechanisms(2).name='currentPulse';
% s.mechanisms(2).equations=currentPulse;
% s.mechanisms(2).name='synAlpha';
% s.mechanisms(2).equations=synAlpha;

%% Connectivity Matrices
ICEnetcon = eye(nFreqs);
EInetcon = eye(nFreqs);
TDInetcon = eye(nFreqs);
ECnetcon = eye(nFreqs);
EEnetcon = zeros(nFreqs);
ESnetcon = eye(nFreqs);
SEnetcon = eye(nFreqs);

IEnetcon = options.IEnetcon;
ECnetcon = options.ECnetcon;

% modeling synaptic gating
tgtChan = options.tgtTDchan;
BFChan = options.BFchan;
EEnetcon(tgtChan,BFChan) = 1; %synaptic gating from target neuron to dendrite of BF neuron


%% synapses
s.connections(1).direction='IC->IC';
s.connections(end).mechanism_list='IC';
s.connections(end).parameters={'g_postIC',10,'tauR',0.5, 'tauD',1.5};

s.connections(end+1).direction='IC->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3,'tauR',0.4, 'tauD',2, 'netcon',ICEnetcon}; 

s.connections(end+1).direction='IC->S';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3,'tauR',0.4, 'tauD',2, 'netcon',ICEnetcon}; 

s.connections(end+1).direction='E->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',0,'tauR',0.4, 'tauD',2,'netcon',EEnetcon}; 

s.connections(end+1).direction='E->XI';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',0.4, 'tauD',2, 'netcon',EInetcon}; 
 
s.connections(end+1).direction='XI->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',1, 'tauD',10, 'netcon',IEnetcon, 'ESYN',-80}; 
%reversal potential ESYN = inhibitory

s.connections(end+1).direction='S->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tau',1, 'tauD',2, 'netcon',SEnetcon, 'ESYN',-80}; 

s.connections(end+1).direction='E->C';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',0.4, 'tauD',2, 'netcon',ECnetcon}; 

% s.connections(end+1).direction='TD->TD';
% s.connections(end).mechanism_list='initI2';
% s.connections(end).parameters={'IappI2',1,'IappI2Varied',0.005,'stSubNetcon',stFreqChanTargetNetcon,'simLen',simLen+200};
% Iapp: current applied to all TD neurons;
% IappI2Varied: current applied to specific neurons; 
%               needs to be used in conjunction with stSubNetcon.

s.connections(end+1).direction='TD->XI';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',3, 'tauR',1, 'tauD',10, 'netcon',TDInetcon, 'ESYN',-80}; 

s.connections(end+1).direction='TD->E';
s.connections(end).mechanism_list='dsSynapse';
s.connections(end).parameters={'gSYN',0.8, 'tauR',1, 'tauD',10, 'netcon',TDInetcon, 'ESYN',-80}; 

if options.vizNetwork, vizNetwork(s); end
%% vary
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end
nVary = calcNumVary(vary);
parfor_flag = double(nVary > 1); % use parfor if multiple sims

%% simulate
tic
compile_flag = 0;
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',0, 'vary',vary, 'verbose_flag',0,...
  'mex_flag',1,'parfor_flag',parfor_flag,'compile_flag',compile_flag);
toc

%% insert spikes
V_spike = 30;
for iData = 1:length(data)
  for pop = {s.populations.name}
    pop = pop{1};
    data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
  end
end

%% plot each population of cells
if options.plotRasters
    set(0, 'DefaultFigureVisible', 'off')
    figure('position',[500 100 700 850]);
    time_end = floor(simLen/fs*1000);
    pops = {data(1).model.specification.populations.name};
    fieldNames = strcat(pops,'_V_spikes');
    imgHeight = 0.8/(length(pops)+1);
    axisInc = nCells/8;
    for j = 1:length(data)
        clf;
        for i = 1:length(pops)
            % dynasim model spikes
            spikes = logical([data(j).(fieldNames{i})])';
            popName = strrep(fieldNames{i},'_V_spikes','-spks');
            ypos = 0.9-imgHeight*(length(pops)+1-i);
            subplot('position',[0.1 ypos 0.6 imgHeight]);
            plotSpikeRasterFs(spikes, 'PlotType','vertline', 'Fs',fs);
            xlim([0 time_end])
            ylabel(popName)
            xticks([])
            
            % PSTH
            subplot('position',[0.755 0.1 0.2 0.875])
            plot(sum(spikes,2),1:nCells); hold on;
        end

        % ic model spikes
        if length(vary{1,3}) > 1
            load([study_dir filesep 'solve' filesep sprintf('IC_spks_t%02i',vary{1,3}(j))],'spk_IC');
        else
            load([study_dir filesep 'solve' filesep sprintf('IC_spks_t%02i',vary{1,3}(1))],'spk_IC');
        end
        icSpikes = logical(spk_IC');
        ypos = 0.9-imgHeight*(length(pops)+1);
        subplot('position',[0.1 ypos 0.6 imgHeight]);
        plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end])
        ylabel('IC model spks')
        yticks(1:axisInc:nCells)
        yticklabels(round(options.cf(1:axisInc:nCells)/1000,1))

        % PSTH
        subplot('position',[0.755 0.1 0.2 0.875])
        plot(sum(icSpikes,2),1:nCells); hold on;
        hl = legend([pops,'IC']);
        hl.Location = 'southeast';
        yticks(1:axisInc:nCells)
        yticklabels(round(options.cf(1:axisInc:nCells)/1000,1))
        ylabel('CF (kHz)')
        xlabel('spike count')
        title('PSTH')
        set(gca,'ydir','reverse')
        axis tight
    
        % save raster plot
        suptitle({});
        titleTxt = {};
        for i = 1:length(data(j).varied)
            titleTxt{i} = [data(j).varied{i} ': ' num2str(eval(sprintf('data(%i).%s',j,data(j).varied{i})))];
            tempidx = strfind(titleTxt{i},'_');
            titleTxt{i}(tempidx) = '-';
        end
        annotation('textbox',[.1 .9 .2 .1],...
               'string',titleTxt(:),...
               'FitBoxToText','on',...
               'LineStyle','none')
        saveLoc = options.saveLoc;
        saveTifName = fullfile(saveLoc,sprintf('%s.tif',options.rasterFileNames{j}));
        if ~exist(saveLoc,'dir'), mkdir(saveLoc); end
        saveas(gcf,saveTifName)
    end
    set(0, 'DefaultFigureVisible', 'on')
end

end