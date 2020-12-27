function [data,s] = runFreqAttenNetwork(study_dir,spk_IC,varies,options)
% adopted from runSAN; originally have 64 freq channel x 5 spatial channel. 
% Updated to only have 1 spatial channel.
%
% May need to edit to enable different netcons

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

noise = 0.01; % low noise

% neuron populations
s.populations(1).name = 'XI';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'Itonic',.03, 'noise',noise}; % 10-20 Hz spiking at rest
% s.populations(1).parameters = {'Itonic',.01, 'noise',0.1}; % 10-20 Hz spiking at rest
%tonic = bias - cells spontaneous firing

s.populations(2).name='E';
s.populations(2).equations = 'chouLIF';
s.populations(2).size = nCells;
s.populations(2).parameters = {'noise',noise};

s.populations(3).name='C';
s.populations(3).equations = 'chouLIF';
s.populations(3).size = nFreqs;
s.populations(3).parameters = {'noise',noise};

s.populations(4).name='TD';
s.populations(4).equations = 'chouLIF';
s.populations(4).size = nCells;
s.populations(4).parameters = {'Itonic',0, 'noise',0};
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

%% Connectivity Matrices
EInetcon = eye(nFreqs);
TDInetcon = eye(nFreqs);
ECnetcon = eye(nFreqs);

% "target channel" that activates during attention
tgtChan = options.tgtTDchan;

% local frequency convergence from E to C
paramH = options.paramH;
x = options.cf/1000; %unit in kHz
hGauss = exp(-0.5*((x-paramH.t0)/paramH.EC_BW).^2);
hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
h = hGauss.*hcos;
h = h(h>0.0001); %remove 0 values
ECnetcon = conv2(ECnetcon,h);
[~,b] = max(h);
ECnetcon = ECnetcon(b:b+nFreqs-1,:);
ECnetcon(tgtChan,15) = options.tgtTDgsyn; %extra connection from E to C cell

% I-E netcon
hGauss = exp(-0.5*((x-paramH.t0)/paramH.IE_BW).^2);
hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
h = hGauss.*hcos;
h = h(h>0.0001); %remove 0 values
IEnetcon = conv2(eye(nFreqs),h);
[~,b] = max(h);
IEnetcon = IEnetcon(b:b+nFreqs-1,:);
IEnetcon = IEnetcon - eye(nFreqs);

% target specific TD neurons to turn off
stFreqChanTargetNetcon = zeros(nFreqs);

% % E->E netcon; specifies applied current
Einitnetcon = zeros(1,nFreqs);
Einitnetcon(tgtChan) = 1;
%% mechanisms
s.connections(1).direction='E->E';
s.connections(end).mechanism_list='IC';
s.connections(end).parameters={'g_postIC',0.017,'Iapp',0,'IappNetcon',Einitnetcon};

s.connections(end+1).direction='E->XI';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.2, 'tauR',0.4, 'tauD',2, 'netcon',EInetcon}; 

s.connections(end+1).direction='XI->E';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.15, 'tauR',0.4, 'tauD',10, 'netcon',IEnetcon, 'ESYN',-80}; 
%reversal potential ESYN = inhibitory

s.connections(end+1).direction='E->C';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.2, 'tauR',0.5, 'tauD',1, 'netcon',ECnetcon}; 

s.connections(end+1).direction='TD->TD';
s.connections(end).mechanism_list='initI2';
s.connections(end).parameters={'IappI2',0.015,'IappI2Varied',0.005,'stSubNetcon',stFreqChanTargetNetcon,'simLen',simLen+200};
% Iapp: current applied to all TD neurons;
% IappI2Varied: current applied to specific neurons; 
%               needs to be used in conjunction with stSubNetcon.

s.connections(end+1).direction='TD->XI';
s.connections(end).mechanism_list='synDoubleExp';
s.connections(end).parameters={'gSYN',.2, 'tauR',0.4, 'tauD',10, 'netcon',TDInetcon, 'ESYN',-80}; 

% vizNetwork(s)
%% vary
vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end
nVary = calcNumVary(vary);
parfor_flag = 1;
% parfor_flag = double(nVary > 1); % use parfor if multiple sims


%% simulate
tic
compile_flag = 0;
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
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

if options.plotRasters
    set(0, 'DefaultFigureVisible', 'off')
    figure;
    %% plot all cells
    for j = 1:length(data)
        clf;
        IVspikes = logical([data(j).XI_V_spikes])';
        EVspikes = logical([data(j).E_V_spikes])';
        TDVspikes = logical([data(j).TD_V_spikes])';
        Cspikes = logical([data(j).C_V_spikes])';
        time_end = floor(simLen/fs*1000)/4;
        if length(vary{1,3}) > 1
            load([study_dir filesep 'solve' filesep sprintf('IC_spks_t%02i',vary{1,3}(j))],'spk_IC');
        else
            load([study_dir filesep 'solve' filesep sprintf('IC_spks_t%02i',vary{1,3}(1))],'spk_IC');
        end
        icSpikes = logical(spk_IC'); 

        % rasters
        plotRasters(Cspikes,fs,time_end,'C spikes',0.8);
        xticks([])
        plotRasters(TDVspikes,fs,time_end,'TD spikes',0.625);
        xticks([])
        plotRasters(IVspikes,fs,time_end,'XI spikes',0.45);
        xticks([])
        plotRasters(EVspikes,fs,time_end,'E spikes',0.275);
        xticks([])
        plotRasters(icSpikes,fs,time_end,'IC spikes',0.1);
        yticks(1:8:64)
        yticklabels(round(options.cf(1:8:64)/1000,1))

        % plot PSTH
        subplot('position',[0.755 0.1 0.2 0.875])
        plot(sum(Cspikes,2),1:64); hold on;
        plot(sum(IVspikes,2),1:64); hold on;
        plot(sum(EVspikes,2),1:64); hold on;
        plot(sum(icSpikes,2),1:64); hold on;
        plot(sum(TDVspikes,2),1:64); hold on;
        hl = legend({'C','XI','E','IC','TD'});
        hl.Location = 'southeast';
        yticks(1:8:64)
        yticklabels(round(options.cf(1:8:64)/1000,1))
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

function plotRasters(spks,fs,time_end,label,ypos)
    subplot('position',[0.1 ypos 0.6 0.175]);
    plotSpikeRasterFs(spks, 'PlotType','vertline', 'Fs',fs);
    xlim([0 time_end+200])
    ylabel(label)
end

end