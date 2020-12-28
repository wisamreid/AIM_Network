% This file solves the full I, C, R, I2 model with IC input
% 
% Load IC spikes and concatenate them to simulate consecutive talkers
% Processes the concatenated spike trains through the network

addpath('dependencies')
addpath('mechs')
addpath(genpath('../dynasim'))

currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', [mfilename currentTime]);
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('PISPA2.0',pathCell)), rmpath(genpath('../PISPA2.0')); end

%% load IC spikes, save to study dir
mkdir(fullfile(study_dir, 'solve'));
fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\simultaneous talkers 64Chan500-20000hz\CRM talker4';
% fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\IC - simultaneous talkers 64Chan 200-8000Hz\TIMIT';
talkerSet = 20;

% load data
mixedName = ls([fileloc filesep sprintf('*set_%02i*mixed.wav',talkerSet)]);
mixed = audioread([fileloc filesep mixedName]);
spkName = ls([fileloc filesep sprintf('*set_%02i*SpkIC.mat',talkerSet)]);
spks = load([fileloc filesep spkName]);
if ~isfield(spks,'fs'), spks.fs = 40000; end
spk_IC = spkTime2Train(spks.spk_IC,spks.fs,length(mixed));
refNames = {'target_conv','masker1_conv','masker2_conv',};
refLabels = {'S2','S3','S1'};
for i = 1:3
    temp = ls([fileloc filesep sprintf('*set_%02i*%s.wav',talkerSet,refNames{i})]);
    ref(i).wav = audioread([fileloc filesep temp]);
end
tgt = ref(1).wav;

save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
simLen = length(spk_IC);
fs = spks.fs;
fcoefs = spks.fcoefs;
cf  = spks.cf;

%% solver params
time_end = floor(simLen/fs*1000); % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
dt = 1/fs*1000; %ms % the IC input is currently at dt=0.025

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
s.populations(1).parameters = {'Itonic',.25, 'noise',noise}; % 10-20 Hz spiking at rest
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
s.connections(end).parameters={'g_postIC',0.015}; % 100 hz spiking

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
    'I->I','g_postIC',0.01;
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
for j = 1:length(data)
    figure;
    IVspikes = logical([data(j).I_V_spikes])';
    RVspikes = logical([data(j).R_V_spikes])';
    I2Vspikes = logical([data(j).st_V_spikes])';
    for i = 1:5
        idx = 1+nFreqs*(i-1):nFreqs*i;
        subplot(4,5,i)
        plotSpikeRasterFs(I2Vspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        if i==1, ylabel('I2 spikes'); end
        subplot(4,5,i+5)
        plotSpikeRasterFs(RVspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        if i==1, ylabel('R spikes'); end
        subplot(4,5,i+10)
        plotSpikeRasterFs(IVspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        if i==1, ylabel('I spikes'); end
        subplot(4,5,i+15)
        icSpikes = logical(squeeze(spk_IC(:,:,i))'); 
        plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200]); 
        if i==1, ylabel('IC spikes'); end
    end
end

%% Evaluate Output Intelligibilty
addpath('eval_scripts')
% targetName = ls(sprintf('%s*set_%02i_*target_conv.wav',spkLoc,talkerSet))
% mixedName = ls(sprintf('%s*set_%02i_*target_conv.wav',spkLoc,talkerSet))
% target = audioread([spkLoc targetName]);
% mixed = audioread([spkLoc mixedName]);
cSpikes = ([data.C_V_spikes])';
params = struct();
params.fcoefs = fcoefs;
params.cf = cf;
params.fs = fs;
params.spatialChan = 3;
params.delay = 0;

% FRmask reconstruction
params.type = 1;
params.maskRatio = 0;
params.tau = 0.018;
[st,rstim] = recon_eval(cSpikes',tgt,mixed,params);

%% calculate NCC
for i = 1:3    
    ref(i).tf = ERBFilterBank(ref(i).wav(:,1),spks.fcoefs)+ERBFilterBank(ref(i).wav(:,2),spks.fcoefs);
    ev(i).cc = movxcorrKC(ref(i).tf,rstim.tf.L+rstim.tf.R,fs*0.02);
    ev(i).scc = smoothdata(ev(i).cc,'movmean',10000);
    toc
end

%%
h1 = figure('position',[200 200 600 300]);
h2 = figure('position',[200 200 600 300]);

for i = 1:3
    set(0, 'CurrentFigure', h1)
    plot(ev(i).cc); hold on;
    set(0, 'CurrentFigure', h2)
    plot(ev(i).scc,'linewidth',1); hold on;
    plot([0 length(ev(i).scc)],[mean(ev(i).scc) mean(ev(i).scc)])
end

set(0, 'CurrentFigure', h1)
legend(refLabels)
ylabel('NCC')

set(0, 'CurrentFigure', h2)
legend(refLabels)
ylabel('Smoothed NCC')
ylim([0,1])

figure('position',[200 200 600 300]);
plot_db(rstim.tf.L+rstim.tf.R,90);


% % % cSpikes = spk_IC(:,:,3)';
% % % maskC = calcSpkMask(cSpikes',40000,'alpha',.02);
% % % maskC = calcSpkMask(squeeze(spk_IC(:,3,:)),40000,'alpha',.2);
% figure;
% imagesc(maskC); title('C mask');
% % % 
% % % 
% % % %% reconstruction
% % % clear out
% % % IC_info = load(fullfile(spkICCopy), 'fcoefs','cf');
% % % wavList = ls([spkLoc sprintf('*%02i*.wav',talkerSet)]);
% % % targetLoc = [spkLoc strtrim(wavList(4,:))];
% % % targetSpatializedLoc = [spkLoc strtrim(wavList(5,:))];
% % % mixedLoc = [spkLoc strtrim(wavList(3,:))];
% % % 
% % % tgt = audioread(targetLoc);
% % % targetFiltmono = ERBFilterBank(tgt,IC_info.fcoefs);
% % % % figure;plot_db(targetFiltmono,90);
% % % 
% % % fs = 40000;
% % % params.fcoefs = IC_info.fcoefs;
% % % params.cf = IC_info.cf;
% % % params.fs = fs;
% % % % params.frgain = 1;
% % % % params.low_freq = 200; %min freq of the filter
% % % % params.high_freq = 8000;
% % % % params.numChannel = length(IC_info.cf);
% % % for j = 1:length(data)
% % %     [out(j,:),rstim1,rstim2,rstim3] = recon_eval(data(j),targetLoc,targetSpatializedLoc,mixedLoc,params);
% % % end
% % % 
% % % %% Plot IC v R cells
% % % figure
% % % plotNum = round(sqrt(length(data)+2));
% % % subplot(plotNum,plotNum,1); 
% % % icSpikes = logical(squeeze(spk_IC(:,:,3))'); 
% % % plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',spk_IC_fs);
% % % xlim([0 time_end]);
% % % title('0 deg IC spikes')
% % % i = 3;
% % % idx = 1+nFreqs*(i-1):nFreqs*i;
% % % for i = 1:length(data)
% % %     subplot(plotNum,plotNum,i+1)
% % %     RVspikes = logical([data(i).R_V_spikes])';
% % %     plotSpikeRasterFs(RVspikes(idx,:), 'PlotType','vertline', 'Fs',spk_IC_fs);
% % %     xlim([0 time_end])
% % %     title(sprintf('R spikes, g IC-R %d',data(i).R_R_g_postIC));
% % % end
% % % subplot(plotNum,plotNum,i+2)
% % % plot_db(targetFiltmono,90);
% % % title('Target')