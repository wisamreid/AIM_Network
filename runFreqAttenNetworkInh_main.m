% This script implements the AIM network that simulates spatial attention
% in the frequency domain, reported in [1] and [2].
%
% Load the appropriate experiment parameters by commenting/uncommenting the
% appropriate lines.
%   params_Atiani2009_3a simulates figure 3a in [1]
%   params_Atiani2009_3b simulates figure 3b in [1]
%   params_Fritz2003 simulates figure 2a in [2]
%   params_Fritz2003_Atiani simulates figure 2a in [2], using parameters
%       from simulating [1]
%   _EE0 simulations removes the extra E-E connection in the AIM network.
%
% Turn the attention mechanism on/off by setting attention = true/false
%
%
% note: it is recommended to set study_dir on a on the local machine,
% otherwise saving data will take forever.
%
% @Kenny F. Chou, Boston Univ, 2020
%
% references:
%   [1] Atiani, S., Elhilali, M., David, S. V., Fritz, J. B. & Shamma, S. A. 
%       Task Difficulty and Performance Induce Diverse Adaptive Patterns in 
%       Gain and Shape of Primary Auditory Cortical Receptive Fields. 
%       Neuron 61, 467¡V480 (2009).
%   [2] Fritz, J., Shamma, S., Elhilali, M. & Klein, D. Rapid task-related
%       plasticity of spectrotemporal receptive fields in primary auditory 
%       cortex. Nat. Neurosci. 6, 1216¡V1223 (2003).

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('BOSSA',pathCell)), rmpath(genpath('../BOSSA')); end
addpath('mechs')
addpath('network_params')
addpath('util')
addpath('util/plotting')
addpath('util/eval_scripts')
addpath(genpath('../dynasim'))
clear datetime

% ========= load experiment parameters ===========%
expPrefix = '053g EE0'; %for logging results
expName = 'Inhibitory strengthening';
attention = true; %<-----------

params_Inh;

% =================== end params ==================

% find varied parameter
varied_param = find(cellfun(@length,{varies.range})>1);
if isempty(varied_param), varied_param = 1; end
expVar = [varies(varied_param).conxn '_' varies(varied_param).param];
expVar = strrep(expVar,'->','_');

% network structure:
nFreqs = numChannel;
nLocs = 1;
if attention
    tgtI2Chan = {[tgtTDidx,1]}; % I2 channels to turn off; pairs of [freq,loc] channels
else
    tgtI2Chan = {[]};
end

% name of raster plot files based on parameters
% name = expVar x vary.range
a = cellstr(num2str(varies(varied_param).range','%.04f'));
b = {expVar};
[Bx,Ax] = ndgrid(1:numel(b),1:numel(a));
rasterFileNames = strcat(b(Bx(:)),'_',a(Ax(:)))
    
% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd-HHmmss'));
study_dir = fullfile(pwd, 'run', ['run' expName currentTime]);
mkdir(fullfile(study_dir, 'solve'));

expFolderName = sprintf('%s_%s_Attend%i',expPrefix,expName,attention);
expRootLoc = pwd; %'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network';
dataFolder = fullfile(expRootLoc,'data',expFolderName); %things are saved to here

% prep IC inputs
% ------------ load spk_IC here -----------------
for trial = trials
    spks = load(fullfile(stimuliRoot,sprintf('IC_spks_tones_%02i.mat',trial)));
    spk_IC = spkTime2Train(spks.spk_IC,spks.fs,0.3*fs);
    spk_IC = squeeze(spk_IC(:,:,3));
    spk_IC(1:400,:) = 0; %remove artifact & have a silent period
    spk_IC(1044:end,:) = 0; %truncate tone to 35 ms
    spk_IC(2867+400:end,:) = []; %truncate stimuli to 65 ms
%     plotSpikeRasterFs(logical(spk_IC'), 'PlotType','vertline', 'Fs',fs);
%     xlim([0 300])
    fs = spks.fs;
    fcoefs = spks.fcoefs;
    cf  = spks.cf;
%     spk_IC = zeros(fs*0.160,64);
%     spk_IC(1000:fs*0.0018:end,tgtTDchan-2:tgtTDchan+1) = 1;
%     spk_IC(1000:fs*0.005:end,tgtTDchan+2:tgtTDchan+5) = 1;
%     sum(spk_IC)
    simLen = length(spk_IC);
    save(fullfile(study_dir, 'solve', sprintf('IC_spks_t%02i.mat',trial)),'spk_IC')
end
% prep I2 input
addpath('mechs');
genI2input([nFreqs,nLocs],tgtI2Chan,simLen,fs,study_dir);

% run network
options.expName = expName;
options.rasterFileNames = rasterFileNames;
options.saveLoc = dataFolder;
options.plotRasters = 1;
options.vizNetwork = 1;
options.nLocs = nLocs;
options.nFreqs = nFreqs;
options.fs = fs;
options.ECnetcon = ECnetcon;
options.IEnetcon = IEnetcon;

options.tgtECgsyn = tgtECgsyn;
options.tgtTDchan = tgtTDidx;
options.BFchan = BFchanidx;
options.attend = attention;

options.rcorr = 1;
options.fcoefs = fcoefs;
options.cf = cf;

tic;
[data,pops] = runFreqAttenNetworkInh(study_dir,spk_IC,varies,options);
toc;

% save parameters & model
model = data(1).model;
save([dataFolder filesep 'options.mat'],'options','model');

%%
% targetChan = 6.8; % kHz
% [~, tgtidx]=min(abs(f0-targetChan*1000)); 
% BFChan = 8.7; % kHz
% [~, BFidx]=min(abs(f0-BFChan*1000)); 
% plotIVTraces_IC_E(data,spk_IC,fs,59,dataFolder,varies,varied_param)
% plotIVTraces_E_E(data(tgtidx),fs,59,42,dataFolder,varies,varied_param)
% 
% for j = 42:59
% plotIVTraces_E_C(data(tgtidx),fs,j,dataFolder,varies,varied_param)
% end
% % plotIVTraces(data,spk_IC,fs,26,dataFolder,varies,varied_param)
% 
% % figure;
% % subplot(2,1,2);plot(t,data.C_E_chouSynapse_I(:,15)); title('E-C current')
% % subplot(2,1,1);plot(t,data.C_V(:,15)); title('C_V')
% % 
% % figure;
% % subplot(3,1,2);plot(t,data.C_E_chouSynapse_I(:,15)); title('iC')
% % subplot(3,1,1);plot(t,data.C_V(:,15)); title('C_V')
% % subplot(3,1,3);plot(t,data.E_V(:,15)); title('E_V')
% 
% % vizNetwork(pops)
% % suptitle([num2str(paramH.BW) ' octaves'])
% 
% if length(trials) < 2
%     return;
% end
%% plot additional results
% addpath('plotting')
% range = 41:60;
% subdata = data(range);
% time_end = floor(simLen/fs*1000);
% 
% plotTightRastersDS(subdata,round(cf(range)/1000,1),200,fs)

range = trials;
faxis = cf;
taxis = subf0;
plotAvgActivitiesDS(data,taxis,faxis,dataFolder);
%% plot receptive field, centered around mean spontaneous FR

% get spiking responses of target neuron
% bestFreq = 4.2; % kHz
[~, idx]=min(abs(cf-targetFreq*1000)); 
% [~, idx]=min(abs(cf-bestFreq*1000)); 
dur = 25; %ms, bin duration
binsize = ceil(5/(1000/fs));
% gaussian = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma.^2));
% kernel = gaussian(0:1000/fs:20,10,1);

clear combinedSpks runningFR
runningFR = zeros(numChannel,length(data(1).time)-binsize);
for trial = 1:length(trials)
    tgtNeuron = data(trial).C_V_spikes(:,idx);     
    
    % calculate FR by calculating a running average
    for i = 1:length(tgtNeuron) - binsize
        runningFR(trials(trial),i) = sum(tgtNeuron(i:i+binsize));
    end
%     runningFR(trial,:) = conv(tgtNeuron,kernel,'same');
end
runningFR = runningFR/0.025; %convert to spikes/sec

% plots
figure
set(gcf,'Position',[300 300 450 500])
subplot('position',[0.1 0.1 0.6 0.8])
t0 = 400;
x = (-400:length(tgtNeuron)-44-400)/fs*1000;
yMinMax = [min(ySpacing) max(ySpacing)];

[~, minf0idx]=min(abs(f0-minf0)); 
[~, maxf0idx]=min(abs(f0-6700));
meanFR = mean(mean(runningFR(minf0idx:maxf0idx,:)));
imagesc(x,log2(f0/1000),runningFR-meanFR);
ylabel('pure tone frequency kHz')
xlabel('time (ms)')
title(['neuron with cf ' num2str(cf(idx))])
set(gca,'ydir','normal')
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
colormap jet
caxis([-150 150])

subplot('position',[0.75 0.1 0.2 0.8])
plot(sum(runningFR,2),log2(f0/1000))
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
xlabel('spike count')
titleTxt = {};
for i = 1:length(data(1).varied)
    titleTxt{i} = [data(1).varied{i} ': ' num2str(eval(sprintf('data(%i).%s',1,data(1).varied{i})))];
    tempidx = strfind(titleTxt{i},'_');
    titleTxt{i}(tempidx) = '-';
end
annotation('textbox',[.45 .8 .2 .1],...
       'string',titleTxt,...
       'FitBoxToText','on',...
       'LineStyle','none')

saveas(gcf,[dataFolder filesep 'STRF.tif'])
save([dataFolder filesep 'STRF.mat'],'runningFR','f0','x','ySpacing','bestFreq','titleTxt')
