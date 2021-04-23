% This scripts offer alternative configurations of the AIM network.
%   - params_alternateNetworkConfigs offers alternative network configurations 
%   that may also explain the oberved sharpening of receptive field
%
% @Kenny Chou, Boston Univ, 2020

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('BOSSA',pathCell)), rmpath(genpath('../BOSSA')); end
addpath('mechs')
addpath('network_params')
addpath('util')
addpath('util\plotting')
addpath('util\eval_scripts')
addpath(genpath('..\dynasim'))
clear options ICEsmall varies

% stimuli information
wgnDir = 'stimuli\WGNs';
sourceLocs = [-80:10:80];
fs = 40000;
stimDur = 0.08*fs;
attendAz = 30; % attended location
attend = 0;
mech = 1; %1,2, or 3 - used by the "alternative configurations" only

% study directory
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run WGN-' currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% -------------------- input information ---------------------
% this simulation uses white gaussian noise, spatialized to originate from 
% the specified locations. The inputs here are model neuron responses to
% white gaussian noise from those locations.
% ------------------------------------------------------------
frReduce = [];
sourceLocIdxs = 1:length(sourceLocs);
% sourceLocIdxs = randperm(length(sourceLocs))
for sourceLocIdx = sourceLocIdxs
    % prep/copy input spikes
    spkFile = ls([wgnDir filesep sprintf('*Position %02i.mat',sourceLocs(sourceLocIdx))]);
    stim = load([wgnDir filesep spkFile]);
    fs = stim.fs;
    spkTrains = spkTime2Train(stim.spk_IC,fs,stimDur);
    spkTrains = spkTrains(1:stimDur,:,:);
    
    % calculate fr & reduce dimension - here I collapse the inputs across 
    % the frequency dimension to reduce the simulation size.
    nAz = length(stim.azList);
    win = ones(0.050*fs,1); %50 ms windows
    fr = zeros(stimDur,length(stim.cf),nAz);
    for i = 1:nAz
        fr(:,:,i) = conv2(spkTrains(:,:,i),win,'same');
    end
    temp = squeeze(sum(fr,2))/nAz/0.05; %avg FR in hz
    frReduce = [frReduce; temp];
end

% new spikes based on combined FR
newTrain = zeros(size(frReduce));
for i = 1:nAz
    [newTrain(:,i),~] = spike_generator_kc(frReduce(:,i),1/fs);
end
newTrainReduced = newTrain(:,1:2:end);
azListReduced = stim.azList(1:2:end);
% figure;
% plotSpikeRasterFs(logical(newTrainReduced)','PlotType','vertline', 'Fs',fs);
% xlim([0 200*length(sourceLocs)])
    
% save modified input. The input must be saved to this folder in order
% for the network to read it.
spk_IC = [zeros(2000,length(azListReduced)); newTrainReduced];
save([study_dir filesep 'solve' filesep sprintf('IC_spks_t%02i.mat',1)],'spk_IC');

% ---------- load network parameters (input-dependent) ----------------

params_alternateNetworkConfigs;

% ---------- save file information -----------------

% name of raster plot files based on parameters
% name = expVar x vary.range
varied_param = find(cellfun(@length,{varies.range})>1);
if isempty(varied_param), varied_param = 1; end
expVar = [varies(varied_param).conxn '_' varies(varied_param).param];
expVar = strrep(expVar,'->','_');
a = cellstr(num2str(varies(varied_param).range','%.5f'));
b = {expVar};
[Bx,Ax] = ndgrid(1:numel(b),1:numel(a));
rasterFileNames = strcat(b(Bx(:)),'_',a(Ax(:)));

% store data here
expRootLoc = pwd;
expPrefix = '048 Middlebrooks cho v inh';
expName = [expPrefix ' Varying ' expVar];
dataFolder = [expRootLoc filesep 'data' filesep expName]; %things are saved to here\

% % % customized Imask for attending to different stimulus locations
% % Imask = ones(size(spk_IC));
% % for n = 1:length(sourceLocIdx)
% %     [~,attendIdx] = min(abs(stim.azList - sourceLocs(sourceLocIdx(n))));
% %     Imask(2001+stimDur*(n-1):2000+stimDur*n,attendIdx) = 0;
% % end
% % figure;imagesc(Imask')
% % save([study_dir filesep 'solve' filesep sprintf('Imask.mat')],'Imask');

% other network parameters, run network
options.expName = expName;
options.rasterFileNames = rasterFileNames;
options.saveLoc = dataFolder;
options.plotRasters = 1;
options.saveRasters = 0;
options.vizNetwork = 0;
options.fs = fs;
options.nLocs = nLocs;
options.nFreqs  = 1;
options.tgtTDchan = tgtChanIdx;
options.attend = attend;
options.constantAttendLoc = 1;

options.rcorr = 1;
options.fcoefs = stim.fcoefs;
options.cf = stim.cf;
tic;
[data,pops] = runSpatialAtten(study_dir,spk_IC,varies,options);
toc
save([study_dir filesep 'options.mat'],'options');

%% voltage plots
if length(sourceLocIdxs) == 1
taxis = (0:length(data.C_V_spikes)-1)/fs*1000; %ms

figure
subplot(3,2,1)
plot(taxis,data.TD_V(:,13),'color',[190 30 45]/255);
ylim([-80 -40])
xlim([50 130])
% xticks([])
% xticklabels([])
set(gca,'Visible','off')
title('attended channel')

subplot(3,2,3)
plot(taxis,data.E_V(:,13),'color',[0 157 220]/255);
ylim([-80 -40])
xlim([50 130])
% xticks([])
% xticklabels([])
set(gca,'Visible','off')

subplot(3,2,5)
plot(taxis,data.XI_V(:,13),'color',[100 100 100]/255);
% legend('TD','I','E');
ylim([-80 -40])
xlim([50 130])
% xticks([])
% xticklabels([])
set(gca,'Visible','off')
ylabel('Voltage Trace (mV)')


% stimuli channels
subplot(3,2,2)
plot(taxis,data.TD_V(:,2),'color',[190 30 45]/255);
ylim([-80 -40])
xlim([50 130])
% xticks([])
% xticklabels([])
set(gca,'Visible','off')
title('stimulus channel')

subplot(3,2,4)
plot(taxis,data.E_V(:,2),'color',[0 157 220]/255);
ylim([-80 -40])
xlim([50 130])
% xticks([])
% xticklabels([])
set(gca,'Visible','off')

subplot(3,2,6)
plot(taxis,data.XI_V(:,2),'color',[100 100 100]/255);
% legend('TD','I','E');
ylim([-80 -40])
xlim([50 130])
% xticks([])
% xticklabels([])
set(gca,'Visible','off')
ylabel('Voltage Trace (mV)')
end
%%
clear fr
h = ones(0.05*fs,1); % x ms filter
for n = 1:length(data)
    Cspikes = zeros(length(sourceLocs),stimDur);
%     figure;
    count = 1;
    for i = sourceLocIdxs %for each stim location
    %     [~,attendIdx] = min(abs(stim.azList - sourceLocs(i)));
        Cspikes(i,:) = data(n).C_V_spikes(2001+stimDur*(count-1):2000+stimDur*count);
        count = count+1;
%         plotSpikeRasterFs(logical(Cspikes),'PlotType','vertline', 'Fs',fs);
%         xlim([0 stimDur/fs*1000])
%         pause;
    end
    fr = conv2(Cspikes',h);
%     figure;
%     subplot('position',[0.1 0.1 0.6 0.8])
%     imagesc(fr')
%     yticks(1:length(sourceLocs))
%     yticklabels(sourceLocs)
%     title([varies(varied_param).conxn varies(varied_param).param ': ' num2str(varies(varied_param).range(n))])
%     colorbar;
%     
%     subplot('position',[0.7 0.1 0.2 0.8])
%     plot(sum(fr),17:-1:1)
%     axis tight
%     yticks([])
%     xticks([])
%     xlim([0 150000])
end
colormap jet

smoothed = imgaussfilt(fr,[50,1]);
figure
subplot('position',[0.1 0.1 0.6 0.8])
imagesc(smoothed');
xlabel('time (ms)')
colormap jet
colorbar
subplot('position',[0.7 0.1 0.2 0.8])
plot(sum(smoothed),17:-1:1)
axis tight
yticks([])
xticks([])
xlim([0 150000])
