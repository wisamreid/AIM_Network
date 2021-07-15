% This scripts simulates the sharpening of receptive field during an
% attention task, as reported by Lee & Middlebrooks, 2011. We offer
% alternative explanations for this observation.
%   - params_Lee2011_inhibitionMech uses top-down selective attention to
%   accomplish the observed sharpening of receptive field
%   - params_Lee2011_ChoMech3 simulates cholinergic mechanisms that may be
%   involved in the observed sharpening of receptive field
%   - params_Lee2011_combinedMech uses both selective attention and
%   cholinergic mechanisms to explain the observed sharpening of receptive
%   field
%
% @Kenny Chou, Boston Univ
% 20200918
% 20200928 - concatenated inputs spikes from various sources into one long
%            stimulus to reduce the number of simulations required
% set up paths

pathCell = regexp(path,pathsep,'split');

if any(strcmpi('BOSSA',pathCell)), rmpath(genpath('../BOSSA')); end
addpath('dependencies')
addpath('mechs')
addpath('util')
addpath('util/plotting')
addpath('util/eval_scripts')
addpath('network_params')
addpath(genpath('../dynasim'))
clear options ICEsmall varies

% stimuli information  - update these!!
wgnDir = [cd '/stimuli/WGNs'];
sourceLocs = [-80:10:80];
fs = 40000;
stimDur = 0.08*fs; %cap stimulus duration to this length
attendAz = 30;
attend = 1; %
mech = 'both'; % 'inh','cho', or 'both'

% study directory
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run WGN-' currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% -------------------- input information ---------------------
frReduce = [];
% sourceLocIdxs = [6,16,15,13,11,14,3,12,1,17,9,2,7,10,8,4,5];
sourceLocIdxs = 1:length(sourceLocs);
% sourceLocIdxs = randperm(length(sourceLocs))
for sourceLocIdx = sourceLocIdxs%1:length(sourceLocs)
    % prep/copy input spikes
    spkFile = ls([wgnDir filesep sprintf('*Position %02i.mat',sourceLocs(sourceLocIdx))]);
    spkFile = erase(spkFile,newline);
    stim = load(spkFile);
    fs = stim.fs;
    spkTrains = spkTime2Train(stim.spk_IC,fs,stimDur);
    spkTrains = spkTrains(1:stimDur,:,:);
    
    % calculate fr & reduce dimension
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
    
% save modified input
spk_IC = [zeros(2000,length(azListReduced)); newTrainReduced];
save([study_dir filesep 'solve' filesep sprintf('IC_spks_t%02i.mat',1)],'spk_IC');

% ---------- load network parameters (input-dependent) ----------------

% params_Lee2011_inhibitionMech;
% params_Lee2011_ChoMech3;
params_Lee2011_combinedMech;

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
dataFolder = [expRootLoc filesep 'data' filesep expName]; %things are saved to here/

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

return;
%% plot tuning and marginals
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
    figure;
    subplot('position',[0.1 0.1 0.6 0.8])
    imagesc(fr')
    yticks(1:length(sourceLocs))
    yticklabels(sourceLocs)
    title([varies(varied_param).conxn varies(varied_param).param ': ' num2str(varies(varied_param).range(n))])
    colorbar;
    
    subplot('position',[0.7 0.1 0.2 0.8])
    plot(sum(fr),17:-1:1)
    axis tight
    yticks([])
    xticks([])
    xlim([0 150000])
end
colormap jet
% save([dataFolder filesep 'STRF.mat'],'sourceLocs','fr');


%% smooth figure
% [U,S,V] = svd(fr);
% nSV = 1;
% reduced = U(:,1:nSV)*S(1:nSV,1:nSV)*V(:,1:nSV)';
smoothed = imgaussfilt(fr,[50,1]);

% figure;
x = [1:5199]/fs;
% imagesc(x,azListReduced,reduced');
% xlabel('time (ms)')
% colormap jet
% caxis([0 25])

figure
subplot('position',[0.1 0.1 0.6 0.8])
imagesc(x,azListReduced,smoothed');
xlabel('time (ms)')
colormap jet
colorbar
subplot('position',[0.7 0.1 0.2 0.8])
plot(sum(smoothed),17:-1:1)
axis tight
yticks([])
xticks([])
xlim([0 150000])

return;
%% plot IC tuning curves
count = 1;
for i = sourceLocIdxs %for each stim location
    spkCounts(:,:,i) = spk_IC(2001+stimDur*(count-1):2000+stimDur*count,:);
    count = count+1;
end
tuningCurves = squeeze(sum(spkCounts));
figure; 
subplot('position',[.1 .1 .6 .6])
imagesc(tuningCurves)
xlabel('stimulus location')
ylabel('Spatial channel')
xticks(1:2:17); xticklabels(sourceLocs(1:2:17))
yticks(1:2:19); yticklabels(azListReduced(1:2:19))
subplot('position',[.1 .7 .6 .2])
plot(sum(tuningCurves,1))
yticks([]); xticks([]);
subplot('position',[.7 .1 .2 .6])
plot(sum(tuningCurves,2),19:-1:1)
yticks([]); xticks([]);
%% adjust tuning curves
figure;
adjustedTC = tuningCurves;

weights = (min(sum(newTrainReduced))./sum(newTrainReduced)).^2;
weights(9) = 1.1;
weights(10) = 1.5;
weights(11) = 1.5;
weights(12) = 1.5;

% weights = ones(1,19);
% weights(1) = 2;
% weights(2) = 3;
% weights(7) = 1.5;
% weights(5) = 1.5;
% weights(8) = 2;
% weights(9) = 2;
% weights(10) = 1.5;
% weights(11) = 2;
% weights(12) = 2;
% weights(18) = 1;

for i = 1:19
    adjustedTC(i,:) = adjustedTC(i,:)*weights(i);
end

subplot('position',[.1 .1 .6 .6])
imagesc(adjustedTC)
xlabel('stimulus location')
ylabel('Spatial channel')
xticks(1:2:17); xticklabels(sourceLocs(1:2:17))
yticks(1:2:19); yticklabels(azListReduced(1:2:19))
subplot('position',[.1 .7 .6 .2])
plot(sum(adjustedTC,1))
yticks([]); xticks([]);
subplot('position',[.7 .1 .2 .6])
plot(sum(adjustedTC,2),19:-1:1)
yticks([]); xticks([]);
%% plot I firing rate
% h = ones(0.05*fs,1); % 50 ms filter
% ISpks = data(1).XI_V_spikes;
% fr = conv2(ISpks,h);
% figure;imagesc(fr'/0.05)