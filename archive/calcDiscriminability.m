% calculate cross correlation for all data
cd('C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork')
addpath(genpath('..\dynasim'))
addpath('eval_scripts')
addpath('dependencies')
stimuliRoot = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\';
expPrefix = '013';
tmr = '2';
stimuliLoc = [stimuliRoot 'simultaneous talkers tmr' tmr ' 64Chan500-20000hz\CRM talker4'];
saveRoot = ['data\' expPrefix '_varyI2Chan3_tmr' tmr '_simultaneous_talker']; %need to be in an existing folder
fs = 40000;
results = struct();

% load simulation data
dataSource = dir(['run' filesep 'run' expPrefix '*']);
figure;
for trial = 1:20
dataFolder = dataSource(trial).name;
data = dsImport(['run' filesep dataFolder]);
variedParam = [data.st_st_IappI2Varied];

% load simuli waveforms
talkerSet = num2str(trial,'%02i');
saveLoc = [saveRoot talkerSet];
targetName = ls(sprintf('%s*set_%s_*target_conv.wav',[stimuliLoc filesep],talkerSet));
mixedName = ls(sprintf('%s*set_%s_*mixed.wav',[stimuliLoc filesep],talkerSet));
spkName = ls(sprintf('%s*set_%s_*SpkIC.mat',[stimuliLoc filesep],talkerSet));
mask1Name = ls(sprintf('%s*set_%s_*masker1_conv.wav',[stimuliLoc filesep],talkerSet));
mask2Name = ls(sprintf('%s*set_%s_*masker2_conv.wav',[stimuliLoc filesep],talkerSet));
target = audioread([stimuliLoc filesep targetName]);
mixed = audioread([stimuliLoc filesep mixedName]);
mask1 = audioread([stimuliLoc filesep mask1Name]);
mask2 = audioread([stimuliLoc filesep mask2Name]);
ic(1).wav = target;
ic(2).wav = mask1;
ic(3).wav = mask2;
refLabels = {'target','masker1','masker2'};
options = load([stimuliLoc filesep spkName],'fcoefs','cf');
options.fs = 40000;

simLen = length(mixed);
time_end = floor(simLen/fs*1000); % ms

% Evaluate Output Intelligibilty
for dataCounter = 1:length(data)
    cSpikes = ([data(dataCounter).C_V_spikes])';
    plotSpikeRasterFs(logical(cSpikes), 'PlotType','vertline', 'Fs',fs);
    xlim([0 time_end+200])
    title(['I_{app} ' num2str(variedParam(dataCounter))])
    saveas(gcf,[saveLoc filesep 'st_Iapp' num2str(variedParam(dataCounter),'%.5f') '.tif'])
    clf;
    
    params = struct();
    params.fcoefs = options.fcoefs;
    params.cf = options.cf;
    params.fs = options.fs;
    params.spatialChan = 3;
    params.delay = 0;

    % FRmask reconstruction
    params.type = 1;
    params.maskRatio = 0;
    params.tau = 0.018;
    [st,rstim] = recon_eval(cSpikes',target,mixed,params);

    % calculate NCC to each of the stimulus
%     ref = struct();
    results(end+1).Iapp = variedParam(dataCounter);
    results(end).talkerSet = talkerSet;
    for i = 1:3
%         ref(i).tf = ERBFilterBank(ic(i).wav(:,1),options.fcoefs)+ERBFilterBank(ic(i).wav(:,2),options.fcoefs);
%         cc(i).movxcorr = movxcorrKC(ref(i).tf,rstim.tf.L+rstim.tf.R,fs*0.02);
%         cc(i).meanxcorr = mean(cc(i).movxcorr);
%         cc(i).smoothed = smoothdata(cc(i).movxcorr,'movmean',10000);
%         cc(i).meansmoothed = mean(cc(i).smoothed);
        results(end).stoi(i) = runStoi(rstim.r1d,ic(i).wav,40000,40000);
    end
%     results(end+1).cc = cc;
end
end
results(1) = []; %remove empty entry
save(['data' filesep expPrefix '_STOIs_tmr' tmr '.mat'],'results')
%% calculate discriminability
h1 = figure(1);
h2 = figure(2);
variedRange = unique([results.Iapp]);
for Iapp = 1:length(variedRange)
    % STOI-based calculation
    currentVar(Iapp).stois = reshape([results([results.Iapp]==variedRange(Iapp)).stoi],3,[]);
    [m,idx]=max(currentVar(Iapp).stois);
    disc(Iapp) = sum(idx==1)./numel(idx);
    
    meanSTOI(Iapp,:) = mean(currentVar(Iapp).stois,2);
    % cross correlation based calculation
%     currentVar(Iapp).mcc= reshape([results([results.Iapp]==variedRange(Iapp)).mcc],3,[]);
%     [m,idx]=max(currentVar(Iapp).mcc);
%     disc(Iapp) = sum(idx==1)./numel(idx);
end
set(0, 'CurrentFigure', h1)
plot(variedRange,disc)
xlabel('I_{app}')
ylabel('Discriminability')
saveas(gcf,sprintf('data\\discriminability exp%s tmr%s.pdf',expPrefix,tmr))

set(0, 'CurrentFigure', h2)
plot(variedRange,meanSTOI)
title(['TMR:' tmr])
ylabel('STOIs')
xlabel('Iapp')
saveas(gcf,sprintf('data\\STOIs exp%s tmr%s.pdf',expPrefix,tmr))

%     % plot ccs
%     h1 = figure('position',[200 200 600 300]);
%     h2 = figure('position',[200 200 600 300]);
% 
%     for i = 1:3
%         set(0, 'CurrentFigure', h1)
%         plot(ev(i).cc); hold on;
%         set(0, 'CurrentFigure', h2)
%         plot(ev(i).scc,'linewidth',2); hold on;
%     end
% 
%     set(0, 'CurrentFigure', h1)
%     legend(refLabels)
%     ylabel('NCC')
% 
%     set(0, 'CurrentFigure', h2)
%     legend(refLabels)
%     ylabel('Smoothed NCC')
%     axis tight
% 
%     figure('position',[200 200 600 300]);
%     plot_db(rstim.tf.L+rstim.tf.R,120);
