% demo: passive switch - use expPrefix '015', simIDs = 3, presentation ='seq'
% demo: tonic firing of SI neurons - use expPrefix '016', presentation ='seq'
% demo: tonic firing of SI neurons - use expPrefix '016', presentation ='sim'
addpath(genpath('../Dynasim'))
addpath('../PISPA2.0/plotting')
expPrefix = '015';
talker = '04';
presentation = 'sim';
dataSource = dir(['run' filesep 'run' expPrefix '*_' presentation '_*talker' talker '*']);
fs = 40000; %crm fs
nFreqs = 64;
for j = 3%1:length(dataSource)
    data = dsImport(['run' filesep dataSource(j).name]);
    load(['run' filesep dataSource(j).name filesep 'options.mat'])
    icSpikeStruct = load(['run' filesep dataSource(j).name filesep 'solve' filesep 'IC_spks.mat']);
    spk_IC = icSpikeStruct.spk_IC;
    for k = 1:length(data)
    simLen = size(data(k).C_V,1);
    time_end = floor(simLen/fs*1000);
    
    figure;
    % output spikes
    subplot(5,3,2);
    cSpikes = ([data(k).C_V_spikes])';
    plotSpikeRasterFs(logical(cSpikes), 'PlotType','vertline', 'Fs',fs);
    xlim([0 time_end+200])
    xticks([]); yticks([]);
    ylabel('C')
%     title(['title placeholder'])
    
    IVspikes = logical([data(k).I_V_spikes])';
    RVspikes = logical([data(k).R_V_spikes])';
    I2Vspikes = logical([data(k).st_V_spikes])';
    channels = [1,3,5];
    for i = 1:3 %number of chennels being plotted
        spatialChan = channels(i);
        idx = 1+nFreqs*(spatialChan-1):nFreqs*spatialChan;
%         if spatialChan == 3
            subplot(5,3,i+3)
            plotSpikeRasterFs(I2Vspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
            xlim([0 time_end+200])
            xticks([]); yticks([]);
            if i==1, ylabel('TD spikes'); end
%         end
        subplot(5,3,i+6)
        plotSpikeRasterFs(RVspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200])
        xticks([]); yticks([]);
        if i==1, ylabel('R spikes'); end
%         if spatialChan == 3
            subplot(5,3,i+9)
            plotSpikeRasterFs(IVspikes(idx,:), 'PlotType','vertline', 'Fs',fs);
            xlim([0 time_end+200])
            xticks([]); yticks([]);
            if i==1, ylabel('SI spikes'); end
%         end
        subplot(5,3,i+12)
        icSpikes = logical(squeeze(spk_IC(:,:,spatialChan))'); 
        plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',fs);
        xlim([0 time_end+200]); 
        xticks([]); yticks([]);
        if i==1, ylabel('IC spikes'); end
    end
    h = gcf;
    h.Renderer='Painters'; % so that pdfs will be saved as vectors
%     saveFileName = sprintf('%s\\%s_%0.4f.tif',options.saveLoc,data(k).varied{2},data(k).st_st_IappI2Varied);
%     saveas(h,saveFileName)
    end
end

% figure;
% channels = [1,3,5];
% for i = 1:3
%     spatialChan = channels(i);
%     subplot(3,1,i)
%     icSpikes = logical(squeeze(spk_IC(:,:,spatialChan))'); 
%     plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',fs);
%     xlim([0 time_end+200]); 
%     xticks([]); yticks([]);
%     if i==1, title('IC spikes'); end
% end