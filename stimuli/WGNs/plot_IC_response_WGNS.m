% vizualize spiking patterns in this folder
files = ls('Spatial Tuning*');
for i = 1:size(files,1)
    load([files(i,:)])
    figure;
    suptitle(files(i,:))
    for neuron = 1:length(azList)
        if sum(~cell2mat(cellfun(@isempty,spk_IC(:,neuron),'UniformOutput',false))) == 0
            continue,
        end
        subplot(5,9,neuron)
        title([num2str(azList(neuron)) 'neuron'])
        plotSpikeRasterFs(spk_IC(:,neuron), 'PlotType','vertline', 'Fs',fs);
        xlim([0 0.150])
    end
    
end

%% response of one set of STN to various stimuli locations
files = ls('Spatial Tuning*');
[~,azIdx] = min(abs(azList+30));
spkTrain30 = zeros(length(stimuliLoc),6000);
for i = 1:size(files,1)
    load(ls(sprintf('*Position %02i.mat',stimuliLoc(i))))
    spkTrains = spkTime2Train(spk_IC,fs,0.15*fs);
    spkTrainTot = squeeze(sum(spkTrains,2));
    spkTrain30(i,:) = spkTrainTot(:,azIdx);
end
imagesc(conv2(spkTrain30',ones(500,1))')
yticks(1:length(stimuliLoc))
yticklabels(stimuliLoc)