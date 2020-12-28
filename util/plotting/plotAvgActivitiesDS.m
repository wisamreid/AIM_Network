function plotAvgActivitiesDS(data,xaxis,faxis,saveLoc)
% data = dynasim output
% range = value of the varied parameter

% formatting
LineFormat.LineWidth = 1;
imgHeight = 0.8/length(data);

% get fieldds to plot
pops = {data(1).model.specification.populations.name};
fieldNames = strcat(pops,'_V_spikes');
    
numTicks = 10;
numTrials = length(data);
nf = size(data(1).(fieldNames{1}),2);
inc = floor(numTrials/numTicks);
incf = floor(nf/numTicks);

xaxis = fliplr(xaxis);

figure;
for i = 1:length(pops)
    clf;
    for j = 1:length(data)
        spikes(j,:) = sum(data(j).(fieldNames{i}));
    end
    imagesc(flipud(spikes));
    ylabel('stimulus frequency (kHz)')
    yticks(1:inc:numTrials)
    yticklabels(round(xaxis(1:inc:numTrials)/1000,1))
    xlabel('Model neuron frequency (kHz)')
    xticks(mod(nf,incf):incf:nf)
    xticklabels(round(faxis(mod(nf,incf):incf:nf)/1000,1))

    title(['avg ' strrep(fieldNames{i},'_','-')])
    caxis([0,50])
    colorbar;
    saveas(gcf,[saveLoc filesep sprintf('spkCount %s.tiff',fieldNames{i})])
end

end