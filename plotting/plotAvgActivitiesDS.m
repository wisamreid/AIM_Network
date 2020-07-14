function plotAvgActivitiesDS(data,xaxis,faxis)
% data = dynasim output
% range = value of the varied parameter

% formatting
LineFormat.LineWidth = 1;
imgHeight = 0.8/length(data);

% get fieldds to plot
pops = {data(1).model.specification.populations.name};
fieldNames = strcat(pops,'_V_spikes');
    

for i = 1:length(pops)
    for j = 1:length(data)
        spikes(:,j) = sum(data(j).(fieldNames{i}));
    end
%     figure;
%     uimagesc(xaxis,flipud(faxis),flipud(spikes))
%     xlabel('stimulus frequency (kHz)')
%     ylabel('Model neuron frequency (kHz)')
%     title(['avg ' strrep(fieldNames{i},'_','-')])
%     set(gca,'ydir','normal')
%     caxis([0,200])
%     colorbar;
%     
    figure;
    imagesc(spikes)
    xticks(1:7:64)
    xticklabels(xaxis(1:7:64))
    yticks(1:7:64)
    yticklabels(faxis(1:7:64))
    xlabel('stimulus frequency (kHz)')
    ylabel('Model neuron frequency (kHz)')
    title(['avg ' strrep(fieldNames{i},'_','-')])
    caxis([0,150])
    colorbar;
end

end