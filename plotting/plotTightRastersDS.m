function plotTightRastersDS(data, labels, xmax, fs)
% plot sets of spike rasters
% data = dynasim output
% range = value of the varied parameter

% formatting
LineFormat.LineWidth = 1;
imgHeight = 0.8/length(data);

% get fieldds to plot
pops = {data(1).model.specification.populations.name};
fieldNames = strcat(pops,'_V_spikes');
    

for i = 1:length(pops)
    figure;
    for j = 1:length(data)
        spikes = logical([data(j).(fieldNames{i})])';
        ypos = 0.9-imgHeight*(j);
        subplot('position',[0.1 ypos 0.8 imgHeight])
        plotSpikeRasterFs(spikes, 'PlotType','vertline', 'Fs',fs,'LineFormat',LineFormat);
        xlim([0 xmax])
        xticks([])
        yticks([])
        ylabel([num2str(labels(j))])
    end
    xticks([0:50:250])
    xticklabels([0:50:250])
    suptitle(strrep(fieldNames{i},'_','-'))
end

end