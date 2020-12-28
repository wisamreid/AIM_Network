% plots many rasters together
% with each row = one layer
% each column = one spatial channel
% rows within each raster = frequency channels

% ------------ network properties --------------
nLocs = options.nLocs;
nFreqs = options.nFreqs;
fs = options.fs;
simLen = length(spk_IC);
time_end = floor(simLen/fs*1000); % ms

pops = {data(1).model.specification.populations.name};
popSizes = [data(1).model.specification.populations.size];
nPops = numel(pops);
fieldNames = strcat(pops,'_V_spikes');
    
% ------------ figure properties ------------------
% set(0, 'DefaultFigureVisible', 'off')
figure('position',[300 50 1300 900]);
xpos0 = 0.1;
ypos0 = 1-0.1;
imgHeight = 0.85/length(pops); % 80% of figure
imgWidth = 0.85/nLocs;
axisInc = nFreqs/8;

% for each varied parameter
for dataNum = 1:length(data)
    clf;
    
    % for each layer of neurons
    for popNum = 1:length(pops)
        spikes = logical([data(dataNum).(fieldNames{popNum})])';
        popName = strrep(fieldNames{popNum},'_V_spikes','-spks');
        ypos = ypos0-imgHeight*(length(pops)+1-popNum);

        % for each location
        for locIdx = 1:popSizes(popNum)/nFreqs
            xpos = xpos0+imgWidth*(locIdx-1);
            subplot('position',[xpos ypos imgWidth*0.95 imgHeight*0.95])
            idx = 1+nFreqs*(locIdx-1):nFreqs*locIdx;
            if isvector(spikes(idx,:))
                plot(spikes(idx,:))
            else
                plotSpikeRasterFs(logical(spikes(idx,:)), 'PlotType','vertline', 'Fs',fs);
                xlim([0 time_end+200])
            end
            if locIdx == 1, ylabel(popName); end
            if popNum == 1, xlabel(num2str(locIdx)); end
            xticks([]); yticks([]);
        end
    end
    
    % save raster plot
    if options.saveRasters
        suptitle({});
        titleTxt = {};
        for i = 1:length(data(dataNum).varied)
            titleTxt{i} = [data(dataNum).varied{i} ': ' num2str(eval(sprintf('data(%i).%s',dataNum,data(dataNum).varied{i})))];
            tempidx = strfind(titleTxt{i},'_');
            titleTxt{i}(tempidx) = '-';
        end
        annotation('textbox',[.1 .9 .2 .1],...
               'string',titleTxt(:),...
               'FitBoxToText','on',...
               'LineStyle','none')
        saveLoc = options.saveLoc;
        saveTifName = fullfile(saveLoc,sprintf('%s.tif',options.rasterFileNames{dataNum}));
        if ~exist(saveLoc,'dir'), mkdir(saveLoc); end
        saveas(gcf,saveTifName)
    end
end
%     set(0, 'DefaultFigureVisible', 'on')
