fileLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\WGNs\AIMnetwork SubICPop NoAttend';
sourceLocs = [-90:5:90];
fs = 44100;
binTime = 1; %ms
binLen = round(fs/1000*binTime);
set(0, 'DefaultFigureVisible', 'on')
clear compiled
k = ones(floor(fs/1000*1),1); %5ms window
for iii = 1:length(sourceLocs)
    file = ls([fileLoc filesep sprintf('WGN%02ideg*.mat',sourceLocs(iii))]);
    load([fileLoc filesep file]);
    totalSpks = sum(cSpikes);
    fr = conv(totalSpks,k,'same');
%     for i = 1:floor(length(totalSpks)/binLen)
%         binnedSpks(i) = sum(totalSpks(i+binLen*(i-1):binLen*i));
%     end
%     compiled(iii,:) = binnedSpks;
    compiled(iii,:) = fr;
end

figure; imagesc(1:size(compiled,2),sourceLocs,compiled)
colormap jet
padded = [zeros(length(sourceLocs),fs/2) compiled zeros(length(sourceLocs),fs)];
t = (0:fs*2)/fs*1000;
figure; imagesc(t,sourceLocs,padded)
colormap jet
figure; plot(sum(padded,2))
ylim([0 4E6])
%%
figure;
simLen = 2*fs;
time_end = floor(simLen/fs*1000); % ms
plotSpikeRasterFs(logical(cSpikes), 'PlotType','vertline', 'Fs',fs);
xlim([0 time_end+200]); 

%%
figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(squeeze(FR(:,:,i)))
end