% plot marginals Fritz2003, Atiani parameters
figure;

load('data\\041_harmonics_AtianiNetwork\spks_t20_F.mat')
plot(sum(spks),256:-1:1); hold on;

load('data\\041_harmonics_AtianiNetwork\spks_t20_M.mat')
plot(sum(spks),256:-1:1); hold on;

load('data\\041_harmonics_AtianiNetwork\ICspks_t20.mat')
plot(sum(spk_IC),256:-1:1); hold on;
axis tight
xlabel('spike count')
set(gca,'box','off')

yticks(1:16:256)
yticklabels(round(cf(256:-16:1)/1000,2))

%%
addpath('..\BOSSA\peripheral')
addpath('..\BOSSA\plotting')
load('data\\041_harmonics_AtianiNetwork\ICspks_t20.mat')
low_freq = 50; %min freq of the filter
high_freq = 8000;
numChannel = 256;
fs = 40000;
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

figure('position',[200 200 1150 550]);
subplot('position',[0.1 0.1 0.7 0.8])
plotSpikeRasterFs(logical(spk_IC'), 'PlotType','vertline', 'Fs',fs);
xlim([0 length(spk_IC)/fs*1000])
title('Passive')
xlabel('time (ms)')
ylabel('frequency (khz)')
axisinc = 1:256/16:256;
yticks(axisinc)
yticklabels(round(cf(axisinc)/1000,2))

