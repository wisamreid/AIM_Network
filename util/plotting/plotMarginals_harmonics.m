% a part of Figure 7A
% plot marginals Fritz2003, Atiani parameters
figure;

load('data\\041_harmonics_AtianiNetwork\spks_t20_M.mat')
plot(sum(spks),256:-1:1); hold on;
attendM_cnt = sum(spks);

load('data\\041_harmonics_AtianiNetwork\spks_t20_F.mat')
plot(sum(spks),256:-1:1); hold on;
attendF_cnt = sum(spks);

% load('data\\041_harmonics_AtianiNetwork\ICspks_t20.mat')
% plot(sum(spk_IC),256:-1:1); hold on;
% noattend_cnt = spk_IC;

load('data\\041_harmonics_AtianiNetwork\spks_t20_none.mat')
plot(sum(spks),256:-1:1); hold on;
noattend_cnt = sum(spks);

axis tight
xlabel('spike count')
set(gca,'box','off')
legend('attend F','attend M','no attend')

yticks(1:16:256)
yticklabels(round(cf(256:-16:1)/1000,2))

%% Suppresstion of information vs f0
figure;
plot(attendM_cnt./noattend_cnt,256:-1:1); hold on;
plot(attendF_cnt./noattend_cnt,256:-1:1); 
yticks(1:16:256)
yticklabels(round(cf(256:-16:1)/1000,2))

%% Passive state, network output spike rasters
addpath('..\BOSSA\peripheral')
addpath('..\BOSSA\plotting')
% load('data\\041_harmonics_AtianiNetwork\ICspks_t20.mat')
% low_freq = 50; %min freq of the filter
% high_freq = 8000;
% numChannel = 256;
fs = 40000;
% [nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
% fcoefs=MakeERBFilters(fs,cf,low_freq);
load('data\\041_harmonics_AtianiNetwork\spks_t20_none.mat')

figure('position',[200 200 1150 550]);
subplot('position',[0.1 0.1 0.7 0.8])
plotSpikeRasterFs(logical(spks'), 'PlotType','vertline', 'Fs',fs);
xlim([0 length(spks)/fs*1000])
title('Passive spike rasters')
xlabel('time (ms)')
ylabel('frequency (khz)')
axisinc = 1:256/16:256;
yticks(axisinc)
yticklabels(round(cf(axisinc)/1000,2))

%% change in spike founds at f0
