% smoothSTRFs
% dataFolder = '038_Aritani_3b_EE0_test_Attend0';
% dataFolder = '038_Aritani_3b_EE0_test_Attend1';
% dataFolder = '038_Aritani_3a_ECBW1.1_EE0_Attend0';
% dataFolder = '038_Aritani_3a_ECBW1.1_EE0_Attend1';

% dataFolder = '040_Fritz_Attend1';
dataFolder = '040_Fritz_Attend0';

load(['C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\data\' dataFolder '\STRF.mat']);

% plot marginals
figure;
plot(sum(runningFR,2),log2(f0/1000)); hold on;
ylabel('total spike count')
xlabel('f0 (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))


% plot STRF from data
figure;
yMinMax = [min(ySpacing) max(ySpacing)];
imagesc(x,log2(f0/1000),runningFR);
ylabel('pure tone frequency kHz')
xlabel('time (ms)')
title(['neuron with cf ' num2str(bestFreq,'%.02f')])
set(gca,'ydir','normal')
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
colormap jet
caxis([-160 160])

%% smooth STRF and plot
figure;
smoothed = imgaussfilt(runningFR,[2,10]);
% smoothed = medfilt2(runningFR,[5 5]);
imagesc(x,log2(f0/1000),smoothed);
ylabel('pure tone frequency kHz')
xlabel('time (ms)')
title(['neuron with cf ' num2str(bestFreq,'%.02f')])
set(gca,'ydir','normal')
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
colormap jet
caxis([-160 160])

%% SVD?
[U,S,V] = svd(runningFR);
nSV = 1;
reduced = U(:,1:nSV)*S(1:nSV,1:nSV)*V(:,1:nSV)';
smoothed = imgaussfilt(reduced,[2,40]);

figure;
imagesc(x,log2(f0/1000),reduced);
ylabel('pure tone frequency kHz')
xlabel('time (ms)')
title({'SVD',dataFolder,['neuron with cf ' num2str(bestFreq,'%.02f')]})
set(gca,'ydir','normal')
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
colormap jet
caxis([-160 160])

figure
imagesc(x,log2(f0/1000),smoothed);
ylabel('pure tone frequency kHz')
xlabel('time (ms)')
title({'SVD then smooth',['neuron with cf ' num2str(bestFreq,'%.02f')]})
set(gca,'ydir','normal')
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
colormap jet
caxis([-220 220])
