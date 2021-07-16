% evaluate IC_spikes
% assuming you've ran gen_IC_spikes
% need the variables
%   wavs
%   fcoefs
addpath('eval_scripts')
% addpath('C:\Users\Kenny\Desktop\GitHub\BOSSA\Plotting')
addpath(fullfile('..','BOSSA','Plotting'))

numTalkers = 2;
% calculate NCC
for j = 1:20

    ic(1).wav = wavs(j).m1;
    ic(2).wav = wavs(j).tgt;
% ic(1).wav = wavs(j).m1;
% ic(2).wav = wavs(j).tgt;
% ic(3).wav = wavs(j).m2;
ref = struct();
% refLabels = {'M1','T','M2'};
refLabels = {'M1','T'};

% correlation(reconstruction of M1, original M1)
% correlation(reconstruction of T, original T)
% correlation(reconstruction of M2, original M2)
windowLen = 0.250; %seconds;

figure;
for i = 1:numTalkers
    tic;
    ref(i).tf = ERBFilterBank(ic(i).wav(:,1),fcoefs)+ERBFilterBank(ic(i).wav(:,2),fcoefs);
    recon = data(i+(j-1)*numTalkers).recon;
    reconTF = ERBFilterBank(recon(:,1),fcoefs)+ERBFilterBank(recon(:,2),fcoefs);
    ev(i).cc = movxcorrKC(ref(i).tf,reconTF,fs*windowLen);
%     ev(i).scc = smoothdata(ev(i).cc,'interp',10000);
    subplot(2,numTalkers,i);
    plot_vspgram((ic(i).wav(:,1)+ic(i).wav(:,2))/mean(rms(ic(i).wav))*0.05,fs)
    title('Ref')
    caxis([-150,10])
    subplot(2,numTalkers,i+numTalkers);
    plot_vspgram(sum(recon,2)/mean(rms(recon))*0.05,fs);
    title('Recon')
    caxis([-150,10])
    toc
end

h1 = figure('position',[200 200 600 300]);
% h2 = figure('position',[200 200 600 300]);

for i = 1:numTalkers
    set(0, 'CurrentFigure', h1)
    stem(ev(i).cc,'LineStyle','none','marker','o'); hold on;
%     set(0, 'CurrentFigure', h2)
%     plot(ev(i).scc,'linewidth',2); hold on;
end

set(0, 'CurrentFigure', h1)
legend(refLabels)
ylabel('NCC')
title(['set ' num2str(j,'%02i')])
saveas(h1,['TM90 NCC results ' num2str(j,'%02i') '.tiff'])
end
% set(0, 'CurrentFigure', h2)
% legend(refLabels)
% ylabel('Smoothed NCC')
% axis tight


% addpath('C:\Users\Kenny\Desktop\GitHub\Plotting')
% load('stimuli\IC_spks_trioSet_01.mat','spk_IC')
% spks = spkTime2Train(spk_IC,fs,length(mixed));
% figure;
% for i = 1:5
%     subplot(1,5,i)
%     icSpikes = logical(squeeze(spks(:,:,i))');
%     plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'Fs',40000);
%     xlim([0 2000])
%     if i==1, ylabel('IC spikes'); end
%     set(gca,'Ytick',[1:64],'YtickLabel',cf)
% end

