% recon_eval_LIN
% for evaluation the output of lateral inhibition networks

%% try reconstruction
dataloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\CISPA2.0\Data\001 IC_spk 64Chan150-8000hz\TM at 0_90 -FreqGainNorm talker4\';
[mixed,fs] = audioread([dataloc 'TM_00_90_set_01_mixed.wav']);
[target,fs] = audioread([dataloc 'TM_00_90_set_01_target.wav']);
target = target/max(abs(target));

low_freq = 150; %min freq of the filter
high_freq = 8000;
numChannel = 64;
 
%apply peripheral filters
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);
targetFilt = ERBFilterBank(target,fcoefs);
mixedFiltL = ERBFilterBank(mixed(:,1),fcoefs);
mixedFiltR = ERBFilterBank(mixed(:,2),fcoefs);
    
clear mixedEnvL mixedEnvR targetEnv targetConvEnv
for i = 1:nf
    mixedEnvL(i,:)= envelope(mixedFiltL(i,:)); %envlope of mixture
    mixedEnvR(i,:)= envelope(mixedFiltR(i,:)); %envlope of mixture
%     targetEnv(i,:) = envelope(targetFiltmono(i,:)); %envelope of target
end
%% IC masks
tempMask = mask(3).Post; %normalize mask
tempMask = tempMask/max(max(abs(tempMask)));
% STOI
[rstim1dual, rstim1mono] = applyMask(tempMask,mixedFiltL,mixedFiltR,1,'filt');
[rstim2dual, rstim2mono] = applyMask(tempMask,mixedEnvL,mixedEnvR,1,'env',cf);
rstim3 = vocode(tempMask,cf,'tone');
st1 = runStoi(double(rstim1mono),target,fs,fs);
st2 = runStoi(double(rstim2mono),target,fs,fs);
st3 = runStoi(double(rstim3),target,fs,fs);
out = [st1 st2 st3]

% audiowrite([dataloc 'BinarizedICmaskRecon_maskFiltered.wav'],rstim1dual/max(max(abs(rstim1dual))),fs);
% audiowrite([dataloc 'BinarizedICmaskRecon_maskFilteredVoc.wav'],rstim2dual/max(max(abs(rstim2dual))),fs);
% audiowrite([dataloc 'BinarizedICmaskRecon_vocode.wav'],rstim3/max(abs(rstim3)),fs);

%% PESQ
audiowrite('temp.wav',rstim2mono,fs)
pq1 = pesqmain('+16000','eval_data\IC_debug\TM_00_90_set_01_target.wav','temp.wav')
%% plot mask
temp1 = db(tempMask.*targetFilt(:,1:76000));
temp1(temp1<-100) = -100;
figure;
imagesc(temp1)
yticks(1:36)
yticklabels(cf)
%% plot target filt
figure;
temp = db(targetFilt);
temp(temp<-80) = -80;
imagesc(temp)
yticks(1:36)
yticklabels(cf)
% set(gca,'yscale','linear')

%% spectrograms of reconstruction
addpath('C:\Users\kfcho\Dropbox\Sen Lab\phase_estimation_temp\voicebox');
plot_voicebox_spectrogram(target,fs,cf); ylim([0,1200]); title('target')
plot_voicebox_spectrogram(rstim1mono,fs,cf); ylim([0,1200]); title('maskfiltered')
plot_voicebox_spectrogram(rstim2mono,fs,cf); ylim([0,1200]); title('vocodedmaskfiltered')
plot_voicebox_spectrogram(rstim3,fs,cf); ylim([0,1200]); title('vocoded')

%% without LIN
tempMask = calcSpkMask(spk_IC,40000,'alpha',.02);
tempMask = tempMask/max(max(abs(tempMask)));
% STOI
[rstim1dual, rstim1mono] = applyMask(tempMask,mixedFiltL,mixedFiltR,1,'filt');
[rstim2dual, rstim2mono] = applyMask(tempMask,mixedEnvL,mixedEnvR,1,'env',cf);
rstim3 = vocode(tempMask,cf,'tone');
st1 = runStoi(double(rstim1mono),target,fs,fs)
st2 = runStoi(double(rstim2mono),target,fs,fs)
st3 = runStoi(double(rstim3),target,fs,fs)
out = [st1 st2 st3]
plot_voicebox_spectrogram(rstim1mono,fs,cf); ylim([0,1200]); title('maskfiltered')
plot_voicebox_spectrogram(rstim2mono,fs,cf); ylim([0,1200]); title('vocodedmaskfiltered')
plot_voicebox_spectrogram(rstim3,fs,cf); ylim([0,1200]); title('vocoded')