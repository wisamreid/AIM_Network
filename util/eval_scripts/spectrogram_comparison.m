addpath('C:\Users\kfcho\Dropbox\Sen Lab\phase_estimation_temp\voicebox');
[test,fs] = audioread('eval_data\IC_debug\TM_00_90_set_01_target.wav');
%%
%peripheral filtering parameters
low_freq = 100; %min freq of the filter
high_freq = 8000;
numChannel = 64;
 
%apply peripheral filters
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

%%
testFilt = ERBFilterBank(mixed,fcoefs);
 
figure;
temp = db(testFilt);
% temp = db(signal.sL);
temp(temp<-80) = -80;
imagesc(temp)
yticks(1:64)
yticklabels(cf)
set(gca,'yscale','linear')

%%
p.fstep=5;        % frequency resolution of initial spectrogram (Hz)
p.fmax=8000;      % maximum frequency of initial spectrogram (Hz)
fres = 20;        %bandwidth, from fxpefac
fmin = 0;
faxis_spec = [fmin p.fstep p.fmax];
%run this line to extract some parameters
[tx,F,B,~,Brawfilt,spec] = spgrambw(test,fs,fres,faxis_spec,[]); 
 
[sf,~]=enframe(test,spec.win,spec.ninc);
frameLen = length(spec.win);
Braw=rfft(sf,frameLen,2);
Braw = Braw';
Bmag = abs(Braw);
figure;
fxraw = [0:frameLen-1]/frameLen*fs;
fxraw = fxraw(1:frameLen/2);
imagesc(tx,fxraw,20*log(Bmag)); set(gca,'ydir','normal')
set(gca,'ydir','normal')
ylim([0,8000])
yticks(flipud(cf))
yticklabels(flipud(cf))
