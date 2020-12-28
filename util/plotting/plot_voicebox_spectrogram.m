function plot_voicebox_spectrogram(input,fs,cf)
% utilizes voicebox's spectrogram calculation to plot the input
p.fstep=5;        % frequency resolution of initial spectrogram (Hz)
p.fmax=8000;      % maximum frequency of initial spectrogram (Hz)
fres = 20;        %bandwidth, from fxpefac
fmin = 0;
faxis_spec = [fmin p.fstep p.fmax];
%run this line to extract some parameters
[tx,F,B,~,Brawfilt,spec] = spgrambw(input,fs,fres,faxis_spec,[]); 
 
[sf,~]=enframe(input,spec.win,spec.ninc);
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