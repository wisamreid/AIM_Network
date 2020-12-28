%% constant-octaves

bw = 1; %octaves = fwhm
sigma = bw/2.3548;
faxis = linspace(-2.3219,4,1000);
[minVal, closestIndex] = min(abs(2.^faxis'-x')); %map to ERB axis

% for all freq channels, do
clear h
i = 1;
for f0 = x'
    hGauss = exp(-0.5*((faxis-log2(f0))/sigma).^2);
    h(i,:) = hGauss(closestIndex);
    i = i+1;
end
figure;
uimagesc(flipud(x),flipud(x),flipud(fliplr(h)));
view([180,90])
xlabel('to, freq, ERB scale (kHz)')
ylabel('from, freq, ERB scale (kHz)')
title('Netcon for bw = 1 octave, interpolated')

figure;
imagesc(h);
xlabel('to, freq chan #')
ylabel('from, freq chan #')
title('Netcon for bw = 1 octave')
%%  IE netcon; constant BW (in terms of kHz)
clear h
paramH.IE_BW = 1;
sigma = paramH.IE_BW/2.3548;
x = cf/1000;
i = 1;
for f0 = 4
    hGauss = exp(-0.5*((x-f0)/sigma).^2);
    hcos = cos(2*pi*paramH.BTM*(x-f0)+paramH.phase*pi);
    h(i,:) = hGauss.*hcos;  
    i = i+1;
end
figure;plot(x,hGauss)

%%

IEnetcon = h - eye(nFreqs);
% figure; uimagesc(flipud(x),flipud(x),flipud(fliplr(IEnetcon)))
% xlabel('frequency (kHz)')
% ylabel('frequency (kHz)')

figure; imagesc(x,x,IEnetcon)
view([180,90])
% axisinc = 1:256/16:256;
% yticks(axisinc)
% yticklabels(round(cf(axisinc)/1000,2))
% xticks(axisinc)
% xticklabels(round(cf(axisinc)/1000,2))
% xtickangle(30)
xlabel('frequency (kHz)')
ylabel('frequency (kHz)')

%% EC netcon;  constant BW (in terms of kHz)
paramH = options.paramH;
paramH.EC_BW = 1.5
x = options.cf/1000; %unit in kHz
i = 1;
for f0 = x'
    hGauss = exp(-0.5*((x-f0)/paramH.EC_BW).^2);
    hcos = cos(2*pi*paramH.BTM*(x-f0)+paramH.phase*pi);
    h(i,:) = hGauss.*hcos;  
    i = i+1;
end
ECnetcon = h;

figure; imagesc(ECnetcon)
axisinc = 1:256/16:256;
yticks(axisinc)
yticklabels(round(cf(axisinc)/1000,2))
xticks(axisinc)
xticklabels(round(cf(axisinc)/1000,2))
xtickangle(30)
xlabel('frequency (kHz)')
ylabel('frequency (kHz)')
%% spacing (number of freq channels) for 1 octave
% cf in Hz
for i = 1:numel(cf)
[~,tgtTDidx1] = min(abs(cf-cf(i)));
[~,tgtTDidx2] = min(abs(cf-cf(i)*2));
d(i) = tgtTDidx1-tgtTDidx2;
end
figure; plot(cf,d)
xlabel('frequency (Hz)')
ylabel('# frequencies for 1 octave')

%% old way of doing things
paramH.IE_BW = 0.35;
hGauss = exp(-0.5*((x-paramH.t0)/paramH.IE_BW).^2);
hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
h = hGauss.*hcos;
h = h(h>0.0001); %remove 0 values
IEnetcon = conv2(eye(nFreqs),h);
[~,b] = max(h);
IEnetcon = IEnetcon(b:b+nFreqs-1,:);
IEnetcon = IEnetcon - eye(nFreqs);

figure;imagesc(IEnetcon)