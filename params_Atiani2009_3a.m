% stimuli path, data folder names
stimuliRoot = 'stimuli\PureTones\256chan 0.2-16khz rms15\';
expPrefix = '036';
expName = 'Aritani_3a_ECBW1.1';

%peripheral filter parameters
low_freq = 200; %min freq of the filter
high_freq = 16000;
numChannel = 256; %dictates the number of freq channels in the model
fs = 40000;
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

% find the subset of stimuli trials to run
f0 = logspace(log10(200),log10(16000),256); % original stimuli vector!
minf0 = 2500;
maxf0 = 8000;
[~, minf0idx]=min(abs(f0-minf0)); 
[~, maxf0idx]=min(abs(f0-maxf0));
trials = minf0idx:maxf0idx;
subf0 = f0(minf0idx:maxf0idx); % subset of stimuli inputs

% attention? target channel? best channel?
attention = false; % <------------------------------------------------
targetFreq = 4.8; % kHz <----------------------------------------------
bestFreq = 4.2; % kHz <----------------------------------------------
[~, tgtTDidx]=min(abs(cf-targetFreq*1000));
[~, BFchanidx]=min(abs(cf-bestFreq*1000)); 
tgtECgsyn = 1;

% network parameters
varies = struct();
varies(1).conxn = 'IC->IC';
varies(end).param = 'trial';
varies(end).range = trials;

varies(end+1).conxn = 'IC->IC';
varies(end).param = 'g_postIC';
varies(end).range = 10;

varies(end+1).conxn = 'IC->E';
varies(end).param = 'gSYN';
varies(end).range = 4;

varies(end+1).conxn = 'E->XI';
varies(end).param = 'gSYN';
varies(end).range = 3;

varies(end+1).conxn = 'E->C';
varies(end).param = 'gSYN';
varies(end).range = 1.25;

varies(end+1).conxn = 'E->E';
varies(end).param = 'gSYN';
varies(end).range = 3;

varies(end+1).conxn = 'XI->E';
varies(end).param = 'gSYN';
varies(end).range = 4;

varies(end+1).conxn = 'TD';
varies(end).param = 'Itonic';
varies(end).range = 8;

varies(end+1).conxn = 'TD->XI';
varies(end).param = 'gSYN';
varies(end).range = 3;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN';
varies(end).range = 4;

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic';
varies(end).range = 3;

varies(end+1).conxn = 'E';
varies(end).param = 'V_init';
varies(end).range = -80;

% varies(end+1).conxn = 'E';
% varies(end).param = 'noise';
% varies(end).range = 7;

% frequency convergence onto C neurons, for EC netcon
paramH.EC_BW= 1.1521; %bandwidth (?) units
% x-frequency inhibition from I to E neurons
paramH.IE_BW= 1.7; %bandwidth
paramH.BTM= 0.01 ; %BTM, modulation
paramH.t0= 2.9; % f0, peak f
paramH.phase= 0; % phase

ECnetcon = netcon_spread(cf/1000,paramH.EC_BW,'q'); % constant q spread
h = netcon_spread(cf/1000,paramH.IE_BW,'bw',-0.3);
IEnetcon = 1 - h - 0.5;
IEnetcon(IEnetcon < 0) = 0;

figure;
title('comparison of EC and IE spread')
plot(cf/1000,ECnetcon(:,BFchanidx)); hold on;
plot(cf/1000,IEnetcon(:,tgtTDidx));
xlabel('frequency (kHz)')
legend('EC spread','IE spread')

% params for plotting STRF
ySpacing = [1:0.25:4]; %powers of 2s
