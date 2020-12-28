% evaluate the tone stimulli & peripheral model

% cosine ramp
N = 20;
r = linspace(-pi,0,N);
w = (cos(r)+1)/2;

% stimuli
fs = 44100;
numTrials = 64;
f0 = logspace(log10(500),log10(16000),numTrials); %hz
tmax = 150/1000; %milliseconds; time axis
t = 0:1/fs:tmax; %
tones = zeros(length(f0),length(t));
for i = 1:length(f0)
    tones(i,:) = cos(2*pi*f0(i)*t);
    tones(i,1:N) = tones(i,1:N).*w;
    tones(i,end-N+1:end) = tones(i,end-N+1:end).*fliplr(w);
end
tones = [tones zeros(size(tones))];

% ERB filterbank parameters
low_freq = 500; %min freq of the filter
high_freq = 16000;
numChannel = 64;
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);
% plot_filterBankSupport(fcoefs,nf);
for i = 46%1:numTrials% [10,46,58]%
    ILDs = gen_ILD(low_freq,high_freq,numChannel,fs,0);
    sig_filt = ERBFilterBank(tones(i,:));
    s_filt.sL = zeros(size(sig_filt));
    % apply ILD additively and zero-pad ITDs
    for j = 1:numChannel
        s_filt.sL(j,:) = sig_filt(j,:)+ILDs(j);
    end
    s_filt.sR = sig_filt;
    s_filt.sL = s_filt.sL;
    
    figure;
    subplot(1,3,1);
    imagesc(sig_filt)
%     ylim([0 64])
    title('ERB-filtered signal')
    subplot(1,3,2);
    plot(cf/1000,sum(sig_filt,2)); view(90,90)
    xlim([0 16])
    title({'ERB-filtered signal', 'collapsed across time'})
    subplot(1,3,3);
    plot(cf/1000,sum(abs(sig_filt),2)); view(90,90)
    xlim([0 16])
    title({'abs(ERB-filtered signal)', 'collapsed across time'})
%     subplot(1,4,4);
    
%     plot(sum(s_filt.sL,2)); view(90,90)
%     xlim([0 64])
%     title('ERB-filtered & added ILD/ITD')
%     xlabel('Frequency channel')
    pause(0.1);
    
    sigStrength(:,i) = sum(abs(sig_filt),2);
end
figure;
imagesc(sigStrength)
xlabel('stimuli frequency (kHz)')
xticks((1:8:64))
xticklabels(round(f0(1:8:64)/1000,2))
ylabel('frequency band cf')
yticks((1:8:64))
yticklabels(round(f0(64:-8:1)/1000,2))

figure; subplot(2,1,1); plot(tones(i,:));
subplot(2,1,2);
imagesc(db(sig_filt))