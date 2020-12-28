colors = {'blue','red'};
h = figure('position',[300 50 1450 900]);

Fs = 9;
Ms = 9;
for F = Fs
    FID = ['F' num2str(F)];
    for M = Ms
        MID = ['M' num2str(M)];
        dataID = {FID,MID};
        for i = 1:numel(dataID)
            load(['stimuli\CRM_pitchData\' dataID{i} '.mat'],'data');
            for n = 1:20
                scatter(data.taxis(data.voiced),n*data.f0(data.voiced),colors{i}); hold on;
            end
        end
        title(dataID)
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        ylim([0 2000])
        saveas(h,['stimuli\CRM_pitchData\' dataID{1} dataID{2} '.tiff'])
        
        pause(0.1)
        clf
    end
end


%% plot harmonics with IC spikes
stimuliRoot = 'stimuli\CRM MF 0deg colocated';

trials = 1:20;
n = 15;
%peripheral filter parameters
low_freq = 50; %min freq of the filter
high_freq = 8000;
numChannel = 256;
fs = 40000;
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

for trial = trials
    % load IC spikes
    spks = load([stimuliRoot filesep sprintf('IC_spks_0DegMF_Set_%02i.mat',trial)]);
    spk_IC = spkTime2Train(spks.spk_IC,spks.fs,0.3*fs);
    spk_IC(1:400,:) = 0; %remove artifact
    
    % load pitch contours
    Fdata = load(sprintf('stimuli\\CRM_pitchData\\F%i.mat',trial),'data'); %attend female
    Mdata = load(sprintf('stimuli\\CRM_pitchData\\M%i.mat',trial),'data'); % attend male
    
    % plot
    figure('position',[50 300 900 400]);
    plotSpikeRasterFs(logical(spk_IC'), 'PlotType','vertline', 'Fs',fs);
    xlim([0 length(spk_IC)/fs*1000])

    hold on;
    for target = {'M','F'}
        if strcmp(target,'M') 
            pitchData = Mdata.data;
            t_axis = Mdata.data.taxis(Mdata.data.voiced);
            color = [0.3010 0.7450 0.9330];
        elseif strcmp(target,'F')
            pitchData = Fdata.data;
            t_axis = Fdata.data.taxis(Fdata.data.voiced);
            color = [239, 154, 154]/255;
        end
        tgtf0 = pitchData.f0(pitchData.voiced);

        harmonics = (1:n)'*tgtf0;
        harmonics(harmonics>=high_freq) = NaN;
        harmonics(harmonics<low_freq) = NaN;

        % find indicies corresponding to harmonics
        tgtTDidx = zeros(size(harmonics));
        for i = 1:numel(harmonics)
            if isnan(harmonics(i))
                tgtTDidx(i) = NaN;
            else
                [d(i),tgtTDidx(i)] = min(abs(cf-harmonics(i)));
            end
        end
        for i = 1:n
            s = scatter(t_axis*1000,1+tgtTDidx(i,:),[],'MarkerEdgeColor',color,'Marker','.');
        end
    end
    xlabel('time (ms)')
    ylabel('frequency (khz)')
    axisinc = 1:256/16:256;
    yticks(axisinc)
    yticklabels(round(cf(axisinc)/1000,2))
    saveas(gcf,['stimuli\CRM_pitchData\harmonics trial' num2str(trial,'%02i') '.tiff'])
end

