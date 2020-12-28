% This scripts simulates the functional uses of AIM network in the spectral
% domain; A listner may want to attend the specific harmonics of a target
% talker. We use a male/female talker mixture as an example. We'll "attend"
% to the f0 of either the male or female talker. The f0 of the talkers were 
% estimated separately.
%
% Make changes to the parameter section before running the simuation.
%
% @Kenny F. Chou, Boston Univ, 2020

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('BOSSA',pathCell)), rmpath(genpath('../BOSSA')); end
addpath('mechs')
addpath('network_params')
addpath('util')
addpath('util\plotting')
addpath('util\eval_scripts')
addpath(genpath('..\dynasim'))

% =============== folders and paths ==============
% for logging results
stimuliRoot = 'stimuli\CRM MF 0deg colocated';
expPrefix = '049';
expName = sprintf('%s_harmonics_AtianiNetwork',expPrefix);
expRootLoc = pwd;
dataFolder = [expRootLoc filesep 'data' filesep expName]; %things are saved to here

% ================ parameters ==============
attention = true; % true or false
attendTarget = 'F'; % 'M' or 'F'
% trials = 2; %select M/F mixture trials

% peripheral filter parameters
low_freq = 50; %min freq of the filter
high_freq = 8000;
numChannel = 256;
fs = 40000;
[nf,cf,bw] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
fcoefs=MakeERBFilters(fs,cf,low_freq);

for trials = 1:19
    for attendTarget = {'M','F'}
% ===== f0 information =====
pitch.winLen = fs*0.0522;
pitch.noverlap = fs*(0.042);
Fdata = load(sprintf('stimuli\\CRM_pitchData\\F%i.mat',trials),'data'); %attend female
Mdata = load(sprintf('stimuli\\CRM_pitchData\\M%i.mat',trials),'data'); % attend male
IE_BWs = ones(1,20)*0.17;
IE_BWs([3,4,5,7,8,10,15,19]) = 0.25;
IE_BWs(3) = 0.3;
IE_BWs(4) = 0.3;
IE_BWs(7) = 0.2;

if strcmp(attendTarget,'M') 
    pitchData = Mdata.data;
    paramH.IE_BW= 0.07; %bandwidth, octaves
elseif strcmp(attendTarget,'F')
    pitchData = Fdata.data;
    paramH.IE_BW= IE_BWs(trials);
end
% tgtf0 = pitchData.f0(pitchData.voiced);
tgtf0 = pitchData.f0; 
unVoicedSeg = find(~pitchData.voiced);
for i = 1:length(unVoicedSeg)
    if unVoicedSeg(i) == 1, continue; end
    tgtf0(unVoicedSeg(i)) = tgtf0(unVoicedSeg(i)-1);
end
tgtf0(~pitchData.voiced);
% define the integer numbner of harmonics to attend
n = 1;

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
varies(end).range = 1.5;

varies(end+1).conxn = 'XI->E';
varies(end).param = 'gSYN';
varies(end).range = 3.5;

varies(end+1).conxn = 'TD';
varies(end).param = 'Itonic';
varies(end).range = 0;

varies(end+1).conxn = 'TD->XI';
varies(end).param = 'gSYN';
varies(end).range = 3;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN';
varies(end).range = 3; 

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic';
varies(end).range = 3;

varies(end+1).conxn = 'E';
varies(end).param = 'V_init';
varies(end).range = -80;

varied_param = find(cellfun(@length,{varies.range})>1);
if isempty(varied_param), varied_param = 1; end
expVar = [varies(varied_param).conxn '_' varies(varied_param).param];
expVar = strrep(expVar,'->','_');

% network structure:
nFreqs = numChannel;
nLocs = 1;

% 
% paramH.IE_BW = 0.1;

% frequency convergence onto C neurons, for EC netcon
paramH.EC_BW= 0.05; %bandwidth
paramH.BTM= 0.01 ; %BTM
paramH.t0= 2.9; % f0, peak f
paramH.phase= 0; % phase

ECnetcon = netcon_spread(cf/1000,paramH.EC_BW,'q'); % constant q spread
paramH.IE_BW = .01;
h = netcon_spread(cf/1000,paramH.IE_BW,'bw');
% h2 = netcon_spread(cf/1000,paramH.IE_BW,'bw');

% IEnetcon = h - eye(nFreqs);
IEnetcon = 1 - h - 0.5;
IEnetcon(IEnetcon < 0) = 0;
IEnetcon = IEnetcon * 0.35;

figure;
title('comparison of EC and IE spread')
plot(cf/1000,ECnetcon(:,124)); hold on;
plot(cf/1000,IEnetcon(:,124));
xlabel('frequency (kHz)')
legend('EC spread','IE spread')
xlim([1 2])

% name of raster plot files based on parameters
% name = expVar x vary.range
a = cellstr(num2str(varies(varied_param).range','%.04f'));
b = {expVar};
[Bx,Ax] = ndgrid(1:numel(b),1:numel(a));
rasterFileNames = strcat(b(Bx(:)),'_',a(Ax(:)))

% =================== end params ==================
    
% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run' expName currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% prep IC inputs
% ------------ load spk_IC here -----------------
for trial = trials
    spks = load([stimuliRoot filesep sprintf('IC_spks_0DegMF_Set_%02i.mat',trial)]);
    spk_IC = spkTime2Train(spks.spk_IC,spks.fs,0.3*fs);
    spk_IC(1:400,:) = 0; %remove artifact
%     spk_IC(3000:end,:) = []; %shorten stimulus
%     spk_IC = zeros(fs*0.160,64);
%     spk_IC(1000:fs*0.0018:end,tgtTDchan-2:tgtTDchan+1) = 1;
%     spk_IC(1000:fs*0.005:end,tgtTDchan+2:tgtTDchan+5) = 1;
%     sum(spk_IC)
    simLen = length(spk_IC);
    save(fullfile(study_dir, 'solve', sprintf('IC_spks_t%02i.mat',trial)),'spk_IC')
    save([dataFolder filesep sprintf('ICspks_t%02i.mat',trials)],'spk_IC');
end

% prep I2 input
addpath('mechs');
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

% time vector for pitch
t_idx = pitchData.taxis*fs;
t_inc = round(t_idx(2)-t_idx(1));
I2input = zeros(nf,length(spk_IC));
for i = 1:size(tgtTDidx,1) %frequency axis
    for j = 1:size(tgtTDidx,2) %time axis
        if ~isnan(tgtTDidx(i,j))
            I2input(tgtTDidx(i,j),1+t_inc*(j-1):t_inc*(j)) = 1;
        end
    end
end
I2input = 1-I2input'; % time in the 1st dimension
save([study_dir filesep 'solve' filesep 'i2input.mat'],'I2input')
% ====================== end input prep ==================================

% run network
options.expName = expName;
options.rasterFileNames = rasterFileNames;
options.saveLoc = dataFolder;
options.plotRasters = 1;
options.vizNetwork = 0;
options.nLocs = nLocs;
options.nFreqs = nFreqs;
options.fs = fs;
options.ECnetcon = ECnetcon;
options.IEnetcon = IEnetcon;

options.tgtTDchan = tgtTDidx;
options.attend = attention;

options.rcorr = 1;
options.fcoefs = fcoefs;
options.cf = cf;

tic;
[data,pops] = runHarmonicAttenNetwork(study_dir,spk_IC,varies,options);
toc;

% save parameters & model
model = data(1).model;
save([dataFolder filesep 'options.mat'],'options','model');
spks = data.C_V_spikes;
save([dataFolder filesep sprintf('spks_t%02i_%s.mat',trials,attendTarget{1})],'spks','cf');
    end
end

return;

%% analysis
n = 2; % number of harmonics to draw
% colorM = [0.3010 0.7450 0.9330];
colorM = [84, 19, 205]/255;
% colorF = [239, 154, 154]/255+25/255;
colorF = [205, 19, 69]/255;
for trial = 1:20
    spkCntTracker = [];
    load([dataFolder filesep sprintf('ICspks_t%02i.mat',trial)],'spk_IC');
    Fdata = load(sprintf('stimuli\\CRM_pitchData\\F%i.mat',trial),'data'); %attend female
    Mdata = load(sprintf('stimuli\\CRM_pitchData\\M%i.mat',trial),'data'); %attend male
    
    for AttendTarget = {'M','F'}
    %rasters
    load([dataFolder filesep sprintf('spks_t%02i_%s.mat',trial,AttendTarget{1})],'spks');
    figure('position',[200 200 1150 550]);
    subplot('position',[0.1 0.1 0.7 0.8])
    plotSpikeRasterFs(logical(spks'), 'PlotType','vertline', 'Fs',fs);
    xlim([0 length(spk_IC)/fs*1000])
    title(sprintf('Attend %s',AttendTarget{1}))
    xlabel('time (ms)')
    ylabel('frequency (khz)')
    axisinc = 1:256/16:256;
    yticks(axisinc)
    yticklabels(round(cf(axisinc)/1000,2))

    % marginals
    subplot('position',[0.8 0.1 0.1 0.8])
    plot(sum(spks),256:-1:1)
    axis tight
    xlabel('spike count')
    yticks([])
    set(gca,'box','off')

    % normalized spkCounts
    normSpkCount = sum(spks)./sum(spk_IC);

    hold on;
    % overlay harmonic traces for each sex
    for sex = {'M','F'}
        if strcmp(sex,'M') 
            pitchData = Mdata.data;
            t_axis = Mdata.data.taxis(Mdata.data.voiced);
            color = colorM;    
        elseif strcmp(sex,'F')
            pitchData = Fdata.data;
            t_axis = Fdata.data.taxis(Fdata.data.voiced);
            color = colorF;
        end
        tgtf0 = pitchData.f0(pitchData.voiced);
        meanf0 = mean(tgtf0);

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
        subplot('position',[0.1 0.1 0.7 0.8])
        hold on;
        for i = 1:n
            s = scatter(t_axis*1000,tgtTDidx(i,:),[],'MarkerEdgeColor',color,'Marker','.');
        end

        % marginals
        subplot('position',[0.8 0.1 0.1 0.8])
        line([100 250],256-ones(n,1)*median(tgtTDidx,2)','color',color)

        % find max normalized spk count at M and F f0
        [~, idx1] = min(abs(cf-meanf0-25));
        [~, idx2] = min(abs(cf-meanf0+25));
        spkCntTracker(end+1) = max(normSpkCount(idx1:idx2));
    end
    saveas(gcf,[dataFolder filesep 'network_output_trial' num2str(trials) 'sex' AttendTarget{1} '.tiff'])
    end
    maxSpkCount(trial).attenM_M = spkCntTracker(1);
    maxSpkCount(trial).attenM_F = spkCntTracker(2);
    maxSpkCount(trial).attenF_M = spkCntTracker(3);
    maxSpkCount(trial).attenF_F = spkCntTracker(4);

end
save([dataFolder filesep 'normSpkCounts.mat'], 'maxSpkCount');

figure;
plot(ones(20,1),[maxSpkCount.attenM_M],'o','color',[164, 202, 227]/255); hold on;
plot(ones(20,1),[maxSpkCount.attenM_F],'o','color',[231, 178, 155]/255);
plot(ones(20,1)*2,[maxSpkCount.attenF_M],'o','color',[164, 202, 227]/255);
plot(ones(20,1)*2,[maxSpkCount.attenF_F],'o','color',[231, 178, 155]/255);
plot([mean([maxSpkCount.attenM_M]) mean([maxSpkCount.attenF_M])],'o-','color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410]);
plot([mean([maxSpkCount.attenM_F]) mean([maxSpkCount.attenF_F])],'o-','color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980]);

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'M','F'})
xlabel('Attend Target')
ylabel('normalized spike count, %')
legend('Male f0','Female f0','location','best')
return;


%% distributions across frequency
figure;
for sex = {'M','F'}
load([dataFolder filesep sprintf('spks_t%02i_%s.mat',trials,sex{1})],'spks');
plot(sum(spks),256:-1:1); hold on;
axis tight
xlabel('spike count')
set(gca,'box','off')
end
plot(sum(spk_IC),256:-1:1); hold on;
axis tight
xlabel('spike count')
set(gca,'box','off')

