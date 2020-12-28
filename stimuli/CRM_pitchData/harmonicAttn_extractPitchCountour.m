% extractPitchCountour
addpath('..\BOSSA\plotting')
load('C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\stimuli\CRM Stimuli TM 0deg.mat','crm_number','talkers')
wavLoc = 'C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\stimuli\CRM MF 0deg colocated';
keyLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\toolboxes\CRM\Keys';
sex = 'F';
if strcmp(sex,'F'), sexIdx = 2; else sexIdx = 1; end

set(0,'DefaultFigureVisible',0)
for speakerSet = 9
[wav,fs] = audioread([wavLoc filesep sprintf('%s_set_%02i.wav',sex,speakerSet)]);

%% voicebox calculation - run this to get pv
addpath(genpath('C:\Users\Kenny\Desktop\GitHub\sap-voicebox'))
[fx,tx,pv,fv]=v_fxpefac(wav(:,1),fs,0.01);

%% matlab calculation
rmpath(genpath('C:\Users\Kenny\Desktop\GitHub\sap-voicebox'))
[~,~,h] = plot_vspgram(wav(:,1),fs,[],1);
hold on;
ylim([0 800])
set(h,'position',[300 300 700 300]);

% use matlab's built-in pitch function to estimate f0
method = {'SRH','NCF','PEF','CEP','LHS'};
clear f0;
for i = 1:numel(method)
    f0(i,:) = pitch(wav(:,1),fs,'method',method{i},'windowlength',fs*0.0522,'overlaplength',fs*(0.042));
end
x = linspace(0,tx(end),length(f0));
plot(x,f0,'linewidth',2)
l = legend(method);

% plot median of the 5 pitch estimates
medf0 = median(f0);
plot(x,medf0,'linewidth',2)

% plot voiced frames
meanf0 = mean(medf0(pv>0.5));
stdf0 = std(medf0(pv>0.5));
voiced = (medf0 > (meanf0 - 2*stdf0)) & (medf0 < (meanf0 + 2*stdf0));
scatter(x(voiced),medf0(voiced),'g')

%retrieve talker timing keys
load([keyLoc filesep 'Talker' num2str(talkers(speakerSet),sexIdx) filesep strrep(crm_number{speakerSet,sexIdx},'bin','mat')])
line([startcall/fs startcall/fs],ylim);
line([endcall/fs endcall/fs],ylim);
line([startw/fs startw/fs],ylim);
line([endw/fs endw/fs],ylim);
l.String(6:end) = [];

%% save data
data.wav = wav;
data.f0 = medf0;
data.voiced = voiced;
data.taxis = x;
data.ID = [sex num2str(speakerSet)];
save(['stimuli\CRM_pitchData\' data.ID '.mat'],'data');
saveas(h,['stimuli\CRM_pitchData\' data.ID '.tiff'])
end
close all
set(0,'DefaultFigureVisible',1)