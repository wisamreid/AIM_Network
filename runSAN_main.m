% runSAN_sims

% in a consecutive talker scenario: concatenate .wav files and spike files,
% then run the network

% note: it is not recommended to set study_dir on a different network (i.e.
% research drive) from the local machine. Saving run data will take
% forever.

% @Kenny Chou, Boston Univ
% 2019 September - created
% 2020-02-05 - combined Consecutive/simultaneous stimuli cases

% to do:
% - move stimuli prep to a new file

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('PISPA2.0',pathCell)), rmpath(genpath('../PISPA2.0')); end
addpath('dependencies')
addpath('mechs')
addpath('..\PISPA2.0\plotting')
addpath(genpath('..\dynasim'))

% ================ parameters ==============
expName = '005 modulate i-itonic singleTalker';
% expVar = 'st_st gpreI2';
expVar = 'i-itonic';
expRootLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network';
dataFolder = [expRootLoc filesep 'data' filesep expName]; %things are saved to here
% scenario = 'consecutive';
% scenario = 'simultaneous';
scenario = 'singular';

talkerSet = 20;

% network parameters
varies = struct();
varies(1).conxn = 'I';
varies(end).param = 'Itonic';
varies(end).range = 0.18:0.005:0.22;

% varies(end+1).conxn = 'st->st';
% varies(end).param = 'g_preI2';
% varies(end).range = 0.009:-0.001:0.001;

% name of raster plot files based on parameters
% name = expVar x vary.range
a = cellstr(num2str(varies(1).range','%.3f'));
b = {expVar};
[Bx,Ax] = ndgrid(1:numel(b),1:numel(a));
rasterFileNames = strcat(b(Bx(:)),'_',a(Ax(:)))
% =================== end ==================


% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run' expName currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% load data, structure to correct format and save to study dir
switch scenario
    case 'consecutive'
        fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\broad mode - staggered talkers 64Chan500-20000hz\CRM talker4';
        chans = {'00','90','-90'};
        refLabels = {'S2','S3','S1'};
        swav123 = [];
        spk123 = [];
        for i = 1:3
            ic(i).wav = audioread([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{i} 'degAz_target_conv.wav']);
            ic(i).data = load([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{i} 'degAz_SpkIC.mat']);
            ic(i).spksbin = spkTime2Train(ic(i).data.spk_IC,ic(i).data.fs,length(ic(i).wav));
            swav123 = [swav123;ic(i).wav];
            spk123 = [spk123;ic(i).spksbin];
        end
        spk_IC = spk123;
        target = swav123;
        mixed = target;
        save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
        fs = ic(1).data.fs;
        fcoefs = ic(1).data.fcoefs;
        cf  = ic(1).data.cf;
    case 'singular'
        fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\broad mode - staggered talkers 64Chan500-20000hz\CRM talker4';
        chans = {'00'};
        tgtwav = audioread([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{1} 'degAz_target_conv.wav']);
        spkData = load([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{1} 'degAz_SpkIC.mat']);
        spk_IC = spkTime2Train(spkData.spk_IC,spkData.fs,length(tgtwav));
        save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
        fs = spkData.fs;
        fcoefs = spkData.fcoefs;
        cf  = spkData.cf;
    case 'simultaneous'
        fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\simultaneous talkers 64Chan500-20000hz\CRM talker4';
        mixedName = ls([fileloc filesep sprintf('*set_%02i*mixed.wav',talkerSet)]);
        mixed = audioread([fileloc filesep mixedName]);
        spkName = ls([fileloc filesep sprintf('*set_%02i*SpkIC.mat',talkerSet)]);
        spks = load([fileloc filesep spkName]);
        if ~isfield(spks,'fs'), spks.fs = 40000; end
        spk_IC = spkTime2Train(spks.spk_IC,spks.fs,length(mixed));
        refNames = {'target_conv','masker1_conv','masker2_conv',};
        refLabels = {'S2','S3','S1'};
        for i = 1:3
            temp = ls([fileloc filesep sprintf('*set_%02i*%s.wav',talkerSet,refNames{i})]);
            ref(i).wav = audioread([fileloc filesep temp]);
        end
        tgt = ref(1).wav;

        save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
        simLen = length(spk_IC);
        fs = spks.fs;
        fcoefs = spks.fcoefs;
        cf  = spks.cf;
end

% run network
options.expName = expName;
options.rasterFileNames = rasterFileNames;
options.saveLoc = dataFolder;
options.plotRasters = 1;
options.fs = fs;

options.rcorr = 1;
options.fcoefs = fcoefs;
options.cf = cf;
data = runSAN(study_dir,spk_IC,varies,options);
save([study_dir filesep 'options.mat'],'options');
% %% Evaluate Output Intelligibilty
% addpath('eval_scripts')
% % targetName = ls(sprintf('%s*set_%02i_*target_conv.wav',spkLoc,talkerSet))
% % mixedName = ls(sprintf('%s*set_%02i_*target_conv.wav',spkLoc,talkerSet))
% % target = audioread([spkLoc targetName]);
% % mixed = audioread([spkLoc mixedName]);
% cSpikes = ([data.C_V_spikes])';
% params = struct();
% params.fcoefs = options.fcoefs;
% params.cf = options.cf;
% params.fs = options.fs;
% params.spatialChan = 3;
% params.delay = 0;
% 
% % FRmask reconstruction
% params.type = 1;
% params.maskRatio = 0;
% params.tau = 0.018;
% [st,rstim] = recon_eval(cSpikes',target,mixed,params);
% 
% %% calculate NCC
% ic(2).wav = [zeros(length(ic(1).wav),2); ic(2).wav];
% ic(3).wav = [zeros(length(ic(2).wav),2); ic(3).wav];
% ref = struct();
% 
% for i = 1:3    
%     ref(i).tf = ERBFilterBank(ic(i).wav(:,1),ic(i).data.fcoefs)+ERBFilterBank(ic(i).wav(:,2),ic(i).data.fcoefs);
%     ev(i).cc = movxcorrKC(ref(i).tf,rstim.tf.L+rstim.tf.R,fs*0.02);
%     ev(i).scc = smoothdata(ev(i).cc,'movmean',10000);
%     toc
% end
% 
% % plot
% h1 = figure('position',[200 200 600 300]);
% h2 = figure('position',[200 200 600 300]);
% 
% for i = 1:3
%     set(0, 'CurrentFigure', h1)
%     plot(ev(i).cc); hold on;
%     set(0, 'CurrentFigure', h2)
%     plot(ev(i).scc,'linewidth',2); hold on;
% end
% 
% set(0, 'CurrentFigure', h1)
% legend(refLabels)
% ylabel('NCC')
% 
% set(0, 'CurrentFigure', h2)
% legend(refLabels)
% ylabel('Smoothed NCC')
% axis tight
% 
% figure('position',[200 200 600 300]);
% plot_db(rstim.tf.L+rstim.tf.R,120);