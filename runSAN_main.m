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

% network structure:
nFreq = 64;
nLocs = 5;
% =================== end ==================

% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', ['run' expName currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% prep IC inputs
fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\broad mode - staggered talkers 64Chan500-20000hz\CRM talker4';
[fs, fcoefs, cf] = prepICinput(scenario,fileloc);

% prep I2 input
addpath('mechs');
genI2input(nFreq,nLocs,simLen,fs,study_dir)

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