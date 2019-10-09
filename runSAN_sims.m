% runSAN_sims

% set up paths
pathCell = regexp(path,pathsep,'split');
if any(strcmpi('PISPA2.0',pathCell)), rmpath(genpath('../PISPA2.0')); end
addpath('dependencies')
addpath('mechs')
addpath('..\PISPA2.0\plotting')
addpath(genpath('../dynasim'))

% ================ parameters ==============
scenario = 'consecutive';
% scenario = 'simultaneous';

fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\broad mode - staggered talkers 64Chan200-20000hz\CRM talker4';
talkerSet = 20;

% network parameters
varies = struct();
varies(1).conxn = 'I';
varies(end).param = 'Itonic';
varies(end).range = 0.36;

% =================== end ==================
% set up directory for simulation data
currentTime = char(datetime('now','Format','yyyyMMdd''-''HHmmss'));
study_dir = fullfile(pwd, 'run', [mfilename currentTime]);
mkdir(fullfile(study_dir, 'solve'));

% load data, structure to correct format and save to study dir
switch scenario
    case 'consecutive'
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
end

% run network
options.plotRasters = 1;
options.fs = fs;

options.rcorr = 1;
options.fcoefs = fcoefs;
options.cf = cf;
data = runSAN_ConsecutiveTalkers(study_dir,spk_IC,varies,options);

%% Evaluate Output Intelligibilty
addpath('eval_scripts')
% targetName = ls(sprintf('%s*set_%02i_*target_conv.wav',spkLoc,talkerSet))
% mixedName = ls(sprintf('%s*set_%02i_*target_conv.wav',spkLoc,talkerSet))
% target = audioread([spkLoc targetName]);
% mixed = audioread([spkLoc mixedName]);
cSpikes = ([data.C_V_spikes])';
params = struct();
params.fcoefs = options.fcoefs;
params.cf = options.cf;
params.fs = options.fs;
params.spatialChan = 3;
params.delay = 0;

target = options.target;
mixed = options.mixed;
% FRmask reconstruction
params.type = 1;
params.maskRatio = 0;
params.tau = 0.018;
[st,rstim] = recon_eval(cSpikes',target,mixed,params);

%% calculate NCC
ic(2).wav = [zeros(length(ic(1).wav),2); ic(2).wav];
ic(3).wav = [zeros(length(ic(2).wav),2); ic(3).wav];
ref = struct();

for i = 1:3    
    ref(i).tf = ERBFilterBank(ic(i).wav(:,1),ic(i).data.fcoefs)+ERBFilterBank(ic(i).wav(:,2),ic(i).data.fcoefs);
    ev(i).cc = movxcorrKC(ref(i).tf,rstim.tf.L+rstim.tf.R,fs*0.02);
    ev(i).scc = smoothdata(ev(i).cc,'movmean',10000);
    toc
end

% plot
h1 = figure('position',[200 200 600 300]);
h2 = figure('position',[200 200 600 300]);

for i = 1:3
    set(0, 'CurrentFigure', h1)
    plot(ev(i).cc); hold on;
    set(0, 'CurrentFigure', h2)
    plot(ev(i).scc,'linewidth',2); hold on;
end

set(0, 'CurrentFigure', h1)
legend(refLabels)
ylabel('NCC')

set(0, 'CurrentFigure', h2)
legend(refLabels)
ylabel('Smoothed NCC')
axis tight

figure('position',[200 200 600 300]);
plot_db(rstim.tf.L+rstim.tf.R,120);