% This file solves the full I, C, R model with IC input
addpath(genpath('C:\Users\kfcho\Documents\GitHub\dynasim'))
%% solver params
time_end = 1.9e3; % ms

solverType = 'euler';
% solverType = 'modified_euler';
% % solverType = 'rk4';
dt = 0.025; %ms % the IC input is currently at dt=0.025

study_dir = fullfile(pwd, 'run', mfilename);
addpath('dependencies')
addpath('mechs')
addpath(genpath('../dynasim'))
%% copy IC spikes to study dir
if exist(study_dir, 'dir')
  rmdir(study_dir, 's');
end
mkdir(fullfile(study_dir, 'solve'));
spkLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\CISPA2.0\Data\001 IC_spk 64Chan150-8000hz\TM at 0_90 -FreqGainNorm talker4\';
spkList = ls([spkLoc '*IC.mat']);
spkFile = spkList(1,:);
spkICCopy = fullfile(study_dir, 'solve', 'IC_spks.mat');
copyfile( fullfile(spkLoc, spkFile), spkICCopy);

%spk_IC should have dimensions [time x freq chan x spatial chan]
spk_IC = load(fullfile(spkICCopy), 'spk_IC','fcoefs','cf');
spk_IC = spk_IC.spk_IC;
spk_IC_fs = 40e3;
spk_IC_t = 1/spk_IC_fs:1/spk_IC_fs:size(spk_IC,1)/spk_IC_fs;

%% neurons
% model structure
s=[];
nCells = size(spk_IC,2);
noise = 0; % low noise

% neuron populations
s.populations(1).name = 'preLI';
s.populations(1).equations = 'chouLIF';
s.populations(1).size = nCells;
s.populations(1).parameters = {'noise',noise};

s.populations(2).name = 'postLI';
s.populations(2).equations = 'chouLIF';
s.populations(2).size = nCells;
s.populations(2).parameters = {'noise',noise};


%% connection mechanisms
% explicitly defined mechs
% synAlpha={
%   'gSYN = 1; ESYN = 0; tauD = 1; delay = 0'
%   'f(x) = (exp(-x/tauD)).*(x>0)'
%   'netcon=ones(N_pre,N_post)'
%   'synAlpha(X,t) = gSYN .* ( f(t - tspike_pre - delay) * netcon ).*(X - ESYN)'
%   '@isyn += synAlpha(V_post,t)'
%   };

synDoubleExp={
  'gSYN_E = 1; ESYN_E = 0; tauD_E = 3; tauR_E = 1; delay = 0'
  'gSYN_I = 1; ESYN_I = 0; tauD_I = 10; tauR_I = 1'
  'f(x) = (exp(-x/tauD_E) - exp(-x/tauR_E)).*(x>0)'
  'g(x) = (exp(-x/tauD_I) - exp(-x/tauR_I)).*(x>0)'
  'E_netcon=eye(N_pre,N_post)'
  'I_netcon=eye(N_pre,N_post)'
  'synDoubleExp_E(X,t) = gSYN_E .* f(t - tspike_pre - delay) * E_netcon .*(X - ESYN_E)'
  'synDoubleExp_I(X,t) = gSYN_I .* g(t - tspike_pre - delay) * I_netcon .*(X - ESYN_I)'
  '@isyn += synDoubleExp_E(V_post,t) + synDoubleExp_I(V_post,t)'
  };

s.mechanisms(1).name='synDoubleExp';
s.mechanisms(1).equations=synDoubleExp;

% s.mechanisms(2).name='synAlpha';
% s.mechanisms(2).equations=synAlpha;

%% Netcons
% build preLI->postLI netcon matrix
% netcons have the size [N_pre,N_post]
inh_width = 5;
LINetcon = conv2(eye(nCells),ones(1,inh_width),'same')-eye(nCells);
% LINetcon = conv2(eye(nCells),ones(1,inh_width),'same');

% LINetcon = conv2(eye(nCells),1-hamming(inh_width),'same');

% xx = (inh_width-1)/2+1;
% for i = 1:xx-1
%     LINetcon(i*2:xx+i-1,i) = LINetcon(i*2:xx+i-1,i)*2;
%     LINetcon(end-xx+2-i:end-i*2+1,end-i+1) = LINetcon(end-xx+2-i:end-i*2+1,end-i+1)*2;    
% end
figure;imagesc(LINetcon);
title('LI Netcon Matrix')
xlabel('to cells')
ylabel('from cells')
       
% LINetcon = zeros(nCells);
% LINetcon(30,:) = 1;
% LINetcon(30,30) = 0;

% E_netcon = zeros(nCells);
% E_netcon(36,36) = 1;
E_netcon = eye(nCells);
%% mechanisms

% initialization
s.connections(1).direction='preLI->preLI';
s.connections(1).mechanism_list='IC_LateralInhibition';
% s.connections(end).parameters={'g_postIC',0.11}; % 1000 hz spiking
s.connections(end).parameters={'g_postIC',0.02}; % 100 hz spiking

% pre-postLI connections
ncm = 1    ; %netcon multiplier
ESYN_I = -80;
s.connections(end+1).direction='preLI->postLI';
s.connections(end).mechanism_list='synDoubleExp';
% s.connections(end).parameters={'gSYN',.07, 'tauR',3, 'tauD',50, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',1, 'tauD',1000, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; % ds iGABAa
s.connections(end).parameters={'gSYN_E',.13,'tauR_E',0.4, 'tauD_E',3, 'E_netcon',E_netcon,'ESYN_E',0,...
                               'gSYN_I',0,'tauR_I',1, 'tauD_I',10, 'I_netcon',LINetcon,'ESYN_I',-80};
%reversal potential ESYN


%% vary
vary = {
    '(preLI->postLI)','gSYN_I',[0.01:0.002:0.02];
%   '(I->I,R->R)', 'SpatialChan', 3;
%   '(I->I,R->R)', 'freqInd', 1;
};
nVary = calcNumVary(vary);
parfor_flag = double(nVary > 1); % use parfor if multiple sims


%% simulate
tic
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'vary',vary, 'debug_flag',1, 'verbose_flag',0,...
  'parfor_flag',parfor_flag);
toc

%%
addpath('eval_scripts')
if numel(size(spk_IC)) == 3
    temp = spk_IC(:,:,3);
elseif numel(size(spk_IC)) == 2
    temp = spk_IC;
end
maskIC = calcSpkMask(temp,40000,'alpha',.002);
figure;
imagesc(maskIC); title('IC mask');
% temp = vocode(maskIC,cf,'tone',fs);
% soundsc(temp,fs);


%% 
figure;
mask(1).Pre = calcSpkMask(data(1).preLI_V_spikes,40000,'alpha',.002);
imagesc(mask(1).Pre); title('preLI mask');
figure;
for jj = 1:nVary
    mask(jj).Post = calcSpkMask(data(jj).postLI_V_spikes,40000,'alpha',.2);

    subplot(2,3,jj)
    imagesc(mask(jj).Post); title(sprintf('postLImask, gSYN_I:%0.3f',data(jj).preLI_postLI_gSYN_I));
    yticks(1:nCells)
%     yticklabels(cf)
%     ylim([36-13+0.5,36+0.5]);

    disp(['preLI spikes: ' num2str(sum(sum(data(jj).preLI_V_spikes)))]);
    disp(['postLI spikes: ' num2str(sum(sum(data(jj).postLI_V_spikes)))]);

end

