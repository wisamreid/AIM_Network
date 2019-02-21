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
study_IC_name = 'data\*SpkIC.mat';
spkLoc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\CISPA2.0\Data\001 64Chan100-8000\TM at 0_90 -FreqGainNorm talker4\';
copyfile( fullfile(spkLoc, study_IC_name), fullfile(study_dir, 'solve', 'IC_spks.mat') );

spk_IC = load(fullfile('eval_data', study_IC_name), 'spk_IC');
% spk_IC = permute(spk_IC.spk_IC,[1 3 2]);
spk_IC = spk_IC.spk_IC;
spk_IC_fs = 40e3;
spk_IC_t = 1/spk_IC_fs:1/spk_IC_fs:size(spk_IC,1)/spk_IC_fs;

%% neurons
% model structure
s=[];
nCells = 36;
noise = 0.01; % low noise

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
  'gSYN_E = 1; ESYN_E = 0; tauD_E = 4; tauR_E = 1; delay = 0'
  'gSYN_I = 1; ESYN_I = 0; tauD_I = 4; tauR_I = 1'
  'f(x) = (exp(-x/tauD_E) - exp(-x/tauR_E)).*(x>0)'
  'g(x) = (exp(-x/tauD_I) - exp(-x/tauR_I)).*(x>0)'
  'E_netcon=ones(N_pre,N_post)'
  'I_netcon=ones(N_pre,N_post)'
  'synDoubleExp_E(X,t) = gSYN_E .* ( f(t - tspike_pre - delay) * E_netcon ).*(X - ESYN_E)'
  'synDoubleExp_I(X,t) = gSYN_I .* ( g(t - tspike_pre - delay) * I_netcon ).*(X - ESYN_I)'
  '@isyn += synDoubleExp_E(V_post,t) + synDoubleExp_I(V_post,t)'
  };

s.mechanisms(1).name='synDoubleExp';
s.mechanisms(1).equations=synDoubleExp;

% s.mechanisms(2).name='synAlpha';
% s.mechanisms(2).equations=synAlpha;

%% Netcons
% build preLI->postLI netcon matrix
% netcons have the size [N_pre,N_post]
inh_width = 17;
inh_kernel = hamming(inh_width);
LINetcon = conv2(eye(nCells),inh_kernel,'same')-eye(nCells);
figure;imagesc(LINetcon)
       
% LINetcon = zeros(nCells);
% LINetcon(30,:) = 1;
% LINetcon(30,30) = 0;

E_netcon = eye(nCells);
%% mechanisms

% initialization
s.connections(1).direction='preLI->preLI';
s.connections(1).mechanism_list='IC_LateralInhibition';
% s.connections(end).parameters={'g_postIC',0.11}; % 1000 hz spiking
s.connections(end).parameters={'g_postIC',0.01}; % 100 hz spiking

% pre-postLI connections
prepostNetcon_E = eye(nCells);
prepostNetcon_I = zeros(nCells); 
s.connections(end+1).direction='preLI->postLI';
s.connections(end).mechanism_list='synDoubleExp';
% s.connections(end).parameters={'gSYN',.07, 'tauR',3, 'tauD',50, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',1, 'tauD',1000, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; % ds iGABAa
s.connections(end).parameters={'gSYN_E',.05,'tauR_E',0.4, 'tauD_E',10, 'E_netcon',prepostNetcon_E,'ESYN_E',0,...
                               'gSYN_I',.05,'tauR_I',0.4, 'tauD_I',10, 'I_netcon',prepostNetcon_I,'ESYN_I',-80};

                          
%reversal potential ESYN
ncm = 1; %netcon multiplier
ESYN_I = -80;
s.connections(end+1).direction='postLI->postLI';
s.connections(end).mechanism_list='synDoubleExp';
% s.connections(end).parameters={'gSYN',.07, 'tauR',3, 'tauD',50, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',1, 'tauD',1000, 'netcon',irNetcon, 'ESYN',-80};
% s.connections(end).parameters={'gSYN',.07, 'tauR',0.4, 'tauD',10, 'netcon',irNetcon, 'ESYN',-80}; % ds iGABAa
s.connections(end).parameters={'gSYN_E',.05,'tauR_E',0.4, 'tauD_E',10, 'E_netcon',E_netcon,'ESYN_E',0,...
                               'gSYN_I',.05,'tauR_I',0.4, 'tauD_I',10, 'I_netcon',LINetcon,'ESYN_I',-80};


%% vary
vary = {
    '(postLI->postLI)','gSYN_I',[0.005];
%   '(I->I,R->R)', 'SpatialChan', 3;
%   '(I->I,R->R)', 'freqInd', 1;
};
nVary = calcNumVary(vary);
parfor_flag = double(nVary > 1); % use parfor if multiple sims


%% simulate
data = dsSimulate(s,'time_limits',[dt time_end], 'solver',solverType, 'dt',dt,...
  'downsample_factor',1, 'save_data_flag',1, 'save_results_flag',1,...
  'study_dir',study_dir, 'debug_flag',1, 'verbose_flag',0,...
  'vary',vary,'parfor_flag',parfor_flag);

%% insert spikes
% V_spike = 50;
% for iData = 1:length(data)
%   for pop = {'I','R','C'}
%     pop = pop{1};
%     data(iData).([pop '_V'])(data(iData).([pop '_V_spikes']) == 1) = V_spike; % insert spike
%   end
% end

%% plot
% dsPlot(data(1:5));
% dsPlot(data(1:5), 'variable','C_V');
%  dsPlot(data(1), 'plot_type','raster');

% %% firing rates
% % calc fr
data2 = dsCalcFR(data, 'bin_size', 100, 'bin_shift',20); % 100ms bin, 20% shift
% 
% figure
% 
% % C firing rate
% sp(1) = subplot(211);
% if nVary == 1
%   plot(data2.time_FR, data2.postLI_V_spikes_FR);
%   xlabel('Time [ms]')
%   ylabel('Firing Rate [Hz]')
%   title('C Firing Rate')
% else
%   cFR = cat(2, data2.postLI_V_spikes_FR);
%   
%   % plot
%   imagesc(data2(1).time_FR, [], cFR');
% %   caxis([0 70])
%   colorbar
%   xlabel('Time [ms]')
%   ylabel('Freqs')
%   title('C Firing Rate')
% end
% 
% % IC input
% sp(2) = subplot(212);
% % freqInd = vary{strcmp(vary(:,2),'freqInd'), 3};
% freqInd = 3;
% if nVary == 1
%   icSpikes = logical(spk_IC(:,:,freqInd)'); % take all inputs for this freq
% else
%   icSpikes = logical(squeeze(spk_IC(:,3,:))'); % take only middle input for all freqs
% end
% plotSpikeRasterFs(icSpikes, 'PlotType','vertline', 'AxHandle',sp(2), 'Fs',spk_IC_fs);
% title('IC Input')
% ylabel('Freqs')
% 
% if nVary > 1
%   colorbar
% end
% 
% linkaxes(sp, 'x');

%%
% temp = [data.postIC_V_spikes];
% length(find(temp)) %number of spikes present

%%
addpath('eval_scripts')
if numel(size(spk_IC)) == 3
    temp = spk_IC(:,:,3);
elseif numel(size(spk_IC)) == 2
    temp = spk_IC;
end
maskIC = calcSpkMask(temp,40000,'alpha',.2);
figure;
imagesc(maskIC); title('IC mask');

for jj = 1:nVary
    mask(jj).Pre = calcSpkMask(data(jj).preLI_V_spikes,40000,'alpha',.01);
    mask(jj).Post = calcSpkMask(data(jj).postLI_V_spikes,40000,'alpha',.01);
end

%%
for jj = 1:nVary
%     subplot(3,3,jj)
    figure(1);imagesc(mask(jj).Pre); title('preLI mask');
    figure(2);imagesc(mask(jj).Post); title(sprintf('postLImask, gSYN_I:%0.3f',data(jj).postLI_postLI_gSYN_I));
end
%%
disp(['preLI spikes: ' num2str(sum(sum(data(jj).preLI_V_spikes)))]);
disp(['postLI spikes: ' num2str(sum(sum(data(jj).postLI_V_spikes)))]);

% figure;
% subplot(121);imagesc(data.preLI_V'); title('Voltages, pre')
% subplot(122);imagesc(data.postLI_V'); title('Voltages, post');


