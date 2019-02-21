% calls the cortical network model. 
% Model is driven by the spike-trains from the IC model, generated previously.


% ==== network parameters =================================
beamLoc = 3;
attendLoc = 3; %3 = 0 deg az
I2_inh = 0; % mA
networkParams = V26_set_parameters(beamLoc,attendLoc,I2_inh);
networkParams.nCortical = 1; %number of cortical neurons in the cortical network model
% networkParams.q = 25; %frequency tuning, Q
% =========================================================

%make backup of current script
FileName=mfilename;
day = today('datetime');
newbackup=sprintf('%s%s_backup_%s.m',saveLoc,mfilename,day);
currentfile=strcat(FileName, '.m');
copyfile(currentfile,newbackup);


load([dataLoc 'IC_spks.mat'],'spk_IC','freqGainNorm','input_gain','cf','fcoefs');
networkParams.cf = cf;
networkParams.q = inf;

trial_id = sprintf('trial_name');
% Cortical Network model
plot_on = 0;
tic
[spk_network26, spk_relay26] = V262RunNetwork(spk_IC, 40000, networkParams, plot_on,sprintf('%s%s',[saveLoc 'v26 '],trial_id));
toc

