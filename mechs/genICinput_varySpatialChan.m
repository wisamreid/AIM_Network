function sigIn = genICinput_varySpatialChan(SpatialChan, tauR, tauD, dt)

% TODO: allow resampling of spk_IC to different dt

try
  fileData = load('IC_spks.mat','spk_IC');
catch
  fileData = load(fullfile('..', 'IC_spks.mat'),'spk_IC');
end

if ndims(fileData.spk_IC) == 3
    sigIn = fileData.spk_IC(:,:,SpatialChan); % time x freqChan x spatialChan
else %2D data, time x freq chan, ignores spatialChan
    sigIn = fileData.spk_IC;
end

%% create EPSP waveform
t_a = max(tauR,tauD)*7; % Max duration of syn conductance
t_vec = 0:dt:t_a;
tau2 = tauR;
tau1 = tauD;
tau_rise = tau1*tau2/(tau1-tau2);
b = ((tau2/tau1)^(tau_rise/tau1) - (tau2/tau1)^(tau_rise/tau2)^-1); % /tau2?
epsc =  - b * ( exp(-t_vec/tau1) - exp(-t_vec/tau2) ); % - to make positive

for iNeuron = 1:size(sigIn, 2)
  sigIn(:,iNeuron) = conv(sigIn(:,iNeuron), epsc, 'same');
end

end