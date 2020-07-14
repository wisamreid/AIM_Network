function sigIn = genICinput(tauR, tauD, dt, trial)
% input:
%   tauR, tauD = rise and fall times of the EPSP waveform
%   dt = sampling frequency of IC data
% output:
%   sigIn = a matrix of EPSP waveforms, with dimensions time x (nfreqs x nlocs)
%
% @ Erik Roberts, Kenny Chou
% Boston Univeristy, 2019
%
% TODO: allow resampling of spk_IC to different dt

IC_file = sprintf('IC_spks_t%02i.mat',trial);
if exist(IC_file,'file')
    fileData = load(IC_file,'spk_IC');
else
    fileData = load(['..' filesep IC_file],'spk_IC');
end

if iscell(fileData.spk_IC)
    fs = 44100;
    spk_IC = spkTime2Train(fileData.spk_IC,fs,2*fs);
else
    spk_IC = fileData.spk_IC;
end

% IC data should have the format of time x freqInd x location
% Permute data if this is misshaped.
if size(spk_IC,2) < size(spk_IC,3)
    spk_IC = permute(spk_IC,[1,3,2]);
end
timeLen = size(spk_IC,1);
nfreqs = size(spk_IC,2);
nlocs = size(spk_IC,3);
% Reshape data into a 2D matrix that is (time) x (frequencies x locations)
sigIn = reshape(spk_IC,timeLen,nfreqs*nlocs);


% ========================= create EPSP waveform =========================
% ============= convolve each time series with epsc waveform ==============
t_a = max(tauR,tauD)*7; % Max duration of syn conductance
t_vec = 0:dt:t_a;
tau2 = tauR;
tau1 = tauD;
tau_rise = tau1*tau2/(tau1-tau2);
b = ((tau2/tau1)^(tau_rise/tau1) - (tau2/tau1)^(tau_rise/tau2)^-1); % /tau2?
epsc =  - b * ( exp(-t_vec/tau1) - exp(-t_vec/tau2) ); % - to make positive
sigIn = conv2(sigIn,epsc','same');

end