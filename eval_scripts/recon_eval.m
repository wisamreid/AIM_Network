function [out,rstim1dual,rstim2dual,rstim3] = recon_eval(data,target_wav_loc,target_spatialized_loc,mix_wav_loc,params)
% performs stimulus reconstruction and objective intelligibility assessment on a set of spike trains
% Inputs:
%	data: output from DynaSim
%	target_wav_loc: file name (including its path) of the target wave (target.wav)
%	target_spatialized_loc: file name of the spatialized wave (target_conv.wav)
%	mix_wav_loc: file name of the sound mixture (target+maskers, mixed.wav)
% Outputs:
%	out: STOI scores of three reconstruction methods

temp = [data.C_V_spikes];

%% set up reference
target = audioread([target_wav_loc]);
target = target/max(abs(target));
targetLR = audioread(target_spatialized_loc);
targetLRmono = targetLR(:,1) + targetLR(:,2);
targetLRmono = targetLRmono/abs(max(targetLRmono));
mixed = audioread([mix_wav_loc]);

%frequency-filtered target & mixtures
if ~isfield(params,'fcoefs')
    fs = params.fs;
    frgain = params.frgain;
    low_freq = params.low_freq; %min freq of the filter
    high_freq = params.high_freq;
    numChannel = params.numChannel;
    [nf,cf,~] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
    fcoefs=MakeERBFilters(40000,cf,low_freq);
else
    fcoefs = params.fcoefs;
    cf = params.cf;
    nf = length(cf);
    fs = params.fs;
    frgain = 1;
end

targetFiltmono = ERBFilterBank(target,fcoefs);
mixedFiltL = ERBFilterBank(mixed(:,1),fcoefs);
mixedFiltR = ERBFilterBank(mixed(:,2),fcoefs);

clear mixedEnvL mixedEnvR targetEnv targetConvEnv
for i = 1:nf
    mixedEnvL(i,:)= envelope(mixedFiltL(i,:)); %envlope of mixture
    mixedEnvR(i,:)= envelope(mixedFiltR(i,:)); %envlope of mixture
    targetEnv(i,:) = envelope(targetFiltmono(i,:)); %envelope of target
end

%% reconstruction
spkMask = calcSpkMask(temp,fs,'alpha');

% mixture carrier reconstruction
[rstim1dual, rstim1mono] = applyMask(spkMask,mixedFiltL,mixedFiltR,frgain,'filt');
st1 = runStoi(rstim1mono,targetLRmono,fs,fs);
%apply to envelope of filtered mixture
[rstim2dual, rstim2mono] = applyMask(spkMask,mixedEnvL,mixedEnvR,frgain,'env',cf);
st2 = runStoi(rstim2mono,targetLRmono,fs,fs);

% Vocoded-SpikeMask
rstim3 = vocode(spkMask,cf,'tone');
st3 = runStoi(rstim3,target,fs,fs);

out = [st1 st2 st3];