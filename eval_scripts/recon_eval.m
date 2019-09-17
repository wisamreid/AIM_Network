function [st,rstim] = recon_eval(spks,target_wav,mix_wav,params)
% performs stimulus reconstruction and objective intelligibility assessment on a set of spike trains
% calculates TF masks. The TF masks are normalized, then directly applied to mix_wav.
% out = recon_eval(spks,target_wav,target_spatialized,mix_wav,params)
% Inputs:
%	spks: a set of spike trains (from the IC for example)
%	target_wav_loc: full path to target.wav (str), or its waveform vector
%       can be single or dual channel. The waveform to be referenced
%       against.
%	mix_wav_loc: full path to mixed.wav, or its waveform vector. The mixed
%       waveform on which spkMask is applied.
%   params: structure with fields
%       .fs - mandatory
%       .tau - mask computation kernel time constant
%       .spatialChan - spatial channel to be reconstructed
%       -------- option 1: ------
%       .cf,
%       .coefs, or
%       -------- option 2: ------
%       .frgain
%       .low_freq
%       .high_freq
%       .numChannel
%       ----- other options -----
%       .maskRatio
%       .type: an array containing the reconstruction conditions
%              see outputs -> rstim for a list of conditions
% Outputs:
%	out: STOI scores of three reconstruction methods
%       ['filt' 'env' 'voc' 'mix' 'mix+pp']
%   rstim: a structure of reconstructed stimuli, with the following fields
%       fields returned depends on the conditions (x) given in the input
%       (1) .r1d = rstim1dual; FRmask
%       (1) .r1m = rstim1mono;
%       (2) .r2d = rstim2dual; FRmask + tone vocode
%       (2) .r2m = rstim2mono;
%       (3) .r3 = rstim3;      Mixed - Vocoded Mask
%       (4) .r4d = rstim4dual; FRmask + mixed vocode
%       (4) .r4m = rstim4mono;
%       (4) .r4pp = rstim4pp; (post-processed)
%
%
% @Kenny Chou
% 2019-04-30
% Boston University
% 20190726 KFC added delay parameter to calcSpkMask
% 20190916 removed rms scaling, removed extra target parameter
%% set up reference
if isa(target_wav,'char')
    target = audioread(target_wav);
    mixed = audioread(mix_wav);
else
    target = target_wav;
    mixed = mix_wav;
end

%frequency-filtered target & mixtures
if ~isfield(params,'fcoefs')
    fs = params.fs;
    frgain = params.frgain;
    low_freq = params.low_freq; %min freq of the filter
    high_freq = params.high_freq;
    numChannel = params.numChannel;
    [nf,cf,~] = getFreqChanInfo('erb',numChannel,low_freq,high_freq);
    fcoefs=MakeERBFilters(fs,cf,low_freq);
else
    fcoefs = params.fcoefs;
    cf = params.cf;
    nf = length(cf);
    fs = params.fs;
    frgain = 1;
end

mixedFiltL = ERBFilterBank(mixed(:,1),fcoefs);
mixedFiltR = ERBFilterBank(mixed(:,2),fcoefs);

%% mask calculation
maskParam.kernel = 'alpha';
maskParam.tau = params.tau;
maskParam.delay = params.delay;
masks = calcSpkMask(spks,fs,maskParam);
masksNorm = zeros(size(masks));
if ndims(masks) == 3
    % normalize masks
    for i = 1:size(masks,3)
        masksNorm(:,:,i) = masks(:,:,i)/max(max(masks(:,:,i)));
    end
    
    % remove side channels
    centerM = masksNorm(:,:,3);
    rightM1 = masksNorm(:,:,5);
    rightM2 = masksNorm(:,:,4);
    leftM1 = masksNorm(:,:,1);
    leftM2 = masksNorm(:,:,2);
    if params.spatialChan == 3
        spkMask = centerM-params.maskRatio.*(rightM1+leftM1+leftM2+rightM2);
        spkMask(spkMask<0)=0;
    else
        error('only the difference mask for the center spatial channel is implemented');
    end
else
    spkMask = masks;
end
spkMask = spkMask./max(max(abs(spkMask))); %normalize to [0,1]

%% Reconstruction
rstim = struct();
st = struct();
% ---------------- mask-filtered reconstruction --------------
if ismember(1,params.type)
    [rstim1dual, maskedWav] = applyMask(spkMask,mixedFiltL,mixedFiltR,frgain,'filt');
    st1 = runStoi(rstim1dual,target,fs,fs);
    
    rstim.r1d = rstim1dual;
    rstim.tf = maskedWav;
    st.r1 = st1;
end

% -------------- vocoded mask-filtered reconstruction ----------------
if ismember(2,params.type)
    %extract envelopes
    mixedEnvL = zeros(nf,length(mixedFiltL));
    mixedEnvR = zeros(nf,length(mixedFiltR));
    for i = 1:nf
        mixedEnvL(i,:)= envelope(mixedFiltL(i,:)); %envlope of mixture
        mixedEnvR(i,:)= envelope(mixedFiltR(i,:)); %envlope of mixture
    end
    
    % perform reconstruction
    [rstim2dual, ~] = applyMask(spkMask,mixedEnvL,mixedEnvR,frgain,'env',cf);
    st2 = runStoi(rstim2dual,target,fs,fs);
    
    rstim.r2d = rstim2dual;
    st.r2 = st2;
end

% -------------- vocoded spike-mask reconstruction ----------------
if ismember(3,params.type)
    % perform reconstruction
    fcutoff = 8000;
    [rstimTone,rstim3t] = vocode(spkMask,cf,'tone');
    [~,rstim3n] = vocode(spkMask,cf,'noise',fs);
    rstim3 = rstim3t;
    rstim3(cf>fcutoff,:) = rstim3n(cf>fcutoff,:);
    rstim3 = sum(rstim3);
    st3 = runStoi(rstim3,target,fs,fs);

    rstim.r3m = rstim3;
    rstim.r3t = rstimTone;
    st.r3 = st3;
end

% ------------- mixed vocoding of filtered mixture envelope ----------
if ismember(4,params.type)
    %extract envelopes
    mixedEnvL = zeros(nf,length(mixedFiltL));
    mixedEnvR = zeros(nf,length(mixedFiltR));
    for i = 1:nf
        mixedEnvL(i,:)= envelope(mixedFiltL(i,:)); %envlope of mixture
        mixedEnvR(i,:)= envelope(mixedFiltR(i,:)); %envlope of mixture
    end
    
    % mixed vocoding of filtered mixture envelope
    [rstim4dual, rstim4mono] = applyMask(spkMask,mixedEnvL,mixedEnvR,frgain,'mixed',cf);
    st4 = runStoi(rstim4mono,target,fs,fs);

    rstim4pp = runF0(rstim4dual,fs);
    st4pp = runStoi(rstim4pp,target,fs,fs);
    
    rstim.r4d = rstim4dual;
    rstim.r4m = rstim4mono;
    rstim.r4pp = rstim4pp;
    st.r4 = st4;
    st.r4pp = st4pp;
end

disp('eval complete')
