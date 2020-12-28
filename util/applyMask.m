function [rstimDual,maskedWav] = applyMask(mask,wavL,wavR,frgain,method,cf,varargin)
% I wrote this function to reduce some lines of code. Assuming fs = 40000;
%
% inputs:
%   mask - any type of time-frequency mask, [channels x time]
%   wavL, wavR - time-frequency representation of signals
%   frgain - vector of frequency gains
%   method:
%       'filt' for spike-mask filtered stimulus. wavL and wavR should be
%       the two channels of s_filt, usually ERB filtered stimulus
%       'env' for vocoded spike-mask filtered stimulus. wavL and wavR
%       should be the envelopes of s_filt. 
%
% outputs:
%   spkMask
%   rstim1_dual, rstim1_mono: mask applied to S(t,f), with fine structure
%   rstim2_dual, rstim2_mono: mask applied to envelope of S(t,f), and vocoded
%
% Kenny Chou
% 4/17/2018
% NSNC Lab, BU Hearing Research Center
%
% 2018-10-16 - removed interleaving support
% 2019-09-16 - replaced rstimMono output with pre-summed freq-filtered output

if iscolumn(wavL), wavL = wavL'; end
if iscolumn(wavR), wavR = wavR'; end

n = min(length(mask),length(wavL));
% interleave = 0;
% cutoff_freq_default = 1200;

P = inputParser;
% addOptional(P,'interleave',interleave,@isnumeric);
% addOptional(P,'cutoff_freq',cutoff_freq_default,@isnumeric);
% parse(P,varargin{:})
% interleave = P.Results.interleave;
% cf_cutoff = P.Results.cutoff_freq;

maskedWavL = mask(:,1:n).*wavL(:,1:n)./frgain;
maskedWavR = mask(:,1:n).*wavR(:,1:n)./frgain;
maskedWav.L = maskedWavL;
maskedWav.R = maskedWavR;
switch method
    case 'filt' %spike-mask filtered stimulus
        %apply to filtered mixture
        rstimL = sum(maskedWavL);
        rstimR = sum(maskedWavR);
        rstimDual = [rstimL' rstimR'];
    case 'env' %vocoded spike-mask filtered stimulus
        %apply to envelope of filtered mixture
        rstimL = vocode(maskedWavL,cf,'tone');
        rstimR = vocode(maskedWavR,cf,'tone');
        rstimDual = [rstimL rstimR];
    case 'noise' %applies mixed/hybrid vocoding to the mask-modulated envelopes
        fs = 40000;
        rstimL = vocode(maskedWavL,cf,'noise',fs);
        rstimR = vocode(maskedWavR,cf,'noise',fs);
        rstimDual = [rstimL rstimR];
    case 'mixed'
        fs = 40000;
        fcutoff = 2000;
        [~,rstimToneL] = vocode(maskedWavL,cf,'tone');
        [~,rstimToneR] = vocode(maskedWavR,cf,'tone');
        [~,rstimNoiseL] = vocode(maskedWavL,cf,'noise',fs);
        [~,rstimNoiseR] = vocode(maskedWavR,cf,'noise',fs);
        rstimL = rstimToneL;
        rstimL(cf>fcutoff,:) = rstimNoiseL(cf>fcutoff,:);
        rstimL = sum(rstimL);
        rstimR = rstimToneR;
        rstimR(cf>fcutoff,:) = rstimNoiseR(cf>fcutoff,:);
        rstimR = sum(rstimR);
        rstimDual = [rstimL; rstimR]';
    otherwise
        disp('method not supported')
end
rstimMono = rstimL + rstimR;
rstimMono = rstimMono/max(abs(rstimMono));