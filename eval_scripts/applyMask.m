function [rstimDual,rstimMono] = applyMask(mask,wavL,wavR,frgain,method,cf,varargin)
% I wrote this function to reduce some lines of code
%
% inputs:
%   mask - any type of time-frequency mask, [channels x time]
%   wavL, wavR
%   frgain - vector of frequency gains
%   method:
%       'filt' for spike-mask filtered stimulus. wavL and wavR should be
%       the two channels of s_filt, usually ERB filtered stimulus
%       'env' for vocoded spike-mask filtered stimulus. wavL and wavR
%       should be the envelopes of s_filt. If interleave is specified, the
%       function implements interleaving as described in Ghitza, 2001.
%
% outputs:
%   spkMask
%   rstim1_dual, rstim1_mono: mask applied to S(t,f), with fine structure
%   rstim2_dual, rstim2_mono: mask applied to envelope of S(t,f), and vocoded
%
% Kenny Chou
% 4/17/2018
% HRC at BU

n = min(length(mask),length(wavL));
% interleave = 0;
% cutoff_freq_default = 1200;

P = inputParser;
% addOptional(P,'interleave',interleave,@isnumeric);
% addOptional(P,'cutoff_freq',cutoff_freq_default,@isnumeric);
% parse(P,varargin{:})
% interleave = P.Results.interleave;
% cf_cutoff = P.Results.cutoff_freq;

switch method
    case 'filt' %spike-mask filtered stimulus
        %apply to filtered mixture
        extractedOriginalL = mask(:,1:n).*wavL(:,1:n)./frgain;
        extractedOriginalR = mask(:,1:n).*wavR(:,1:n)./frgain;
        rstimL = sum(extractedOriginalL);
        rstimR = sum(extractedOriginalR);
        rstimDual = [rstimL' rstimR'];
        rstimMono = rstimL + rstimR;
        rstimMono = rstimMono/max(abs(rstimMono));
    case 'env' %vocoded spike-mask filtered stimulus
        %apply to envelope of filtered mixture
        extractedEnvL = mask(:,1:n).*wavL(:,1:n)./frgain;
        extractedEnvR = mask(:,1:n).*wavR(:,1:n)./frgain;
        rstimL = vocode(extractedEnvL,cf,'tone');
        rstimR = vocode(extractedEnvR,cf,'tone');
        rstimDual = [rstimL rstimR];
        rstimMono = rstimL + rstimR;
        rstimMono = rstimMono/max(abs(rstimMono));
    case 'maskOriginalPhase'
        
    otherwise
        disp('method not supported')
end
            

