function [out, reconMat] = vocode(estenv,cf,method,varargin)
% performs vocoding on the estimated 2D envelope ESTENV
% cf is an array of center frequencies
% method - either 'tone', 'noise', or 'carrier' is supported
%   
%   vocode(estenv,cf,'tone',fs);
%   fs = sampling frequency
%
%   vocode(estenv,cf,'tone',fs,interleave,cf_cutoff);
%   interleave: 0 or 1, implements interleaving based on Ghitza 2001;
%   returns a dichotic signal.
%   cf_cutff: specifies the lowpass frequency to not interleave
%
%   vocode(estenv,cf,'carrier',original,fcoefs);
%   original = original signal, used for calculating carrier
%   fcoefs = ERB filter coefficients
%
% Assumes ERB filterbank is used.
%
% Kenny Chou
% 12/19/2016
% Boston Univ. HRC
%
% 2017 - added 'carrier' support
% 20180508 - added interleaving option
% 20180525 - slight edits to 'carrier' option

nf = length(cf);
envLen = max(size(estenv));
recon = 0;
interleave = 0;
cutoff_freq_default = 1200;

if min(size(estenv)) ~= nf
    error('mismatch in the number of channels between CF and estEnv')
end

if size(estenv,2) ~= nf
    estenv = estenv';
end

P = inputParser;
addOptional(P,'interleave',interleave,@isnumeric);
addOptional(P,'cutoff_freq',cutoff_freq_default,@isnumeric);
parse(P,varargin{:})
interleave = P.Results.interleave;
cf_cutoff = P.Results.cutoff_freq;

switch method
    case 'tone'        
        for f = 1:nf
            t = (0:1/40000:(envLen-1)/40000)';
            carrier = sin(2*pi*cf(f)*t);
            envshapenoise = carrier.*estenv(:,f);
            reconMat(f,:) = envshapenoise;
        end
        
        if interleave
            reconE = zeros(size(reconMat));
            cf_cutoff_idx = find(cf<cf_cutoff,1);

            reconE(cf_cutoff_idx:end,:) = reconMat(cf_cutoff_idx:end,:);
            reconO = reconE;

            reconE(2:2:cf_cutoff_idx,:) = reconMat(2:2:cf_cutoff_idx,:);
            reconO(1:2:cf_cutoff_idx,:) = reconMat(1:2:cf_cutoff_idx,:);
            recon = [sum(reconE); sum(reconO)]';
        else
            recon = sum(reconMat)';
        end
        
    case 'noise'
        if nargin < 4
            error('must pass in fs')
        else
            fs = varargin{1};
        end
        temp = randn(envLen,1);
        fcoefs = MakeERBFilters(fs,cf,200);
        noise = ERBFilterBank(temp,fcoefs);
        for f = 1:nf
             reconMat(f,:) = noise(f,1:envLen)'.*estenv(:,f);
        end
        recon = sum(reconMat)';
    case 'carrier'
        if nargin < 5
            error('must pass in original sound mixture and fcoefs')
        else
            original = varargin{1};
            fcoefs = varargin{2};
        end
        filtered = ERBFilterBank(original,fcoefs);
        for kk = 1:nf
            senvs(kk,:) = envelope(filtered(kk,:),10000,'analytic'); %envelope of the original
            carrier(kk,:) = filtered(kk,:)./(senvs(kk,:)); %carrier of the original
            n1 = min(length(carrier(1,:)),length(estenv(:,1)));
            temp(kk,:) = carrier(kk,1:n1)' .* estenv(1:n1,kk);
        end
        reconMat = temp;
        recon = sum(temp);
end

out = recon;


    