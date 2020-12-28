function [nf,cf,bandwidth] = getFreqChanInfo(type,nf,low_freq,high_freq)
% get the center freqs and bandwidths of the channels

switch type
    case 'erb'
        if ~exist('high_freq','var')
            high_freq = 5000;
            disp('using 5kHz as ceiling of filterbank')
        end
        cf = ERBSpace(low_freq,high_freq,nf);
        bandwidth = 24.7*(4.37*cf+1)/1000;
    case 'nsl'
        shft = log2(2.5);
        cf = 440 * 2 .^ ((-31:97)/24 + shft); %cfs used by the NSL filterbank
        cf = cf(1:nf); %using the first 40 channels
        bandwidth = 24.7*(4.37*cf+1)/1000;
    otherwise
        error('type not supported')
end