function spks = spkTime2Train(spkTime,fs,timeLen)
% Converts cell array of spike times into a matrix of spike trains
%
% inputs:
%   spkTime = cell array of spike times. One cell for each freq channel
%   fs      = sampling freq
%   timeLen = length of time vector. Optional

if ~exist('timeLen','var')
    nTime = max(max(cellfun(@max,spkTime)))*fs; %find max time vector length
else
    nTime = timeLen;
end

spks = zeros(nTime,size(spkTime,1),size(spkTime,2));
for i = 1:size(spkTime,1)
    for j = 1:size(spkTime,2)
        spks(round(spkTime{i,j}*fs),i,j) = 1;
    end
end