function i2InputCurrent = genI2input(networkDim,targetChans,simLen,fs,study_dir)
% out = genI2input(nFreq,nLocs,simLen,fs)
% input current to I2 cells; defaulted to be on
% control current amplitude with IappI2 parameter (in initI2.mech)
%     1 = on, 0 = off
% I2 input should have the format of time x freqInd x location
% input:
%   networkDim = [nFreqs, nLocs];
%   targetChans = {[freq1,loc1],[freq2,loc2]...}; channels to toggle current
% output: (time)x(nFreq*nLocs)
%
% @ Kenny Chou
% Boston Univeristy, 2019

if ~iscell(targetChans), warning('TARGETCHANS must be a cell array'); end

nFreq = networkDim(1);
nLocs = networkDim(2);
taps = simLen;
% dur = 1; %seconds

% define when input current to I2 will be turned off and on
i2input = ones(taps,nFreq,nLocs);
if ~cellfun('isempty',targetChans)
    chanSwitch = cell(nFreq,nLocs);
    for i = 1:numel(targetChans)
        freq = targetChans{i}(1);
        loc = targetChans{i}(2);
        chanSwitch(freq,loc) = {0}; %parameterize this later

        % toggle i2 input off/on at the appropriate index
        switchIdx = round(chanSwitch{freq,loc}*(fs/1000));
        for j = 1:length(switchIdx)
            i2input(switchIdx(j)+1:end,freq,loc) = abs(i2input(switchIdx(j)+1:end,freq,loc)-1);
        end
    end
end

% define time-varying function for input current here
% i2input(xtaps,freq,loc) = some function of t

% reshape to collapse dimensions
i2InputCurrent = reshape(i2input,taps,nFreq*nLocs);

save([study_dir filesep 'solve' filesep 'i2input.mat'],'i2InputCurrent')