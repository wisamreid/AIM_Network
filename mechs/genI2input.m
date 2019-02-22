function out = genI2input(nFreq)
% input to I2 cells; defaulted to be on
%     1 = on, 0 = off
% I2 input should have the format of time x freqInd x location
% output: (time)x(nFreq*nLocs)
%
% @ Kenny Chou
% Boston Univeristy, 2019


% taps = size(fileData.spk_IC,1);
taps = 76000;
% dur = 1; %seconds
fs = 40000;
nLocs = 5;

%define when input current to I2 will be turned off and on
chanSwitch{1} = []; %ms
chanSwitch{2} = [];
chanSwitch{3} = [0];
chanSwitch{4} = [];
chanSwitch{5} = [];

i2input = ones(taps,nFreq,nLocs);
for i = 1:nLocs
    switchIdx = chanSwitch{i}*(fs/1000);
    for j = 1:length(switchIdx)
        i2input(switchIdx(j)+1:end,:,i) = abs(i2input(switchIdx(j)+1:end,:,i)-1);
    end
end
out = reshape(i2input,taps,nFreq*nLocs);