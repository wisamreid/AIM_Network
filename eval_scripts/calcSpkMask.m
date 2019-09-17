function mask = calcSpkMask(spks,fs,param)
% mask = calcSpkMask(spks,fs,param)
%
% Convolves spikes SPKS with an alpha function, which has parameter tau and
% fs. Default value for tau is 0.01 (10 ms). SPKS is assumed to have the
% dimensions [time x freq x neurons]
%
% param: parameter structure, containing
%   .kernel
%   .tau (optional)
%   .delay (optional)
%
% output is [freq x time]. Dimensions are swapped for easy plotting
%
% Kenny F Chou
% updated 2018/9/18
% 20190726 KFC added delay parameter to calcSpkMask
%              charged varargin to param structure

if ~isfield(param,'tau'), tau = 0.01; else, tau = param.tau; end
if isfield(param,'delay'), delayLen = param.delay; else, delayLen = 0; end

switch param.kernel
    case 'alpha'
        t = 0:1/fs:100/1000;
        A = t.*exp(-t/tau);
    case 'hamming'
        A = hamming(2000); %2000/fs = 0.05 second duration
end

if delayLen > 0 
    A = [zeros(1,delayLen) A];
elseif delayLen < 0
    A = [A zeros(1,abs(delayLen))];
end

numNeurons = size(spks,3);
tempc = zeros(size(spks));
for n = 1:numNeurons
    tempc(:,:,n) = conv2(spks(:,:,n),A','same');
end
mask = permute(tempc,[2 1 3]);