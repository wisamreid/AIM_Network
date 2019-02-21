function mask = calcSpkMask(spks,fs,kernel,varargin)
% mask = calcSpkMask(spks,fs,kernel,varargin)
%
% Convolves spikes SPKS with an alpha function, which has parameter tau and
% fs. Default value for tau is 0.01 (10 ms). SPKS is assumed to have the
% dimensions [time x freq x neurons]
%
% output is [freq x time]. Dimensions are swapped for easy plotting
%
% Kenny F Chou
% updated 2018/9/18

switch kernel
    case 'alpha'
        if nargin < 4, tau = 0.01; else, tau = varargin{1}; end
        t = 0:1/fs:100/1000;
        A = t.*exp(-t/tau);
    case 'hamming'
        A = hamming(2000); %2000/fs = 0.05 second duration
end

numNeurons = size(spks,3);
tempc = zeros(size(spks));
for n = 1:numNeurons
    tempc(:,:,n) = conv2(spks(:,:,n),A','same');
end
mask = tempc';