function h = netcon_spread(x,bw,type,f0Shift)
% returns netcon matrix H, where each frequency channel in X has a spread
% of BW octaves or BW kHz
%   x - frequency axis vector, in kHz. Center frequencies of each channel.
%   bw - bandwidth in octaves
%   type - 'q' for constant Q, or 'bw' for constant bandwidth

if ~exist('f0Shift','var'), f0Shift = 0; end
sigma = bw/2.3548;

switch type
    case 'q'
        i = 1;
        faxis = linspace(-4.3219,4,1500); % 0.05~16.0 kHz
        [~, closestIndex] = min(abs(2.^faxis'-x')); %map to ERB axis
        for f0 = x'
            hGauss = exp(-0.5*((faxis-log2(f0+f0Shift))/sigma).^2);
            h(i,:) = hGauss(closestIndex);
            i = i+1;
        end

    case 'bw'
        i = 1;
        for f0 = x'
            hGauss = exp(-0.5*((x-f0+f0Shift)/sigma).^2);
            h(i,:) = hGauss;  
            i = i+1;
        end
end