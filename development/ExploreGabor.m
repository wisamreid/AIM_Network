% gabor functions
addpath('C:\Users\Kenny\Dropbox\Sen Lab\m-toolboxes\plotting utils')
addpath('C:\Users\Kenny\Dropbox\Sen Lab\m-toolboxes\plotting utils\MatPlotLib2.0 Colormaps')
color1 = brewermap(3,'set1');

x = options.cf/1000; %unit in kHz

% parameter values to vary over
range = cell(4);
range{1} = 8.9; %f0, kHz
range{2} = 0.001:0.05:1; %BFM, kHz/cycle?
range{3} = 0; %phase
range{4} = 0.2:0.2:1; %BW
variedParam = {'t_0','BTM','phase','BW'};
figure;
for j = 1:length(range)
    subplot(2,2,j);
    colormap = inferno(length(range{j}));
    i = 1;
    % default parameters:
    paramH.BW= 0.5; %bandwidth
    paramH.BTM= 0.01 ; %BTM, modulation
    paramH.t0= 8.9; % t0, peak latency (s)
    paramH.phase= 0; % phase

    for var = range{j}
        eval(['paramH.' variedParam{j} ' = var;']);
        hGauss = exp(-0.5*((x-paramH.t0)/paramH.BW).^2);
        hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
        h = hGauss.*hcos;
        plot(x,h,'color',colormap(i,:)); hold on;
        i = i+1;
    end
    hleg = legend(cellstr(num2str((range{j})')));
    title(variedParam{j})
end

% hGauss = exp(-0.5*((x-paramH.t0)/paramH.EC_BW).^2);
% hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
% h = hGauss.*hcos;
% paramH.EC_BW= .2; %bandwidth, octaves
% paramH.BTM= 0.01 ; %BTM, modulation
% paramH.t0= 8.9; % f0, peak f
% paramH.phase= 0; % phase

%%
% default parameters:
paramH.BW= 0.25; %bandwidth
paramH.BTM= 0.01 ; %BTM, modulation
paramH.t0= 8.9; % t0, peak latency (s)
paramH.phase= 0; % phase
hGauss = exp(-0.5*((x-paramH.t0)/paramH.BW).^2);
hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
h = hGauss.*hcos;
plot(x,h)
flat = conv(h,ones(1,100));
hold on;
plot(x,flat(1:length(x))/max(flat));