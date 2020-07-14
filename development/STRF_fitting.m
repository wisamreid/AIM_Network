x = 0:0.001:16; %kHz
center = 8.1;
sigma = 1;
g=gaussmf(x,[sigma,center]);

plot(log2(x),g)
xlim([0,4])
view(90,-90)

addpath('C:\Users\Kenny\Dropbox\Sen Lab\m-toolboxes\plotting utils')
addpath('C:\Users\Kenny\Dropbox\Sen Lab\m-toolboxes\plotting utils\MatPlotLib2.0 Colormaps')
color1 = brewermap(3,'set1');

%% investigate a range of parameters
% parameter values to vary over
x = options.cf/1000; %unit in kHz
range = cell(4);
range{1} = 8.1; %x0
range{2} = [0.005:0.005:0.1]; %BTM
range{3} = [-1:.1:1]; %phase
range{4} = [0.125,0.25,0.5,1,2,4,8,16]; %BW
variedParam = {'t0','BTM','phase','BW'};
figure;
for j = 1:length(range)
    subplot(2,2,j);
    colormap = inferno(length(range{j}));
    i = 1;
    % default parameters:
    paramH.BW= 0.35; %bandwidth
    paramH.BTM= 0.01 ; %BTM, modulation
    paramH.t0= 3; % t0, peak latency (s)
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

%% visualize a single set of parameters
%xLin = 0:0.001:16; %kHz
%x = log2(xLin);
x = options.cf/1000; %unit in kHz
paramH.BW= 0.35; %bandwidth
paramH.BTM= 0.01 ; %BTM, modulation
paramH.t0= 3; % t0, peak latency (s)
paramH.phase= 0; % phase
hGauss = exp(-0.5*((x-paramH.t0)/paramH.BW).^2);
hcos = cos(2*pi*paramH.BTM*(x-paramH.t0)+paramH.phase*pi);
h = hGauss.*hcos;
plot(x,h)
xlim([0,4])
view(90,-90)