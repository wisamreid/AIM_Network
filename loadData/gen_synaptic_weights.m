%calculate synaptic weights for each frequency channel, based on a gaussian
%distribution of frequency "leakage"
function weights = gen_synaptic_weights(qs,cf,plots)
% Q factors
% qs = [2,3,5,8,12,18,25,55];
%
% There are two ways to define Q: 
% Q = CF/sigma, or
% Q = CF/FWHM
%
% This function uses the first definition. The two definitons differ by a
% factor of 2sqrt(2ln(2)) = 2.355. FWHM = 2.355*sigma
% When labeling, Qs are divided by 2.355 to be consistent with the
% eNeuro manuscript.

if ~exist('plots','var'), plots = 0; end
if qs == inf
    temp = diag(ones(1,length(cf)));
    for i = 1:length(cf)
    	weights(i,1).weights = temp(i,:);
    end
    if plots
        figure;
        imagesc(reshape([weights.weights],length(cf),[]));
        title('weights for all CFs')
    end
    return;
end

%define a gaussian function
gaussian = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma.^2));

%calculate indices on the x axis corresponding to CFs
x = 0:0.1:5000;
for i = 1:length(cf)
    cfidx(i) = find(x==round(cf(i)));
end
temp = zeros(size(x));
temp(cfidx) = 1;


%% cf weighing calculation
clear gauss
for q = 1:length(qs)
    sigma = cf/qs(q);
    for i = 1:length(cf)
        gauss(i,:,q) = gaussian(x,cf(i),sigma(i));
        weights(i,q).weights = gauss(i,cfidx,q);
        weights(i,q).cf = cf(i);
        weights(i,q).q = qs(q);
    end
end

%% plotting
if plots
    i = 18; %18th cf
    figure;
    plot(x,temp,'color',[165,181,206]/255);
    hold on;
%     cmap=pmkmp(length(qs),'Edge');
    cmap = [0.9298    0.0789    0.2230;
    0.8119    0.2614         0;
    0.4608    0.4943         0;
    0.1023    0.6669         0;
         0    0.7030    0.2964;
    0.0206    0.6124    0.6525;
    0.2669    0.4759    0.8102;
    0.4989    0.3875    0.7428];
    clear gauss
    for q = 1:length(qs)
        sigma = cf/qs(q);
        gauss(i,:) = gaussian(x,cf(i),sigma(i));
        plot(x,gauss(i,:),'color',cmap(q,:),'linewidth',1.5);
    end
    set(gca, 'XScale', 'log')
    xlim([150 5000])
    xlabel('frequency (hz)')
    title(sprintf('gussian for various Qs, for CF = %d',round(cf(i))))
    temp = round(qs/2.355,2);
    legend(['cf'; cellstr(num2str((temp)'))])
    
    figure;
    imagesc(reshape([weights.weights],length(cf),[]));
    title('weights for all CFs')
    xlabel('CF of relay neurons, for different Qs')
    ylabel('weight of neighboring CF')
    xticks(length(cf)/2:length(cf):length(cf)*length(qs));
    xticklabels(temp)
end