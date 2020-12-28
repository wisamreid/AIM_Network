% plot marginals Fritz2003, Atiani parameters
figure;

load('data\039_Fritz_AtianiParams_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\039_Fritz_AtianiParams_Attend1\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

ylabel('total spike count')
xlabel('f0 (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))

%% Atiani parameters, higher EC gsyn
figure;

load('data\040_Fritz_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\040_Fritz_Attend1\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

xlabel('cumulative firing rate')
ylabel('f0 (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))