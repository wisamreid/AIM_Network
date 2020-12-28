% aritani2009_plotMarginals
% with E-E connection
load('data\036_Aritani_3a_ECBW1.1_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\036_Aritani_3a_ECBW1.1_Attend1\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

ylabel('total spike count')
xlabel('f0 (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))

%% without E-E connection
figure;
load('data\038_Aritani_3a_ECBW1.1_EE0_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\038_Aritani_3a_ECBW1.1_EE0_Attend1\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

ylabel('total spike count')
xlabel('f0 (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
