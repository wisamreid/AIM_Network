% aritani2009_plotMarginals
figure;

load('data\037_Aritani_3b_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\037_Aritani_3b_Attend1\STRF.mat')
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

load('data\037_Aritani_3b_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\038_Aritani_3b_EE0_Attend1\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

ylabel('total spike count')
xlabel('f0 (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))

%% matching cnx strength of 3a (no E-E)
figure;

load('data\038_Aritani_3b_EE0_test_Attend0\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

load('data\038_Aritani_3b_EE0_test_Attend1\STRF.mat')
plot(sum(runningFR,2),log2(f0/1000)); hold on;

xlabel('total spike count')
ylabel('stimulus frequency (Hz)')
legend('no attend','attend');
yMinMax = [min(ySpacing) max(ySpacing)];
ylim(yMinMax)
yticks(ySpacing)
yticklabels(round(2.^ySpacing,1))
