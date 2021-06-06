% 2 talker CRM stimuli for harmonic attention network
% monaural waveforms; assume co-location

% paths
crmpath = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\toolboxes\CRM';
hrtfpath = 'C:\Users\Kenny\Desktop\GitHub\BOSSA\HRTF_40k';
addpath(crmpath);

% names of CRM files
load([crmpath filesep 'CRMnumber_3source_40pairs_Nov-22_18.mat'],'crm_number')
fs = 40000;

% Select M/F talkers for the simulation. M/F are numbered 0-3/4-7
talkers = [randi([0 3],20,1), randi([4 7],20,1)];
numtalkers = 2;

% load CRM sentences
wavs = struct();
for i = 1:20
    for j = 1:numtalkers
    signals{j} = readbin([crmpath filesep 'Talker ' num2str(talkers(i,j)) filesep crm_number{i,j}]);
    end

    % pad with zeros to make them the same length
    lens = cellfun(@length,signals);
    maxlen = max(lens);
    for j = 1:numtalkers
    x = maxlen - lens(j);
    signals{j} = [signals{j} ; zeros(x,1)];
    end
    C = cell2mat(signals);

    wavs(i).M = [C(:,1) C(:,1)];
    wavs(i).F = [C(:,2) C(:,2)];
    wavs(i).mixed = sum(C,2);
end
save('stimuli\CRM Stimuli TM 0deg.mat','wavs','talkers','crm_number');
