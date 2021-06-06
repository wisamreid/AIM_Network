% create sets of stimuli for top-down circuit demos
% talker1 @ 0 deg azimuth, from the CRM corpus
% talker2 @ 15:15:90 deg azimuth, from CRM corpus
% 2021-05-18

% paths
crmpath = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\toolboxes\CRM';
hrtfpath = 'C:\Users\Kenny\Desktop\GitHub\BOSSA\HRTF_40k';
addpath(crmpath);

% names of CRM files
load([crmpath filesep 'CRMnumber_3source_40pairs_Nov-22_18.mat'],'crm_number')
fs = 40000;

% randomly select talkers for each trio - use MFM or FMF scheme for
% greatest separataion. M/F talkers are numbered 0-3/4-7
talker1 = randi([0 3],20,1);
talker2 = randi([4 7],20,1);
talker3 = randi([0 3],20,1);

% load impulse responses
talker1az = 0;
talker2az = 15:15:90;
talker_azs = [talker1az talker2az];
numLocations = length(talker_azs);

for i = 1:numLocations
    hrtf_name = ['kemar_small_horiz_' num2str(talker_azs(i)) '_0.mat'];
    hrir = load([hrtfpath filesep hrtf_name]);
    hrtfL(:,i) = hrir.hrir_left;
    hrtfR(:,i) = hrir.hrir_right;
end

% spatialize CRM sentences with HRIRs
talkerKey = [1 ones(size(talker2az))*2];
for i = 2:20
    % load CRM sentences
    signals{j} = readbin([crmpath filesep 'Talker ' num2str(talker1(i)) filesep crm_number{i,1}]);
    for j = 2:numLocations
        signals{j} = readbin([crmpath filesep 'Talker ' num2str(talker2(j)) filesep crm_number{i,2}]);
    end
    
    % pad with zeros to make them the same length
    lens = cellfun(@length,signals);
    maxlen = max(lens);
    for j = 1:numLocations
        x = maxlen - lens(j);
        signals{j} = [signals{j} ; zeros(x,1)];
    end
    C = cell2mat(signals);
    
    % apply impulse responses
    sentencesL = fftfilt(hrtfL,C);
    sentencesR = fftfilt(hrtfR,C);
    
    wavs.tgt = [sentencesL(:,1) sentencesR(:,1)];
    for j = 1:numLocations-1
        wavs.(['m' num2str(talker2az(j))]) = [sentencesL(:,j+1), sentencesR(:,j+1)];
        wavs.(['mixed' num2str(talker2az(j))]) = [sentencesL(:,1)+sentencesL(:,j+1),...
                                                     sentencesR(:,1)+sentencesR(:,j+1)];
    end
    save(['stimuli' filesep 'CRM Stimuli TM Xdeg' filesep 'set' num2str(i,'%02i') '.mat'],'wavs','talker_azs','talker1','talker2','fs');
end
