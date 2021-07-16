% create sets of stimuli for top-down circuit demos
% 3 talkers @ -90, 0 , and 90 deg az, from the CRM corpus
% 2020-05-26

% paths
crmpath = fullfile('..','CRM');
hrtfpath = fullfile('..','BOSSA','HRTF_40k');
addpath(crmpath);

% names of CRM files
load([crmpath filesep 'CRMnumber_3source_40pairs_Nov-22_18.mat'],'crm_number')
fs = 40000;

% randomly select talkers for each trio - use MFM or FMF scheme for
% greatest separataion. M/F talkers are numbered 0-3/4-7
talkers = [randi([0 3],20,1), randi([4 7],20,1), randi([0 3],20,1)];

% load impulse responses
talker_azs = [0 90];
numtalkers = length(talker_azs);
for i = 1:numtalkers
    hrtf_name = ['kemar_small_horiz_' num2str(talker_azs(i)) '_0.mat'];
    hrir = load([hrtfpath filesep hrtf_name]);
    hrtfL(:,i) = hrir.hrir_left;
    hrtfR(:,i) = hrir.hrir_right;
end

for i = 1:20
    % load CRM sentences
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
    
    % apply impulse responses
    sentencesL = fftfilt(hrtfL,C);
    sentencesR = fftfilt(hrtfR,C);
    
    wavs(i).tgt = [sentencesL(:,1) sentencesR(:,1)];
    wavs(i).m1 = [sentencesL(:,2) sentencesR(:,2)];
%     wavs(i).m2 = [sentencesL(:,3) sentencesR(:,3)];
    wavs(i).mixed = [sum(sentencesL,2) sum(sentencesR,2)];
    
end
save(fullfile('stimuli','CRM Stimuli TM 90deg.mat'),'wavs','talkers','talker_azs');