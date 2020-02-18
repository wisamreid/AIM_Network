function [spk_IC, fs, fcoefs, cf] = prepICinput(scenario,talkerSet,fileloc,study_dir)
% load data, structure to correct format and save to study dir


switch scenario
    case 'consecutive'
%         fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\broad mode - staggered talkers 64Chan500-20000hz\CRM talker4';
        chans = {'00','90','-90'};
        refLabels = {'S2','S3','S1'};
        swav123 = [];
        spk123 = [];
        for i = 1:3
            ic(i).wav = audioread([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{i} 'degAz_target_conv.wav']);
            ic(i).data = load([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{i} 'degAz_SpkIC.mat']);
            ic(i).spksbin = spkTime2Train(ic(i).data.spk_IC,ic(i).data.fs,length(ic(i).wav));
            swav123 = [swav123;ic(i).wav];
            spk123 = [spk123;ic(i).spksbin];
        end
        spk_IC = spk123;
        target = swav123;
        mixed = target;
        save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
        fs = ic(1).data.fs;
        fcoefs = ic(1).data.fcoefs;
        cf  = ic(1).data.cf;
    case 'singular'
%         fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\broad mode - staggered talkers 64Chan500-20000hz\CRM talker4';
        chans = {'00'};
        tgtwav = audioread([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{1} 'degAz_target_conv.wav']);
        spkData = load([fileloc filesep '0_masker_set_' num2str(talkerSet) '_pos_' chans{1} 'degAz_SpkIC.mat']);
        spk_IC = spkTime2Train(spkData.spk_IC,spkData.fs,length(tgtwav));
        save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
        fs = spkData.fs;
        fcoefs = spkData.fcoefs;
        cf  = spkData.cf;
    case 'simultaneous'
%         fileloc = 'Z:\eng_research_hrc_binauralhearinglab\kfchou\ActiveProjects\Top-Down AIM network\Sims\simultaneous talkers 64Chan500-20000hz\CRM talker4';
        mixedName = ls([fileloc filesep sprintf('*set_%02i*mixed.wav',talkerSet)]);
        mixed = audioread([fileloc filesep mixedName]);
        spkName = ls([fileloc filesep sprintf('*set_%02i*SpkIC.mat',talkerSet)]);
        spks = load([fileloc filesep spkName]);
        if ~isfield(spks,'fs'), spks.fs = 40000; end
        spk_IC = spkTime2Train(spks.spk_IC,spks.fs,length(mixed));
        refNames = {'target_conv','masker1_conv','masker2_conv',};
        refLabels = {'S2','S3','S1'};
        for i = 1:3
            temp = ls([fileloc filesep sprintf('*set_%02i*%s.wav',talkerSet,refNames{i})]);
            ref(i).wav = audioread([fileloc filesep temp]);
        end
        tgt = ref(1).wav;
        fs = spks.fs;
        fcoefs = spks.fcoefs;
        cf  = spks.cf;
        save(fullfile(study_dir, 'solve', 'IC_spks.mat'),'spk_IC')
end

end