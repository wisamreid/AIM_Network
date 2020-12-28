function plotIVTraces_E_E(data,fs,tgtChan,BFchan,dataFolder,varies,varied_param,xrange)
% data = dynasim output
% range = value of the varied parameter
% E-E from tgtChan to BFchan

if ~exist('xrange','var'), xrange = [0 length(data(1).IC_V(:,tgtChan))/fs]; end

t = (0:length(data(1).IC_V(:,tgtChan))-1)/fs;

h = figure('position',[450 200 1200 700]);
foldername = [dataFolder filesep 'chan' num2str(tgtChan)];
if ~exist(foldername,'dir'), mkdir(foldername); end

range = varies(varied_param).range;
for i = 1:length(data)
    clf; 
    subj(1).trace = data(i).E_V(:,tgtChan);
    subj(1).label = 'E Voltage (tgt chan)';
    subj(2).trace = data(i).E_E_dsSynapse_Isyn(:,BFchan);
    subj(2).label = 'E-E current (tgt->BF)';
    subj(3).trace = data(i).E_V(:,BFchan);
    subj(3).label = 'E voltage (BF chan)';
    subj(4).trace = data(i).E_XI_dsSynapse_Isyn(:,tgtChan);
    subj(4).label = 'Inhibition to E (tgt chan)';
%     subj(5).trace = spk_IC(:,tgtChan);
%     subj(5).label = 'IC model spikes';
    
    for j = 1:length(subj)
        nsubplots = length(subj);
        height = 0.9/nsubplots;
        width = 0.85;
        xpos = (1-width)/2;
        ypos = (1-height*nsubplots)/2;
        subplot('position',[xpos ypos+height*(nsubplots-j) width height*0.85])
        
        plot(t,subj(j).trace)
        ylabel(subj(j).label)
        set(gca,'fontsize',12)
        if j ~= nsubplots, xticks([]); end
        xlim(xrange)
    end
    suptitle([varies(varied_param).conxn ' ' varies(varied_param).param ' ' num2str(range(i))])
    saveas(gcf,[dataFolder filesep 'chan' num2str(tgtChan) filesep 'trial' num2str(range(i),'%.04f') '.tif'])

end


end