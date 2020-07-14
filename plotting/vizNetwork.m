function vizNetwork(s)
% vizNetwork(s)
% s is the structure containing information for dynasim simulations
% assumes there are 64 frequencys, source and sink of graph are 'C' and 'R'
% respectively
%
% @ Kenny F Chou, BU
% 2020-06-09

nFreqs = 64;

% calculate indices for populating adjacency matrix
numCons = length(s.connections);
numPops = length(s.populations);
popLabels = {s.populations.name};
adjuMtx = zeros(sum([s.populations.size])+nFreqs);
currentIdx = 1;
for i = 1:numPops
    adjMtxIdx(i).pop = popLabels{i};
    adjMtxIdx(i).size = s.populations(i).size;
    adjMtxIdx(i).start = currentIdx;
    adjMtxIdx(i).end = currentIdx+adjMtxIdx(i).size-1;
    currentIdx = adjMtxIdx(i).end+1;
end

% populate adjacency matrix with netcons   
for i = 1:numCons
    direction = s.connections(i).direction;
    popPrePost = strsplit(direction,'->');
    popPre = popPrePost{1};
    popPost = popPrePost{2};
    popPreLabelIdx = find(ismember(popLabels,popPre));
    popPostLabelIdx = find(ismember(popLabels,popPost));
    popPreStart = adjMtxIdx(popPreLabelIdx).start;
    popPreEnd = adjMtxIdx(popPreLabelIdx).end;
    popPostStart = adjMtxIdx(popPostLabelIdx).start;
    popPostEnd = adjMtxIdx(popPostLabelIdx).end;
    
    userDefinedNetcon = ismember(s.connections(i).parameters(1:2:end),'netcon');
    if ~userDefinedNetcon
        % netcon matrix is ones(size_pre,size_post)
        disp('use zeros')
        nPre = s.populations(popPreLabelIdx).size;
        nPost = s.populations(popPostLabelIdx).size;
        adjMtx(popPreStart:popPreEnd,popPostStart:popPostEnd) = zeros(nPre,nPost);
    else
        % find netcon matrix within the parameters field
        disp('use netcon')
        netConIdx = find(userDefinedNetcon)*2; %hacky; but works
        netCon = s.connections(i).parameters{netConIdx};
        adjMtx(popPreStart:popPreEnd,popPostStart:popPostEnd) = netCon;
    end
end

% labels
clear names
namesIdx = 1;
for p = 1:numPops
    popName = adjMtxIdx(p).pop;
    popSize = adjMtxIdx(p).size;
    if popSize < nFreqs
        l = 1;
        for f = 1:popSize
            names{namesIdx} = sprintf('%s_%i_%i',popName,l,f);
            namesIdx = namesIdx + 1;
        end
    else
        for l = 1:popSize/nFreqs
            for f = 1:nFreqs
                names{namesIdx} = sprintf('%s_%i_%i',popName,l,f);
                namesIdx = namesIdx + 1;
            end
        end 
    end
end

% display adj matrix
figure;
imagesc(adjMtx);
title('Consolidated NetCon Matrix')
xlabel('post-synaptic population')
ylabel('pre-synaptic population')
xticks([adjMtxIdx.start]+[adjMtxIdx.size]/2)
xticklabels(popLabels)
set(gca,'xaxisLocation','top')
yticks([adjMtxIdx.start]+[adjMtxIdx.size]/2)
yticklabels(popLabels)
for i = 1:numPops-1
    line([1,size(adjMtx,2)],[adjMtxIdx(i).end+0.5 adjMtxIdx(i).end+0.5],'Color', 'r', 'LineWidth', 0.5);
    line([adjMtxIdx(i).end+0.5 adjMtxIdx(i).end+0.5], [1,size(adjMtx,1)],'Color', 'r', 'LineWidth', 0.5);
end

% graph
sources = 'C';
sinks = 'E';
sourcesIdx = find(ismember(popLabels,sources));
sinksIdx = find(ismember(popLabels,sinks));
sourceStart = adjMtxIdx(sourcesIdx).start;
sourceEnd = adjMtxIdx(sourcesIdx).end;
sinkStart = adjMtxIdx(sinksIdx).start;
sinkEnd = adjMtxIdx(sinksIdx).end;

g = digraph(adjMtx,names);
figure; p = plot(g,'layout','layered','Marker','o',...
                   'MarkerSize',8,'NodeLabel',names);
set(p,'ArrowSize',12)
layout(p,'layered','direction','down',...
         'sources',[sourceStart:sourceEnd],'sinks',[sinkStart:sinkEnd]);
