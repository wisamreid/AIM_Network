% vizNetwork
% pops = {data.model.specification.populations.name};
% sz = {data.model.specification.populations.size};
%
%
% 2019-03-06 added IC->I and IC->R connections

% construct cell array of names
clear names
i2Size = s.populations(4).size;
for l = 1:nLocs
    for f = 1:nFreqs
        i2Names{f+(l-1)*nFreqs} = sprintf('I2_%i_%i',l,f);
    end
end 
iSize = s.populations(1).size;
for l = 1:nLocs
    for f = 1:nFreqs
        iNames{f+(l-1)*nFreqs} = sprintf('I_%i_%i',l,f);
    end
end
rSize = s.populations(2).size;
for l = 1:nLocs
    for f = 1:nFreqs
        rNames{f+(l-1)*nFreqs} = sprintf('R_%i_%i',l,f);
    end
end
cSize = s.populations(3).size;
for l = 1
    for f = 1:nFreqs
        cNames{f+(l-1)*nFreqs} = sprintf('C_%i_%i',l,f);
    end
end
icSize = cSize;
for l = 1
    for f = 1:nFreqs
        icNames{f+(l-1)*nFreqs} = sprintf('IC_%i_%i',l,f);
    end
end
names2 = [icNames i2Names iNames rNames cNames]; %order matters!

%% IC NetCon to I and R neurons, not included in s.populations

icrNetconCell = repmat({ones(1,5)},1,nFreqs);
icrNetconGroupByLoc = blkdiag(icrNetconCell{:});
icrNetcon = regroup(icrNetconGroupByLoc', [nFreqs,nLocs])';
iciNetconCell = repmat({ones(1,5)},1,nFreqs);
iciNetconGroupByLoc = blkdiag(iciNetconCell{:});
iciNetcon = regroup(iciNetconGroupByLoc', [nFreqs,nLocs])';
%% combine connectivity matrices, in the order of I2,I,R,C
icStart = 1;
icEnd = icStart+icSize-1;
i2Start = icEnd+1;
i2End = i2Start+i2Size-1;
iStart = i2End+1;
iEnd = iStart+iSize-1;
rStart = iEnd+1;
rEnd = rStart+rSize-1;
cStart = rEnd+1;
cEnd = cStart+cSize-1;


adjMtx = zeros(sum([s.populations.size])+nFreqs);
adjMtx(i2Start:i2End,iStart:iEnd) = i2iNetcon;
adjMtx(iStart:iEnd,rStart:rEnd) = irNetcon;
adjMtx(rStart:rEnd,cStart:cEnd) = rcNetcon;
adjMtx(icStart:icEnd,rStart:rEnd) = icrNetcon;
adjMtx(icStart:icEnd,iStart:iEnd) = iciNetcon;
%% graph
g = digraph(adjMtx,names);
figure; p = plot(g,'layout','layered','Marker','o',...
                   'MarkerSize',8,'NodeLabel',names);
set(p,'ArrowSize',12)
layout(p,'layered','direction','down',...
         'sources',[cStart:cEnd],'sinks',[icStart:icEnd]);
