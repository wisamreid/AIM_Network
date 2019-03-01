% vizNetwork
% pops = {data.model.specification.populations.name};
% sz = {data.model.specification.populations.size};

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
names = [i2Names iNames rNames cNames]; %order matters!

%% combine connectivity matrices, in the order of I2,I,R,C
i2Start = 1;
i2End = i2Size;
iStart = i2Size+1;
iEnd = iStart+iSize-1;
rStart = iEnd+1;
rEnd = rStart+rSize-1;
cStart = rEnd+1;
cEnd = cStart+cSize-1;

adjMtx = zeros(sum([s.populations.size]));
adjMtx(i2Start:i2End,iStart:iEnd) = i2iNetcon;
adjMtx(iStart:iEnd,rStart:rEnd) = irNetcon;
adjMtx(rStart:rEnd,cStart:cEnd) = rcNetcon;

%% graph
g = digraph(adjMtx,names);
figure; p = plot(g,'layout','layered','Marker','o',...
                   'MarkerSize',8,'NodeLabel',names);
set(p,'ArrowSize',12)
layout(p,'layered','direction','down',...
         'sources',[i2Start:i2End],'sinks',[cStart:cEnd]);