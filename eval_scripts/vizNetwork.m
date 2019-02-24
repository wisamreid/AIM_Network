% vizNetwork
% pops = {data.model.specification.populations.name};
% sz = {data.model.specification.populations.size};

% construct cell array of names
clear names
for l = 1:5
    for f = 1:36
        names{f+(l-1)*36} = sprintf('I2_%i_%i',l,f);
    end
end
for l = 1:5
    for f = 1:36
        itemp{f+(l-1)*36} = sprintf('I_%i_%i',l,f);
    end
end
names = [names itemp];
for l = 1:5
    for f = 1:36
        rtemp{f+(l-1)*36} = sprintf('R_%i_%i',l,f);
    end
end
names = [names rtemp];
for l = 1
    for f = 1:36
        ctemp{f+(l-1)*36} = sprintf('C_%i_%i',l,f);
    end
end
names = [names ctemp];

%% combine connectivity matrices
adjMtx = zeros(180+180+180+36);
adjMtx(1:180,181:180*2) = i2iNetcon;
adjMtx(181:180*2,180*2+1:180*3) = irNetcon;
adjMtx(180*2+1:180*3,180*3+1:end) = rcNetcon;

%%
g = digraph(adjMtx,names);
plot(g);