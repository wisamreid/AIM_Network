% compile data & calculate discriminability
fileLoc = 'C:\Users\Kenny\Desktop\GitHub\SpatialAttentionNetwork\data';
resultFiles = dir([fileLoc filesep '*tmr*.mat'])
for i = 1:length(resultFiles)
    data(i) = load([fileLoc filesep resultFiles(i).name])
end

%% calculate discriminability
h1 = figure(1);
for i = 1:length(data)
results = data(i).results;
variedRange = unique([results.Iapp]);
    for Iapp = 1:length(variedRange)
        % STOI-based calculation
        currentVar(Iapp).stois = reshape([results([results.Iapp]==variedRange(Iapp)).stoi],3,[]);
        [m,idx]=max(currentVar(Iapp).stois);
        disc(Iapp) = sum(idx==1)./numel(idx);
        meanSTOI(Iapp,:) = mean(currentVar(Iapp).stois,2);
    end
set(0, 'CurrentFigure', h1)
plot(variedRange,disc); hold on;
xlabel('I_{app}')
ylabel('Discriminability')
end
tmrs = [0,5,2,-2];
a = cellstr(num2str(tmrs'));
% b = {'TMR'};
% [Bx,Ax] = ndgrid(1:numel(b),1:numel(a));
% legendLabels = strcat(b(Bx(:)),':',a(Ax(:)))
hl = legend(a)
hl.Title.String = 'TMR (dB)';