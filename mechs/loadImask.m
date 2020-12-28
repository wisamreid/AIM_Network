function Imask = loadImask()
ImaskFile = 'Imask.mat';
if exist(ImaskFile,'file')
    load(ImaskFile,'Imask');
else
    load(['..' filesep ImaskFile],'Imask');
end