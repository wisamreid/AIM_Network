function I2input = loadI2input()

disp(pwd)
if exist('i2input.mat','file')
    fileData = load('i2input.mat');
else
    fileData = load('..\i2input.mat');
end
I2input = fileData.i2InputCurrent;