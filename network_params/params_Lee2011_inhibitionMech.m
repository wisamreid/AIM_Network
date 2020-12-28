% Lee & Middlebrooks 2011 parameters
% same model as Atiani et al., and Fritz et al.

varies = struct();
varies(1).conxn = 'IC->IC';
varies(end).param = 'trial'; 
varies(end).range = 1;
% trials -> source location index

varies(end+1).conxn = 'IC->E';
varies(end).param = 'gSYN'; 
varies(end).range = 2.5;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN'; 
varies(end).range = 2.4;

varies(end+1).conxn = 'XI->E';
varies(end).param = 'gSYN'; 
varies(end).range = 3;

varies(end+1).conxn = 'E->C';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'XI';
varies(end).param = 'Itonic'; 
varies(end).range = 3;

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic'; 
varies(end).range = 0;

% off-target E->C strength: uniform
offTargetECGain = 1;

% IC->E strength = uniform
% weights = (min(sum(newTrainReduced))./sum(newTrainReduced)).^2;
% weights(9) = 1.1;
% weights(10) = 1.5;
% weights(11) = 1.5;
% weights(12) = 1.5;
% weights = weights/max(weights)*1.5;
% ICEsmall = diag(weights);
ICEsmall = eye(19);
options.ICEnetcon = ICEsmall;

% TD neuron current activation - turn off tgt TD neuron if attending.
nLocs = length(azListReduced);
[~,tdChanIdx] = min(abs(azListReduced - attendAz));
TD_Imask = ones(nLocs,1);
if attend, TD_Imask(tdChanIdx) = 0; end
options.TD_Imask = TD_Imask;