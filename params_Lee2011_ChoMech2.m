% Lee & Middlebrooks 2011 parameters
% Models the cholinergic mechanism of attention
% No inhibitory connections
% modify off-target EC strength
% model works but we can try to match the parameters form inhibitionMech

varies = struct();
varies(1).conxn = 'IC->IC';
varies(end).param = 'trial'; 
varies(end).range = 1;
% trials -> source location index

varies(end+1).conxn = 'IC->E';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN'; 
varies(end).range = 0.5;

varies(end+1).conxn = 'XI->E';
varies(end).param = 'gSYN'; 
varies(end).range = 1;

varies(end+1).conxn = 'E->C';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'E->XI';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'XI';
varies(end).param = 'Itonic'; 
varies(end).range = 2;

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic'; 
varies(end).range = 1;

% varies(end+1).conxn = 'TD';
% varies(end).param = 'Itonic'; 
% varies(end).range = 4:0.5:8;

% off-target E->C strength: uniform
if attend
    offTargetECGain = 0;
else
    offTargetECGain = 1;
end

% IC->E strength = non-uniform
weights = (min(sum(newTrainReduced))./sum(newTrainReduced)).^2;
weights(9) = 1.1;
weights(10) = 1.5;
weights(11) = 1.5;
weights(12) = 1.5;
weights = weights/max(weights)*1.5;
ICEsmall = diag(weights);
% ICEsmall = eye(19);
options.ICEnetcon = ICEsmall;

% TD neuron current activation - always on.
nLocs = length(azListReduced);
[~,tdChanIdx] = min(abs(azListReduced - attendAz));
TD_Imask = ones(nLocs,1);
options.TD_Imask = TD_Imask;