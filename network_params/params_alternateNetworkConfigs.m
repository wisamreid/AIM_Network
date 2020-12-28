% scenario 1: Suppression action of I comes from release of inhibition ¡V TD-I is important
% scenario 2: Feed back inhibition ¡V tonic firing of E neurons ¡V E-I is important
% scenario 3: Direct injection into I due to CHO neuromodulators

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
varies(end).range = 2.8;

varies(end+1).conxn = 'XI->E';
varies(end).param = 'gSYN'; 
varies(end).range = 3;

varies(end+1).conxn = 'E->C';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'E->XI';
varies(end).param = 'gSYN'; 
varies(end).range = 2.5;

% varies(end+1).conxn = 'TD';
% varies(end).param = 'Itonic'; 
% varies(end).range = 0:1:8;

nLocs = length(azListReduced);
[~,tgtChanIdx] = min(abs(azListReduced - attendAz));
TD_Imask = ones(nLocs,1);
offTargetECGain = 1;
nicotinicGain = 1;
if mech==1
    if attend, TD_Imask(tgtChanIdx) = 0; end
    varies(end+1).conxn = 'TD->XI';
    varies(end).param = 'gSYN'; 
    varies(end).range = 4;
    
    varies(end+1).conxn = 'E';
    varies(end).param = 'Itonic'; 
    varies(end).range = 0;
    
    varies(end+1).conxn = 'XI';
    varies(end).param = 'Itonic'; 
    varies(end).range = 3;
elseif mech == 2
    if attend, TD_Imask(tgtChanIdx) = 0; end
    varies(end+1).conxn = 'TD->XI';
    varies(end).param = 'gSYN'; 
    varies(end).range = 0;
    
    varies(end+1).conxn = 'E';
    varies(end).param = 'Itonic'; 
    varies(end).range = 2.5;
    
    varies(end+1).conxn = 'XI';
    varies(end).param = 'Itonic'; 
    varies(end).range = 0;
elseif mech == 3
    if attend
        TD_Imask(tgtChanIdx) = 0;
        varies(end+1).conxn = 'XI';
        varies(end).param = 'Itonic'; 
        varies(end).range = 3;
    else
        varies(end+1).conxn = 'XI';
        varies(end).param = 'Itonic'; 
        varies(end).range = 0;
    end
    varies(end+1).conxn = 'TD->XI';
    varies(end).param = 'gSYN'; 
    varies(end).range = 0;
    
    varies(end+1).conxn = 'E';
    varies(end).param = 'Itonic'; 
    varies(end).range = 0;
end
options.TD_Imask = TD_Imask;


% IC->E strength = non-uniform
% weights = (min(sum(newTrainReduced))./sum(newTrainReduced)).^1;
weights = ones(19);
weights(1:5) = weights(1:5)*0.9;
weights(9) = weights(9)*1.7;
weights(10) = weights(10)*1.7;
weights(11) = weights(11)*1.7;
weights(13:19) = weights(13:19)*0.6;
weights = weights/max(weights)*nicotinicGain;
ICEsmall = diag(weights);
% ICEsmall = eye(19);
options.ICEnetcon = ICEsmall;

% adjust off-target E-C strength
nLocs = numel(azListReduced);
ECsmall = ones(nLocs,1)*offTargetECGain;
ECsmall(tgtChanIdx) = 1.2;
options.ECnetcon = ECsmall;
