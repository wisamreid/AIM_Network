% Lee & Middlebrooks 2011 parameters
% Models the cholinergic mechanism of attention
% No inhibitory connections
% modify off-target EC strength

varies(end+1).conxn = 'IC->E';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN'; 
varies(end).range = 0;

varies(end+1).conxn = 'XI->E';
varies(end).param = 'gSYN'; 
varies(end).range = 0;

varies(end+1).conxn = 'E->C';
varies(end).param = 'gSYN'; 
varies(end).range = 2;

varies(end+1).conxn = 'XI';
varies(end).param = 'Itonic'; 
varies(end).range = 3.5;

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic'; 
varies(end).range = 0;

if attend == 1
    offTargetECGain = 0.05;
elseif attend == 0
    offTargetECGain = 0.8;
end
ICEsmall = eye(19);
ICEsmall(8,8) = 3;
ICEsmall(9,9) = 7;
ICEsmall(10,10) = 3;
options.ICEnetcon = ICEsmall;