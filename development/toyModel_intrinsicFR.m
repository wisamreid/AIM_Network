% 3 neuron network

dynasimPath = 'C:\Users\Kenny\Desktop\GitHub\DynaSim';
addpath(genpath(dynasimPath));
addpath('mechs')

% parameters
% Neuron populations
nPops = 1;
itonic = 8; % defualt is 1,2, and 3;
gsyn = 3; % default is 3;

% vary parameters
varies = struct();
varies(1).conxn = 'E'; 
varies(end).param = 'Itonic';
varies(end).range = 3;

varies(end+1).conxn = 'I'; 
varies(end).param = 'Itonic';
varies(end).range = 6;

varies(end+1).conxn = 'TD'; 
varies(end).param = 'Itonic';
varies(end).range = 7;

varies(end+1).conxn = 'TD'; 
varies(end).param = 'noise';
varies(end).range = 0;

varies(end+1).conxn = 'TD->I';
varies(end).param = 'gSYN';
varies(end).range = 3;

varies(end+1).conxn = 'TD->E';
varies(end).param = 'gSYN';
varies(end).range = 3;

varies(end+1).conxn = 'E';
varies(end).param = 'Itonic';
varies(end).range = 9:16;

varies(end+1).conxn = 'E';
varies(end).param = 'noise';
varies(end).range = 0;

varied_param = find(cellfun(@length,{varies.range})>1);
expVar = [varies(varied_param).conxn ' ' varies(varied_param).param];

vary = cell(length(varies),3);
for i = 1:length(varies)
    vary{i,1} = varies(i).conxn;
    vary{i,2} = varies(i).param;
    vary{i,3} = varies(i).range;
end

% netcons
ABnetcon = eye(nPops)*1;
% ABnetcon(2,2) = 2;
% ABnetcon(3,3) = 3;

t_end = 50; %ms

s=[];
s.pops(1).name='TD';
s.pops(1).size=nPops;
s.pops(1).equations='chouLIF';
s.pops(1).parameters={'Itonic',itonic};

s.pops(2).name='E';
s.pops(2).size=nPops;
s.pops(2).equations='chouLIF';
s.pops(2).parameters={'toff',t_end};

s.pops(end+1).name='I';
s.pops(end).size=nPops;
s.pops(end).equations='chouLIF';
s.pops(end).parameters={'Itonic',itonic};

s.cons(1).direction='TD->I';
s.cons(end).mechanism_list='dsSynapse';
s.cons(end).parameters={'delay',0,'gSYN',gsyn, 'tauR',1, 'tauD',10,'ESYN',-80};

s.cons(end+1).direction='TD->E';
s.cons(end).mechanism_list='dsSynapse';
s.cons(end).parameters={'delay',0,'gSYN',gsyn, 'tauR',1, 'tauD',10,'ESYN',-80}; 

% simulation
data=dsSimulate(s,'vary',vary,'time_limits',[0 t_end],'solver','rk1','dt',.01);
% dsPlot(data,'plot_type','rastergram');

%% plot rasters

xstart = 0;
xend = 67.0;
imgHeight = 0.8/length(data);
pops = {data(1).model.specification.populations.name};
fieldNames = strcat(pops,'_V_spikes');
for j = 1:length(pops)
    figure
    for i = 1:length(data)
        ypos = 0.9-imgHeight*(i);
        subplot('position',[0.1 ypos 0.8 imgHeight])
        spikes = [data(i).(fieldNames{j})];
        plot(data(i).time, spikes);
        fr = sum(spikes)/(t_end/1000);
        ylabel(varies(varied_param).range(i))
        yticks([])
        xticks([])
        xlim([xstart xend]);
        text(52.0,0.5,['FR: ' num2str(round(fr,2))])
    end
    suptitle([pops{j} ' spikes; vary' expVar])
end
return;
%%
if isempty(varied_param) == 1
    yspace = 0.95/5;
    height = yspace*0.75;
    width = 0.95/2*0.9;
    xpadSingle = (1-1*width)/2;
    xpadDouble = (1-2*width)/2;

    
    subData = data([data.TD_Itonic]==8);
    for i = 1:length(subData)
        figure;
        subplot('position',[xpadSingle 0.05+4*yspace width height]);
        plot(subData(i).time, subData(i).TD_V(:,1)); title('TD voltage - current pulse input')
        xlim([xstart xend]);
        subplot('position',[xpadSingle 0.05+3*yspace width height]); 
        plot(subData(i).time, subData(i).E_TD_dsSynapse_Isyn(:,1)); title('synaptic current: dsSynapse')
        xlim([xstart xend]);
        subplot('position',[xpadSingle 0.05+2*yspace width height]); 
        plot(subData(i).time, subData(i).E_V(:,1)); title('E voltage')
        xlim([xstart xend]);
        subplot('position',[xpadSingle 0.05+1*yspace width height]); 
        plot(subData(i).time, subData(i).I_TD_dsSynapse_Isyn(:,1)); title('synaptic current: dsSynapse')
        xlim([xstart xend]);
        subplot('position',[xpadSingle 0.05+0*yspace width height]); 
        plot(subData(i).time, subData(i).I_V(:,1)); title('I voltage')
        xlim([xstart xend]);
    end
end

return;

%% 1d plots
figure;
if length(varied_param) == 1
   for i = 1:length(data)
      ICspks(i) = sum(data(i).IC_V_spikes); 
      Espks(i) = sum(data(i).E_V_spikes);
      TDspks(i) = sum(data(i).TD_V_spikes);
   end
   plot(varies(varied_param).range,ICspks); hold on;
   plot(varies(varied_param).range,Espks);
   plot(varies(varied_param).range,TDspks);
   legend('IC','E','TD')
   xlabel(expVar)
   ylabel('spkCount')
end

%% 2d plots
if length(varied_param) == 2
    clear AFR BFR CFR
    nvar1 = length(varies(1).range);
    nvar2 = length(varies(2).range);
    for i = 1:length(data)
        AFR(i) = sum(data(i).A_V_spikes)/0.1;
        BFR(i) = sum(data(i).B_V_spikes)/0.1;
        CFR(i) = sum(data(i).C_V_spikes)/0.1;
    end
    AFR = reshape(AFR,nvar2,nvar1);
    BFR = reshape(BFR,nvar2,nvar1);
    CFR = reshape(CFR,nvar2,nvar1);

    figure;
    imagesc(varies(1).range,varies(2).range,AFR)
    xlabel('A I tonic')
    ylabel('A->B gSYN')
    title('neuron A FR')

    figure;
%     imagesc(AFR(1,:),CFR(:,1),BFR)
    imagesc(varies(1).range,varies(2).range,BFR)
    xlabel('A I tonic')
    ylabel('C I tonic')
    title('neuron B FR')

    figure;
    cmap = brewermap(50,'RdBu');
    imagesc(AFR(1,:),varies(2).range,AFR-BFR)
    xlabel('A Firing Rate')
    ylabel('C Firing Rate')
    title('Difference in FR (A-B)')
%     caxis([-300 300])
    colormap(cmap)
end

%%
clear Bspks
subData = data([data.A_Itonic]==8);
for i = 1:length(subData)
     Bspks(i) = sum(subData(i).B_V_spikes)/0.1;
end
figure;
plot([subData.C_Itonic],Bspks)
xlabel('C Itonic')
ylabel('B spikes')

% figure;
% plot(varies(varied_param).range,AFR,varies(varied_param).range,BFR)
% xlabel(expVar)
% ylabel('FR (Hz)')
% hl = legend('neuron A','neuron B');
% hl.Location = 'northwest';

%% Definitions, just in case we need it
% % % basic cell definition
% % LIF={
% %     'tau=10; tref=.5; E_leak=-70; V_thresh=-55; V_reset=-75; R=10; noise=0'
% %     'dV/dt=( (E_leak-V) + noise*randn -@isyn + R*I(t) )/tau; V(0)=-70'
% %     'if( any(t<tspike + tref, 1) )(V=V_reset)'
% %     'if(V >= V_thresh)(V = V_reset)'
% %     'g_leak = 1/10; C=1;'
% %     'Itonic = 0;         % injected current, [nA]'
% %     'ton = 100;            % [ms]'
% %     'toff = 900;        % [ms]'
% %     'I(t)=(Itonic+Itonic*square(2*pi*f*t))*(t>ton&t<toff)'
% %     'f=0.005'
% %     'monitor V.spikes(V_thresh,2)'
% %      };
% % 
% % % synapse 1: DS definition
% % iampa={
% %   'gSYN=0.5; ESYN=0; tauD=2; tauR=0.4; delay=0'
% %   'netcon=eye(N_pre,N_post)'
% %   'f(x) = (exp(-x/tauD)-exp(-x/tauR)).*(x>0)'
% %   'Isyn(X) = gSYN.*(sum(f(t-tspike_pre-delay))*netcon).*(X-ESYN)'
% %   '@isyn += Isyn(V_post)'
% %   'monitor Isyn'
% % };