function networkParams = V26_set_parameters(beamLoc,attendLoc,I2_inh)
% Modified by Kenny Chou
% 2017/8/1

%% Network parameters: Conductance
networkParams.g_direct = [1 1 1 1 1]*.07;
networkParams.g_preinter = 0.11; 
networkParams.g_postinter = 0;
networkParams.g_cortex0 = .07; 
networkParams.g_cortex = [1 1 1 1 1]*networkParams.g_cortex0; 

if ~exist('attendLoc','var')
    disp('missing attend location variable, assume spacing out')
end
    
%5x5, for the relay neurons in each spatial channel (row),
%inhibition is coming from this spatial channel (col)
networkParams.g_latinh=zeros(5);
if beamLoc ~= 0,
    networkParams.g_latinh(:,beamLoc) = 0.2;
end
networkParams.g_I2inh = networkParams.g_latinh;

% %Alternatively, a more complex network can be set up this way
% switch beamLoc
%     case 1 %inhibition from -90 deg, to all locations
%         networkParams.g_latinh(:,1) = 0.2;
%     case 3 %inhibition from 0 deg, to all locations
%         networkParams.g_latinh(:,3) = 0.2;
%     case 5 %inhibition from 90 deg, to all locations
%         networkParams.g_latinh(:,5) = 0.2;
%     otherwise
%         disp('location not supported')
% end

% Inhibitory current from I2 to I neurons
networkParams.I_I2 = zeros(1,5);
networkParams.I_I2(attendLoc) = I2_inh;

    
%% cell properties: Timing
networkParams.t_exc_rise = 1;
networkParams.t_exc_fall = 3;
networkParams.t_inh_rise = 1; %4
networkParams.t_inh_fall = 1000;

% common
networkParams.ref_max = 3; %ms
% interneuron
networkParams.t_preinter = 1;
%% adaptation
% cortical
networkParams.G_inc_l1=0;%0.5/40=.0125
networkParams.tau_ad1=400;

% relay
networkParams.G_inc_l2=0;
networkParams.tau_ad2=400;

% common
%% synaptic depression in cortical (theneuron)
networkParams.tau_sd1=80;
networkParams.d1=.95;
networkParams.d1=1; %why are there repeated numbers here?
networkParams.tau_sd2=1000;
networkParams.d2=.995;
networkParams.d2=1;
networkParams.tau_f=15;
networkParams.f=1.05;
networkParams.f=1;