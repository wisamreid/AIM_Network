function [V_trace,I_syn_inh,I_syn_relay,g_relay_hist,g_inh_hist]=V26_NEURON_relay(spike_postinter,spike_train_relay,qWeights,t,...
    g_postinter,g_relay,g_latinh,t_inh_rise,t_inh_fall,t_exc_rise,t_exc_fall,ref_max,tau_ad,G_inc,neuron_id)
% V26: enables neuron to accept input from multiple IC neurons (from
% different CFs). Inputs correcponding to other CFs would have a different
% g_relay weighting. Edits are made under the "direct inputs" section.

% d=fraction of d consumed by each spike = free parameter
n_spatial_channels = length(g_latinh); % 120916 edit

dt=t(2)-t(1);
h = dt; % step size, Euler method, = dt ms
tstop = length(t); % number of time steps
ref = 0; % refractory period counter
t_a = 1000; % Max duration of syn conductance (ms)
% t = 0:h:t_a;
t_vec = 0:h:t_a;
ncf = length(qWeights);

% EPSC and IPSC
tau2=t_exc_rise;tau1=t_exc_fall;tau_rise=tau1*tau2/(tau1-tau2);
b=([(tau2/tau1)^(tau_rise/tau1)-(tau2/tau1)^(tau_rise/tau2)]^-1)/tau2;
epsc_relay=b*[exp(-t_vec/tau1)-exp(-t_vec/tau2)];

tau2=t_inh_rise;tau1=t_inh_fall;tau_rise=tau1*tau2/(tau1-tau2);
b=([(tau2/tau1)^(tau_rise/tau1)-(tau2/tau1)^(tau_rise/tau2)]^-1)/tau2;
ipsc_latinh=b*[exp(-t_vec/tau1)-exp(-t_vec/tau2)];

% capacitance and leak resistance
C = 0.5; % nF
R = 40; % M ohms
g_ad = 0; % uS
D0 = 1; % maximum proportion of g_peak 
D = D0; % fraction of g_peak remaining (initialized to 1)

% 
E_leak = -60; % mV, equilibrium potential
E_k=-70; % adaptation channel
E_exc = 0;
E_inh = -70; % inhibition channel
V_th = -40; % spike threshold mV
V_spike = 50; % spike value mV

% Initialize 
V = E_leak;
t_list_inter = cell(1,n_spatial_channels); % possible inhibitory input from 4 inh channels
% t_list_inter = num2cell(ones(1,n_spatial_channels)*5000);
t_list_relay = cell(1,ncf); % possible inputs from ncf neighboring frequency channels

V_trace = V*ones(tstop,1);
% t_trace = 0:dt:(tstop-1)*dt';
% t_list_inter = cell(1,4); % possible inhibitory input from 4 inh channels
% t_list_relay = [];
I_syn_relay=zeros(tstop,1);
I_syn_inh=zeros(tstop,1);
g_relay_hist=zeros(tstop,1);
g_inh_hist=zeros(tstop,1);
g_syn_inh=zeros(n_spatial_channels,1);
g_syn_relay=zeros(ncf,1);
for t = 1:tstop
    %% Direct relay inputs
    for nn = 1:ncf
        if (spike_train_relay(t,neuron_id,nn) > 0)
            t_list_relay{nn} = [t_list_relay{nn}; 1];
        end
        g_syn_relay(nn) = sum(g_relay*qWeights(nn)*epsc_relay(t_list_relay{nn}));
    end
    I_exc = sum(g_syn_relay)*(E_exc - V);

    I_syn_relay(t)=I_exc;
    g_relay_hist(t) = sum(g_syn_relay);

    %% Check for inhibitory inputs (direct and lateral) from interneurons
    for i = 1:n_spatial_channels
        if (spike_postinter(t,i) > 0) % check for input spike
            t_list_inter{i} = [t_list_inter{i}; 1];
        end
        if i==neuron_id
            g_syn_inh(i,1) = sum(g_postinter*ipsc_latinh(t_list_inter{i}));
        else
            g_syn_inh(i,1) = sum(g_latinh(i)*ipsc_latinh(t_list_inter{i}));
        end
    end
    
    g_syn_inh_sum=sum(g_syn_inh);
    I_inh = g_syn_inh_sum*(E_inh - V);
    if neuron_id ~= find(g_latinh), I_inh = I_inh-0.24; end
    
    g_inh_hist(t)=g_syn_inh_sum;
    I_syn_inh(t)=I_inh;
    %% Sum all synaptic inputs and integrate
    I_syn=I_exc+I_inh;
    
    %this block is here for debugging purposes
%     if t >= 61060
%         [V, I_inh, I_exc]
%         delta_v = I_syn/C + (E_leak-V)/(R*C)
%         disp('hello')
%     end
    
    % Compute membrane voltage
    % Euler method: V(t+h) = V(t) + h*dV/dt
    if ~ref
        V = V + h*I_syn/C + (E_leak-V)/(R*C)*h+(E_k-V)*g_ad/C*h;
    elseif ref==ref_max % If AP just happened, set V to V_th-10
        ref = ref - 1;
        V = V_th-10; % reset voltage after spike
    else
        ref = ref - 1;
    end
    g_ad = g_ad + h*(- g_ad/tau_ad); % spike rate adaptation
    
    % Generate spike
    if V > V_th
        V = V_spike;
        g_ad = g_ad + G_inc;
        ref=ref_max;
    end
    
    V_trace(t)= V;
    
    %increment time
    t_list_relay = incrementTime(t_list_relay,length(t_vec),ncf);
    t_list_inter = incrementTime(t_list_inter,length(t_vec),n_spatial_channels);
    
end
end

function out = incrementTime(tList,max_time,n)
    for i=1:n
        if tList{i}
            tList{i} = tList{i} + 1;
            if (tList{i}(1) == max_time)  % Reached max duration of syn conductance
                tList{i} = tList{i}(2:max(size(tList{i})));
            end
        end
    end
    out = tList;
end