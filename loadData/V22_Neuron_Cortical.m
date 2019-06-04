function V_trace = V22_Neuron_Cortical(spike_relay,dt,tstop,locs,corticalParams,E_leak)

% Input Parameters
tau_sd1 = corticalParams.tau_sd1;
tau_sd2 = corticalParams.tau_sd2;
tau_f = corticalParams.tau_f;
d1 = corticalParams.d1;
d2 = corticalParams.d2;
f = corticalParams.f;
tau_ad1 = corticalParams.tau_ad1;
G_inc_l1 = corticalParams.G_inc_l1;
g_cortex = corticalParams.g_cortex;


% Timing Parameters
ref = 0; % refractory period counter
ref_max = corticalParams.ref_max/dt;

% Voltage diffeq
% alpha func synaptic conductance
t_a = 100; % Max duration of syn conductance=100ms
t_vec = 0:dt:t_a;
tau2 = corticalParams.t_exc_rise;
tau1 = corticalParams.t_exc_fall;
tau_rise = tau1*tau2/(tau1-tau2);
b = ([(tau2/tau1)^(tau_rise/tau1)-(tau2/tau1)^(tau_rise/tau2)]^-1)/tau2;
epsc_cortex = b*[exp(-t_vec/tau1)-exp(-t_vec/tau2)];

% capacitance and leak resistance
C = 0.5; % nF
R = 40; % M ohms
g_ad = 0; % uS
D0 = 1; % maximum proportion of g_peak
D1 = zeros(tstop,locs);
D1(1,:) = D0*ones(1,locs); % fraction of g_peak remaining (initialized to 1)
D2  = zeros(tstop,locs);
D2(1,:) = D0*ones(1,locs); % fraction of g_peak remaining (initialized to 1)
F = zeros(tstop,locs);
F(1,:) = D0*ones(1,locs); % fraction of g_peak remaining (initialized to 1)

% Initialize basic parameters
% E_leak = -60; % mV, equilibrium potential
E_k = -70;
E_syn_exc = 0;
V_reset = -80; %reset voltage after spike, mV
V_th = -40; % spike threshold mV
V_spike = 50; % spike value mV
t_list_cortical = cell(1,locs);
V = E_leak;
V_trace = V*ones(tstop,1);
t_trace = 0:dt:(tstop-1)*dt';
g_cortex_real = zeros(locs,1);
I_syn_cortex = zeros(locs,tstop);
g_syn_cortical = zeros(locs,tstop);
I_cortical = zeros(locs,1);

%%
for t = 1:tstop
    %% check for inputs from relay neurons
    for loc = 1:locs
        if (spike_relay(t,loc) > 0) % check for input spike
            t_list_cortical{loc} = [t_list_cortical{loc}; 1];
            g_cortex_real(loc) = g_cortex(loc)*D1(t,loc)*D2(t,loc)*F(t,loc);
            D1(t+1,loc) = D1(t,loc)*d1;
            D2(t+1,loc) = D2(t,loc)*d2;
            F(t+1,loc) = F(t,loc)*f;
        else
            D1(t+1,loc) = D1(t,loc) + dt*(1/tau_sd1)*(D0-D1(t,loc)); % pre-syn recovery
            D2(t+1,loc) = D2(t,loc) + dt*(1/tau_sd2)*(D0-D2(t,loc)); % pre-syn recovery
            F(t+1,loc) = F(t,loc) + dt*(1/tau_f)*(D0-F(t,loc)); % pre-syn recovery
        end
        g_syn_cortical(loc,t) = sum(g_cortex_real(loc)*epsc_cortex(t_list_cortical{loc}));
        I_cortical(loc) = g_syn_cortical(loc,t)*(E_syn_exc - V);
        I_syn_cortex(loc,t) = I_cortical(loc);

        if t_list_cortical{loc}
            t_list_cortical{loc} = t_list_cortical{loc} + 1;
            if (t_list_cortical{loc}(1) == length(t_vec))  % Reached max duration of syn conductance
                t_list_cortical{loc} = t_list_cortical{loc}(2:max(size(t_list_cortical{loc})));
            end
        end
    end
    
    %% Compute membrane voltage
    if ~ref
        % Euler method: V(t+h) = V(t) + h*dV/dt
        V = V  + sum(I_cortical)/C*dt + (E_leak-V)/(R*C)*dt + (E_k-V)*g_ad/C*dt; 
        % Dayan analytical solution
%         V = E_leak + R_m*I_e + (V - E_L - R_m*I_e )*exp(-dt/tau_m);
    elseif ref==ref_max % If AP just happened, set V to V_th-10
        ref = ref - 1;
%         V = V_th - 10;
        V = V_reset; % reset voltage after spike
    else
        ref = ref - 1;
    end
    g_ad = g_ad + dt*(- g_ad/tau_ad1); % spike rate adaptation
    
    % Generate spike
    if V > V_th
        V = V_spike;
        g_ad = g_ad + G_inc_l1;
        ref=ref_max;
    end
    
    V_trace(t) = V;
end
