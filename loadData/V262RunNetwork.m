function [spk_network,spk_relay_out] = V262RunNetwork(spk_subcortical,fs,neuronParams,plot_on,varargin)
% input:
%   spk_subcortical = spikes from the IC model, time x nSpatial x nFreqChan
%   fs = sampling freq
%   nCortical = number of cortical neurons
%   plot_on = plots of spikes for each frequency and spatial channel
%   saveLoc = location for plots. Must be present if plot_on = 1
% Modified by Kenny Chou
% V26 - added cross-freq inhibitory connections by variable frequency
%       tuning in the input layer
% 12/25/2017 on a grounded airplane
%
% 2018/04/06 V26.2 revert back to parallel structure

% initialize
[n_spk, locs, nf] = size(spk_subcortical);
% n_spk: number of spikes
% locs: number locations
% nf: number of frequency channels
weights = gen_synaptic_weights(neuronParams.q,neuronParams.cf);

spk_network = zeros(n_spk,neuronParams.nCortical,nf);  
spk_relay_out = zeros(n_spk,locs,nf);
% V_cortical = zeros(n_spk,5,nf);
dt = 1000/fs;
t_end = n_spk*1000/fs;

if plot_on
    saveLoc = varargin{1};
    savePlotLoc = sprintf('%s_%s rise%.2f fall%d',saveLoc,'plots',neuronParams.t_inh_rise,neuronParams.t_inh_fall);
    if ~exist(savePlotLoc,'dir'), mkdir(savePlotLoc); end
else
    savePlotLoc = '';
end

h = dt; % step size, Euler method, = dt ms
tstop=round(t_end/h);% number of time steps
t=0:h:h*(tstop-1);
thr = 0; % threshold for determining whether pre-synaptic input has spikes, NOT the same as V_th
g_preinter = neuronParams.g_preinter;
t_preinter = neuronParams.t_preinter;
g_postinter = neuronParams.g_postinter;
g_relay = neuronParams.g_direct; %<------ synaptic weight dependent
g_latinh = neuronParams.g_latinh;
g_I2inh = neuronParams.g_I2inh;
t_exc_rise = neuronParams.t_exc_rise;
t_exc_fall = neuronParams.t_exc_fall;
t_inh_rise = neuronParams.t_inh_rise;
t_inh_fall = neuronParams.t_inh_fall;
ref_max = neuronParams.ref_max;
G_inc_l2 = neuronParams.G_inc_l2;
tau_ad2 = neuronParams.tau_ad2;
I_I2 = neuronParams.I_I2;

% Convert time constants from msec to time steps
ref_max = ref_max/dt;

%% Run CISPA Network - execute each network layer consecutively
% interneurons
length_inter=length(t);
spike_inter=zeros(length_inter,locs,nf);
V_inter = zeros(length_inter,locs,nf);
I_syn_exc = zeros(length_inter,locs,nf);
I_syn_inh = zeros(length_inter,locs,nf);
freqs = 1:nf; %for debugging purposes
parfor freq = freqs
    for i=1:locs
%         V_inter(:,i,freq) = V19_NEURON_inter(spk_subcortical(:,i,freq),t,g_preinter,t_preinter,ref_max);
        [V_inter(:,i,freq),I_syn_inh(:,i,freq),g_syn_inh,I_syn_exc(:,i,freq),g_relay_exc]=V25_neuron_inter(I_I2(i),spk_subcortical(:,i,freq),t,...
                g_preinter,g_I2inh,t_inh_rise,t_inh_fall,t_exc_rise,t_exc_fall,ref_max,tau_ad2,G_inc_l2,i);
    end

%   debugging purposes
%     figure;
%     freq = 25;
%     for i = 1:locs
%         subplot(3,5,i); plot(V_inter(:,i,freq)); title('voltage out');
%         n = sum(V_inter(:,i,freq)>thr);
%         title(sprintf('v_{out}, %d spikes',n))
%         ylim([-70,75])
%         subplot(3,5,i+5); plot(I_syn_exc(:,i,freq)); title('I_{exc}');
%         subplot(3,5,i+10); plot(I_syn_exc(:,i,freq)+I_syn_inh(:,i,freq)); title('I_{tot}');
%     end
%     cellstr(num2str((1:5)'))
end

spike_inter(V_inter > thr) = 1;
spk_inh = spike_inter;
 
%relay neurons
length_relay=length(spike_inter);
parfor freq = freqs
    freqSpreadIdx = find(weights(freq).weights>0.0001);
    freqWeights = weights(freq).weights(freqSpreadIdx);
    spike_relay=zeros(length_relay,locs);
    V_relay=zeros(length_relay,locs);
    I_inh=zeros(length_relay,locs);
    I_relay=zeros(length_relay,locs);
    g_relay_hist=zeros(length_relay,locs);
    g_inh_hist=zeros(length_relay,locs);
    for i=1:locs
        [V_relay(:,i),I_inh(:,i),I_relay(:,i),g_relay_hist(:,i),g_inh_hist(:,i)] = V26_NEURON_relay(spk_inh(:,:,freq),spk_subcortical(:,:,freqSpreadIdx),freqWeights,t,...
            g_postinter,g_relay(i),g_latinh(locs+1-i,:),t_inh_rise,t_inh_fall,t_exc_rise,t_exc_fall,ref_max,tau_ad2,G_inc_l2,i);    
    end
    spike_relay(V_relay > thr) = 1;
    spk_relay_out(:,:,freq) = spike_relay; 
    
% figure; plot(I_inh); title('I_{inh relay}');
% figure;plot(g_inh_hist); title('conductance inh hist');
% figure;plot(V_relay); title('V_{relay}');
% cellstr(num2str((1:5)'))

%cortical neurons
    V_cortical = zeros(size(V_relay));
    r = -45:2:-35; %V_th between -35 and -45 mV
    V_th = ones(1,neuronParams.nCortical);
    for i = 1:neuronParams.nCortical %n cortical neurons
        idx = mod(i,length(r));
        V_th(i) = r(idx+1);
        V_cortical(:,i) = V22_Neuron_Cortical([zeros((i-1)*200,5);spike_relay(1:end-(i-1)*200,:)],dt,tstop,locs,neuronParams,V_th(i));
    end
    temp = zeros(n_spk,neuronParams.nCortical);
    temp(V_cortical > 0) = 1;
    spk_network(:,:,freq) = temp;
 
    
    %plots
    if plot_on
        plotNetworkActivities(n_spk,fs,locs,freq,V_cortical,V_relay,V_inter(:,:,freq),spk_subcortical(:,:,freq),savePlotLoc);
    end
end