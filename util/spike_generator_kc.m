function [spike_train,spike_times]=spike_generator_kc(spike_rate,dt,frMax)
% poisson spiking generator
% Inputs:
%   [spike_rate]
%   [time]          time vector for spike rate (unit: s)
%   [frMax]         maximum spiking rate
% Outputs:
%   [spike_train]        array with 1's for all spike times
%   [spike_times]        unit in s

if ~exist('frMax','var') frMax = 500; end

%% Spike generation
n=length(spike_rate);
tw=0:dt:6/1000; % (ms) time vector for recovery rate
w=1-1./(1+exp(4000*(tw-3.8/1000))); %weighting vector; not used
n_refrac=6/1000/dt; %refractory period 6ms

% If random, have random spiking. Otherwise, always spike
spike_train=zeros(1,n);
spike_times=[];

x = frMax;
spike_rate = spike_rate/150*x;
i = 1;
while i < n
    if mod(i,1) ~= 0
        i = round(i);
    end
    
    if spike_rate(i)*rand > 50
        spike_train(i) = 1;
        spike_times = [spike_times; i];
        i = i + min(round(1/spike_rate(i)/dt),n_refrac);
    else
        i = i+1;
    end
end

spike_times=spike_times*dt;
