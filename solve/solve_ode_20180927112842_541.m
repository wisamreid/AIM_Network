function [T,E_V,I_V,E_V_spikes,I_V_spikes]=solve_ode
% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);
% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng_wrapper(p.random_seed);
t=0; k=1;

% STATE_VARIABLES:
E_V = zeros(nsamp,p.E_Npop);
  E_V(1,:) = -65 * ones(1,p.E_Npop);
I_V = zeros(nsamp,p.I_Npop);
  I_V(1,:) = -65 * ones(1,p.I_Npop);

% MONITORS:
E_tspike = -1e6*ones(2,p.E_Npop);
E_buffer_index = ones(1,p.E_Npop);
E_V_spikes = zeros(nsamp,p.E_Npop);
I_tspike = -1e6*ones(2,p.I_Npop);
I_buffer_index = ones(1,p.I_Npop);
I_V_spikes = zeros(nsamp,p.I_Npop);
% ###########################################################
% Numerical integration:
% ###########################################################
% seed the random number generator
rng_wrapper(p.random_seed);
n=2;
for k=2:ntime
  t=T(k-1);
  E_V_k1=(p.E_E-E_V(n-1)+p.E_R*p.E_I+p.E_noise*randn-0)/p.E_tau;
  I_V_k1=(p.I_E-I_V(n-1)+p.I_R*p.I_I+p.I_noise*randn-(((( p.I_E_iampa_gSYN.*sum((( (exp(-(t-E_tspike-p.I_E_iampa_delay)/p.I_E_iampa_tauD)-exp(-(t-E_tspike-p.I_E_iampa_delay)/p.I_E_iampa_tauR)).*((t-E_tspike-p.I_E_iampa_delay)>0)))).*((I_V(n-1))-p.I_E_iampa_ESYN))))))/p.I_tau;
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  E_V(n)=E_V(n-1)+dt*E_V_k1;
  I_V(n)=I_V(n-1)+dt*I_V_k1;
  % ------------------------------------------------------------
  % Conditional actions:
  % ------------------------------------------------------------
  conditional_test=(I_V(n,:)>=p.I_thresh&I_V(n-1,:)<p.I_thresh);
  if any(conditional_test), I_V_spikes(n,conditional_test)=1;inds=find(conditional_test); for j=1:length(inds), i=inds(j); I_tspike(I_buffer_index(i),i)=t; I_buffer_index(i)=mod(-1+(I_buffer_index(i)+1),2)+1; end; end
  conditional_test=(E_V(n,:)>=p.E_thresh&E_V(n-1,:)<p.E_thresh);
  if any(conditional_test), E_V_spikes(n,conditional_test)=1;inds=find(conditional_test); for j=1:length(inds), i=inds(j); E_tspike(E_buffer_index(i),i)=t; E_buffer_index(i)=mod(-1+(E_buffer_index(i)+1),2)+1; end; end
  conditional_test=( any(t<E_tspike + p.E_tref, 1) );
  if any(conditional_test), E_V(n,conditional_test)=p.E_reset; end
  conditional_test=( any(t<I_tspike + p.I_tref, 1) );
  if any(conditional_test), I_V(n,conditional_test)=p.I_reset; end
  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  n=n+1;
end

T=T(1:downsample_factor:ntime);

end

