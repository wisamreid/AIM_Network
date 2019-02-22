function [T,E_pre_V,E_post_V,E_pre_V_spikes,E_post_V_spikes]=solve_ode
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
E_pre_V = zeros(nsamp,p.E_pre_Npop);
  E_pre_V(1,:) = -65 * ones(1,p.E_pre_Npop);
E_post_V = zeros(nsamp,p.E_post_Npop);
  E_post_V(1,:) = -65 * ones(1,p.E_post_Npop);

% MONITORS:
E_tspike = -1e6*ones(2,p.E_Npop);
E_buffer_index = ones(1,p.E_Npop);
E_pre_V_spikes = zeros(nsamp,p.E_pre_Npop);
E_tspike = -1e6*ones(2,p.E_Npop);
E_buffer_index = ones(1,p.E_Npop);
E_post_V_spikes = zeros(nsamp,p.E_post_Npop);
% ###########################################################
% Numerical integration:
% ###########################################################
% seed the random number generator
rng_wrapper(p.random_seed);
n=2;
for k=2:ntime
  t=T(k-1);
  E_pre_V_k1=(p.E_pre_E-E_pre_V(n-1)+p.E_pre_R*p.E_pre_I+p.E_pre_noise*randn-0)/p.E_pre_tau;
  E_post_V_k1=(p.E_post_E-E_post_V(n-1)+p.E_post_R*p.E_post_I+p.E_post_noise*randn-(((( p.E_post_E_pre_iampa_gSYN.*sum((( (exp(-(t-E_pre_tspike-p.E_post_E_pre_iampa_delay)/p.E_post_E_pre_iampa_tauD)-exp(-(t-E_pre_tspike-p.E_post_E_pre_iampa_delay)/p.E_post_E_pre_iampa_tauR)).*((t-E_pre_tspike-p.E_post_E_pre_iampa_delay)>0)))).*((E_post_V(n-1))-p.E_post_E_pre_iampa_ESYN))))))/p.E_post_tau;
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  E_pre_V(n)=E_pre_V(n-1)+dt*E_pre_V_k1;
  E_post_V(n)=E_post_V(n-1)+dt*E_post_V_k1;
  % ------------------------------------------------------------
  % Conditional actions:
  % ------------------------------------------------------------
  conditional_test=(E_post_V(n,:)>=p.E_post_thresh&E_post_V(n-1,:)<p.E_post_thresh);
  if any(conditional_test), E_post_V_spikes(n,conditional_test)=1;inds=find(conditional_test); for j=1:length(inds), i=inds(j); E_tspike(E_buffer_index(i),i)=t; E_buffer_index(i)=mod(-1+(E_buffer_index(i)+1),2)+1; end; end
  conditional_test=(E_pre_V(n,:)>=p.E_pre_thresh&E_pre_V(n-1,:)<p.E_pre_thresh);
  if any(conditional_test), E_pre_V_spikes(n,conditional_test)=1;inds=find(conditional_test); for j=1:length(inds), i=inds(j); E_tspike(E_buffer_index(i),i)=t; E_buffer_index(i)=mod(-1+(E_buffer_index(i)+1),2)+1; end; end
  conditional_test=( any(t<E_pre_tspike + p.E_pre_tref, 1) );
  if any(conditional_test), E_pre_V(n,conditional_test)=p.E_pre_reset; end
  conditional_test=( any(t<E_post_tspike + p.E_post_tref, 1) );
  if any(conditional_test), E_post_V(n,conditional_test)=p.E_post_reset; end
  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  n=n+1;
end

T=T(1:downsample_factor:ntime);

end

