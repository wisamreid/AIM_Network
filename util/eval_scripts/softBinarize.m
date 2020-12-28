function out = softBinarize(in)

sigmoid = @(x,max,k,x0) max./(1+exp(-k.*(x-x0)));
% x = 0:.01:1; %membrane potential range
L = 1; %max firing rate
k = 50;
x0 = (max(max(in))-min(min(in)))*0.2;


L = L*ones(size(in));
k = k*ones(size(in));
x0 = x0*ones(size(in));
out = sigmoid(in,L,k,x0);
figure;imagesc(out)