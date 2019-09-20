function cc = xcorrKC(a,b)
% Normalized cross correlation for small lags in the x-direction
% a and b can be vectors or 2d matrices (have not tested higher dim)
%
% a = rstim.tf.L+rstim.tf.R;
% b = ref(jj).tf;

lags = 0;
abig = zeros(size(a,1),size(a,2)+length(lags)-1);
abig(:,abs(min(lags))+1:abs(min(lags))+size(a,2)) = a;
minn = min(size(a,2),size(b,2));
for i = 1:length(lags)
    subabig = abig(:,i:i+minn-1);
    temp = subabig.*b(:,1:minn);
    cc(i) = sum(temp(:))/sqrt(sum(a(:).^2)*sum(b(:).^2));
end