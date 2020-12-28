function cc = movxcorrKC(a,b,framelen)
% running cross correlation between matrices a and b
% framelen - window length to compare cross correlation of

%pad to match length
if size(a,2) < size(b,2)
    a = [a zeros(size(a,1),size(b,2)-size(a,2))];
elseif size(a,2) > size(b,2)
    b = [b zeros(size(b,1),size(a,2)-size(b,2))];
end
        
i = 1;
skip = 50; % samples to skip (kinda like n-overlap)
while i+framelen < size(a,2)
    smalla = a(:,i:i+framelen-1);
    smallb = b(:,i:i+framelen-1);
    if sum(smalla(:)) == 0 || sum(smallb(:)) == 0
        cc(i) = 0;
    else
        temp = (smalla.*smallb);
        cc(i) = sum(temp(:))/sqrt(sum(smalla(:).^2)*sum(smallb(:).^2)+10E-10);
    end
    i = i+skip;
end
