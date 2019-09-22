function cc = movxcorrKC(a,b,framelen)

%pad to match length
if size(a,2) < size(b,2)
    a = [a zeros(size(a,1),size(b,2)-size(a,2))];
elseif size(a,2) > size(b,2)
    b = [b zeros(size(b,1),size(a,2)-size(b,2))];
end
        
i = 1;
while i+framelen < size(a,2)
    smalla = a(:,i:i+framelen-1);
    smallb = b(:,i:i+framelen-1);
    if sum(smalla(:)) == 0 || sum(smallb(:)) == 0
        cc(i) = 0;
    else
        cc(i) = xcorrKC(smalla,smallb);
    end
    i = i+1;
end
