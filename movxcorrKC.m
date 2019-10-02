function cc = movxcorrKC(a,b,framelen)
% compute normalized cross-correlation between a and b by taking a windowed
% frame with length framelen for 
i = 1;
while i+framelen < size(a,2)
    smalla = a(:,i:i+framelen-1);
    smallb = b(:,i:i+framelen-1);
    cc(i) = xcorrKC(smalla,smallb);
    i = i+1;
end
