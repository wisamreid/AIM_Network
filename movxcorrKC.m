function cc = movxcorrKC(a,b,framelen)

i = 1;
while i+framelen < size(a,2)
    smalla = a(:,i:i+framelen-1);
    smallb = b(:,i:i+framelen-1);
    cc(i) = xcorrKC(smalla,smallb);
    i = i+1;
end
