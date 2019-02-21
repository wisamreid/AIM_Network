function out = runStoi(a,b,fs_a,fs_b,warp)
%out = RUNSTOI(a,b,fs)
% warp = option to timewarp.
%matches the lenths of the input vectors before calling stoi()
%zeropadding is applied to the end of the shorter vector

%check row or vector
if ~isrow(a)
    a = a';
end
if ~isrow(b)
    b = b';
end

%timewarp
if ~exist('warp','var'), warp = 0; end
if warp; b = timeMatch(a,b,fs_a,fs_b,0,0); end

%zero pad if necessary
d = length(a) - length(b);
if d > 0
    b = [b zeros(1,d)];
elseif d < 0
    a = [a zeros(1,abs(d))];
end

out = stoi(a,b,fs_a,fs_b);