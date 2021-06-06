function out = runStoi(a,b,fs_a,fs_b)
% out = RUNSTOI(a,b,fs)
%   Calculates STOI on two soundwaves
%   matches the lengths & normalize the input vectors before calling stoi()
%   zeropadding is applied to the end of the shorter vector.
% Inputs:
%   a, b = soundwaves. If dual channel, combine the two.
%
% STOI operates on single channel audio inputs.
%
% @Kenny Chou, BU Hearing Reserach Center
% 20190916 - removed outdated 'warp' option

if ~exist('fs_b','var'), fs_b = fs_a; end

%combine channels
if ~isvector(a)
    if size(a,1) > size(a,2)
        a = sum(a,2); 
    else
        a = sum(a);
    end
end
if ~isvector(b)
    if size(b,1) > size(b,2)
        b = sum(b,2); 
    else
        b = sum(b);
    end
end

%check row or vector
if ~isrow(a)
    a = a';
end
if ~isrow(b)
    b = b';
end

%normalize inputs
a = a/max(abs(a));
b = b/max(abs(b));

%zero pad if necessary
d = length(a) - length(b);
if d > 0
    b = [b zeros(1,d)];
elseif d < 0
    a = [a zeros(1,abs(d))];
end

out = stoi(a,b,fs_a,fs_b);