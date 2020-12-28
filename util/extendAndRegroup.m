function out = extendAndRegroup(netcon,nFreqs,nLocs)
% Takes in a NETCON for a small network (e.g., one frequency channel,
% channels are independent), repeats the netcon along another dimension
% (e.g., locations), and regroups the netcon so that connectivities are
% grouped by the 2nd dimension.
%
% @ Kenny Chou
% Boston Univeristy 09/2020

netconPerFreq = netcon;
netconCell = repmat({netconPerFreq},1,nFreqs);
netconGroupBy1 = blkdiag(netconCell{:});
netconGroupBy2 = regroup(netconGroupBy1, [nFreqs,nLocs]);
out = netconGroupBy2;