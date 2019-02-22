function out = regroup(in, groupSizes)
% out = regroup(in, groupSizes)
%   groupSizes - a vector of various group lengths. Should have more than
%                one element. For example, to group the data first by
%                frequency, then locations (i.e. [1:36,1:36,1:36,...]), then
%                groupSizes should be [36 5].
% @ Kenny Chou
% Boston Univeristy 2/2019
temp = in;
out = zeros(size(in));
group1 = groupSizes(1); %36;
group2 = groupSizes(2); %5;
if size(in,2) >= size(in,1)
    for i = 1:group2
        idxIn = (1:group2:group1*group2)+i-1;
        idxOut = 1+(group1)*(i-1):group1*i;
        temp(:,idxOut) = in(:,idxIn);
    end
end
if size(in,1) >= size(in,2)
    for i = 1:group2
        idxIn = (1:group2:group1*group2)+i-1;
        idxOut = 1+(group1)*(i-1):group1*i;
        out(idxOut,:) = temp(idxIn,:);
    end
end
% prefFreq = repmat((1:36),1,5)';
% prefLocs = [ones(1,36) ones(1,36)*2 ones(1,36)*3 ones(1,36)*4 ones(1,36)*5]';
% figure;
% subplot(131); imagesc(in); title('initial connectivity matrix')
% subplot(132); imagesc(temp); title('intermediate')
% subplot(133); imagesc(out); title('regrouped connectivity matrix')