function p=relRankIn_includeValue_lowerBound(set,value)

% returns the relative rank of value within set.
% the relative rank is the proportion of set smaller than value.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council
set=[set(:);value];
p=sum(set(:)<value)/numel(set);
