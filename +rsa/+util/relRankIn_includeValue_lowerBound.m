function p=relRankIn_includeValue_lowerBound(set,value)

% returns the relative rank of value within set.
% the relative rank is the proportion of set smaller than value.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

set=[set(:);value];
p=sum(set(:)<value)/numel(set);

end%function
