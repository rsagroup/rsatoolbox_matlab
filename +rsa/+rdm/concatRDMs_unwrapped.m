function RDMs=concatRDMs_unwrapped(varargin)
% concatenates dissimilarity sets (ltv or square RDMs). All inputs should
% have the same number of dissimilarity entries.
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

RDMs=[];
for RDMvarI=1:nargin
    RDMs=cat(3,RDMs,unwrap(varargin{RDMvarI}));
end

end%function
