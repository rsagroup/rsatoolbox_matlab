function Xscaled=scale01(X,range)
% USAGE
%       Xscaled=scale01(X[,range])
%
% FUNCTION
%       ..linearly scales matrix X into the range [0,1].
%
%       if range is given, min(range) and max(range) define the scaling and
%       shifting, rather than min(X(:)) and max(X(:)).
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

if ~exist('range','var')
    mn=min(X(:)); mx=max(X(:));
else
    mn=min(range(:)); mx=max(range(:));
end    
Xscaled=(X-mn)./(mx-mn);

end%function
