function Xscaled=scale01(X,range)
% USAGE
%       Xscaled=scale01(X[,range])
%
% FUNCTION
%       ..linearly scales matrix X into the range [0,1].
%
%       if range is given, min(range) and max(range) define the scaling and
%       shifting, rather than min(X(:)) and max(X(:)).
%
%       if all entries in the matrix are the same, the result is an
%       all-zeros matrix of the same size.
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
    % Get the top and bottom of the scale from X
    mini = min(X(:));
    maxi = max(X(:));
else
    % Get the top and bottom of the scale from 'range'
    mini = min(range(:));
    maxi = max(range(:));
end%if

% If the matrix is constant (which is possible), or the specified range is
% constant (which shouldn't be allowed), then we can't scale by the width
% of the range, so we put every entry at 0
if maxi == mini
    Xscaled = zeros(size(X));
else
    Xscaled = (X - mini) ./ (maxi - mini);
end%if

end%function
