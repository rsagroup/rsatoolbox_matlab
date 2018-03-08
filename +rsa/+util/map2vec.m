function vec=map2vec(map,mask)
% USAGE
%       vec=map2vec(map,mask)
%
% FUNCTION
%       to convert a spatial map of dimensions of array mask into space
%       column vector vec. mask is a logical array whose true values will
%       correspond to vec in number and order.
%
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

mask=logical(mask);
vec=nan(sum(mask(:)),1);
vec=map(mask);

end%function
