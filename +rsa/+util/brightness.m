function b=brightness(RGBrows)
% given triples of RGB values, returns the overall brightness vector. This
% is performed on each row of the input (RGBrows) independently.
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

RGBweights=[.241 .691 .068]';

b=sqrt(RGBrows*RGBweights);

end%function
