function RDM=squareRDM(RDM)
% converts RDM to square form
%
% if they are in square form already they are vectorized
% using the lower triangle (matlab squareform convention) and resquared,
% thus fixing inconsistencies between upper and lower triangle (the lower
% is used).
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

RDM=squareform(vectorizeRDM(RDM));

end%function
