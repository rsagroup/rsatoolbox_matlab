function v = randomlyPermute(v)
% randomly permutes the elements of the input vector (v).
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

v = v(randomPermutation(numel(v)));

end%function
