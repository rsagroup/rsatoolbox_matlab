function RDMs=concatRDMs(varargin)
% concatenates representational dissimilarity matrices and returns the
% concatenated RDMs in a wrapped structure.
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

for RDMI=1:nargin
    
    if ~isstruct(varargin{RDMI})
        varargin{RDMI}=wrapRDMs(varargin{RDMI});
    end
    
    RDMs(RDMI)=varargin{RDMI};
end

end%function