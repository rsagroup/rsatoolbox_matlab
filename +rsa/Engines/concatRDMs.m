function RDMs=concatRDMs(varargin)
% concatenates representational dissimilarity matrices and returns the
% concatenated RDMs in a wrapped structure.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council



for RDMI=1:nargin
    
    if ~isstruct(varargin{RDMI})
        varargin{RDMI}=wrapRDMs(varargin{RDMI});
    end
    
    RDMs(RDMI)=varargin{RDMI};
end