function RDMs=concatRDMs_unwrapped(varargin)
% concatenates dissimilarity sets (ltv or square RDMs). All inputs should
% have the same number of dissimilarity entries.
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

RDMs=[];
for RDMvarI=1:nargin
    RDMs=cat(3,RDMs,unwrap(varargin{RDMvarI}));
end