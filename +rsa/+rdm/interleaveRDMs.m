function RDMs = interleaveRDMs(RDMs, varargin)
%
%  RDMs = interleaveRDMs(RDMs[, transpose?])
%
%  interleaveRDMs is a function which accepts an arbitrarily-high-dimensional
%  struct of RDMs and recursively "interleaves" the final (non-singleton)
%  dimensions out by shuffling and stacking, finally returning a 1 (or 2, with 1
%  singleton dimension) dimensional struct of RDMs all ready for display.
%  
%  Cai Wingfield 11-2009
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

if numel(varargin) == 0
	transposed = false;
else
	transposed = varargin{1};
end

if transposed
	RDMs = reshape(RDMs, [], 1);
else
	RDMs = reshape(RDMs, 1, []); 
end

end%function
