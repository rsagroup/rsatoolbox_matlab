function RDMsOut = concatenateRDMs(varargin)
% RDMsOut = concatenateRDMs(RDMs1[, RDMs2[, ...]])
% concatenates a number of wrapped RDMs (RDM structures). 
% the main difference between this function and the similar one called
% 'concatRDMs' is that this function requires the input to be wrapped RDM
% st
% CW 5-2010, 6-2010
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

	if nargin == 1
		RDMsIn = varargin{1};
		RDMsOut = interleaveRDMs(RDMsIn);
	else
		arg1 = catTwoRDMs(varargin{1}, varargin{2});
		arg1 = {arg1};
		arg2 = varargin(3:end);
		args = [arg1 arg2];
		RDMsOut = concatenateRDMs(args{:});
	end

end%function

%% Subfunctions:

function RDMs = catTwoRDMs(RDMs1, RDMs2)
	import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*
	RDMs1 = interleaveRDMs(RDMs1, true);
	RDMs2 = interleaveRDMs(RDMs2, true);
	RDMs2 = orderfields(RDMs2, fieldnames(RDMs1));
	cRDMs1 = struct2cell(RDMs1);
	cRDMs2 = struct2cell(RDMs2);
	cRDMs = [cRDMs1 cRDMs2];
	RDMs = cell2struct(cRDMs, fieldnames(RDMs1), 1);
end%function
