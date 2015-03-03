function stringOut = removeSpaces(stringIn);
%
% removeSpaces will remove all spaces from a string
%
% Cai Wingfield 1-2010
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

if ~isempty(stringIn)
	stringOut = strrep(stringIn, ' ', '');
end

end%function
