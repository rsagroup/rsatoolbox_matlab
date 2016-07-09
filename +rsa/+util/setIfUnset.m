function options=setIfUnset(options,field,value)
% if options.(field) is empty or doesn't exist, this function sets options.(field) to value.
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

if ~isfield(options, field) || isempty(options.(field))
       options.(field)=value;
end

end%function
