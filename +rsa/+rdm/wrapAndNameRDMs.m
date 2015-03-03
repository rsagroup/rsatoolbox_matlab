function RDMs=wrapAndNameRDMs(RDMs,names)

% wraps the RDMs in argument RDMs into a structured array and assigns the
% names in arguments names.
%__________________________________________________________________________
% Copyright (C) 2012 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

RDMs=wrapRDMs(RDMs);

for nameI=1:numel(names)
    RDMs(nameI).name=names{nameI};
end

end%function
