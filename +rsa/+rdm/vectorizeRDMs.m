function RDMs_utv=vectorizeRDMs(RDMs)
% converts set of RDMs (stacked along the 3rd dimension)
% to lower-triangular form (set of row vectors)
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

if isstruct(RDMs)
    % wrapped
    RDMs_struct=RDMs;
    RDMs=unwrapRDMs(RDMs_struct);
    
    nRDMs=size(RDMs,3);
    RDMs_utv=[];
    for RDMI=1:nRDMs
        RDMs_utv=cat(3,RDMs_utv,vectorizeRDM(RDMs(:,:,RDMI)));
    end
    
    RDMs_utv=wrapRDMs(RDMs_utv,RDMs_struct);
else
    % bare
    nRDMs=size(RDMs,3);
    RDMs_utv=[];
    for RDMI=1:nRDMs
        RDMs_utv=cat(3,RDMs_utv,vectorizeRDM(RDMs(:,:,RDMI)));
    end
end

end%function
